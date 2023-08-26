#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
#![deny(unreachable_pub)]
#![warn(rust_2018_idioms)]

//! This crate provides fast DNA sequence extraction from 2bit files, a
//! [standard format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) in bioinformatics.
//!
//! The motivation for this crate is speed.
//! Extracting sequences is 1.5-30x faster than the best alternative, depending on use case.
//!
//! Tested on Linux, Mac, and Windows:
//! <a href="https://dl.circleci.com/status-badge/redirect/gh/andrewdelong/twobitreader-rust/tree/main">
//!   <img style="vertical-align:middle; margin-top:-2px;"
//!        src="https://dl.circleci.com/status-badge/img/gh/andrewdelong/twobitreader-rust/tree/main.svg?style=shield&circle-token=866d66445adcd45b6a135f83a6211987fa1c4cf3"/>
//! </a>
//!
//! # Examples
//!
//! **Extracting sequences** is straightforward with [`TwobitReader`]:
//! ```no_run
//! # use std::io;
//! # use twobitreader::TwobitReader;
//! let tbr = TwobitReader::open("hg38.2bit")?;  // Human genome, build 38
//! let seq = tbr.get("chr1", 10000, 10005);     // String ("TAACC")
//! # Ok::<(), io::Error>(())
//! ```
//!
//! **Concatenation** works by iterating over (start, end) pairs.
//! For example, assembling a spliced transcript:
//! ```no_run
//! # use std::io;
//! # use std::iter::zip;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! // Exon ranges for human FKHL6 gene transcript (Gencode v43)
//! let exons = [(1389575, 1391118),             // Exon 1 (start, end)
//!              (1394695, 1395603)];            // Exon 2 (start, end)
//! let transcript = tbr.concat("chr6", exons);  // String
//! # //
//! # // Same, but with zipped parallel arrays
//! # let starts = [1389575, 1394695];
//! # let ends = [1391118, 1395603];
//! # let transcript = tbr.concat("chr6", zip(starts, ends));
//! # Ok::<(), io::Error>(())
//! ```
//!
//! **Parallelism** is easy with crates like [`rayon`](https://docs.rs/rayon/latest/rayon/).
//! For example, batch extraction of sequences:
//! ```no_run
//! # use std::io;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! use rayon::prelude::*;
//! let args = [("chr1", 10000, 15000),
//!             ("chr1", 30000, 35000), /* ... */ ];
//! let seqs = args.into_par_iter()
//!     .map(|(chrom, start, end)| tbr.get(chrom, start, end))
//!     .collect::<Vec<_>>();  // Vec<String>
//! # Ok::<(), io::Error>(())
//! ```
//! Or, assembling a batch of spliced transcripts in parallel:
//! ```no_run
//! # use std::io;
//! # use std::collections::HashMap;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! use rayon::prelude::*;
//! let transcripts = [                           // (transcript_id, chromosome, exons)
//!     ("ENST00000407983.7", "chr2", vec![(264899, 265007),     // Exon 1 (start, end)
//!                                        (271865, 271939),     // Exon 2 (start, end)
//!                                        (272036, 272557)]),   // Exon 3 (start, end)
//!     ("ENST00000319331.4", "chr3", vec![(3799430, 3799919),   // Exon 1 (start, end)
//!                                        (3844363, 3849834)]), // Exon 2 (start, end)
//!     /* ... */
//! ];
//! let seqs = transcripts.into_par_iter()
//!     .map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons)))
//!     .collect::<HashMap<_, _>>();       // HashMap<&str, String>
//! let seq = &seqs["ENST00000407983.7"];  // Look up transcript sequence
//! # Ok::<(), io::Error>(())
//! ```
//!
//! # Speed
//!
//! Two tasks were benchmarked:
//! - **exons**: extract 133,388 distinct human exon sequences;
//! - **transcripts**: concatenate 319,468 exons into 29,180 human spliced transcript sequences.
//!
//! Speed depends on parallelism and page cache (hot vs cold):
//! - **hot** represents computations that are repeated, or are run interactively, on the same server;
//! - **cold** represents a first run of a computational pipeline, and is bound by disk speed.
//!
//! The table below shows running times in milliseconds. Experimental details are `doc/comparison-details.md`.
#![doc = include_str!("../doc/comparison-result.html")]
//!
//! # Dependencies
//!
//! * `byteorder` for handling endian-ness
//! * `memmap2` for memory mapping the 2bit file
//! * `seq-macro` for generating 2bit decoder lookup table

// Standard library
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Cursor, ErrorKind::InvalidData, ErrorKind::UnexpectedEof, Read};
use std::iter::{repeat_with, zip};
use std::mem::size_of;
use std::panic;
use std::path::Path;
use std::vec;

// Crate modules
mod decode;
use decode::{decode, NUCS_PER_U8};
mod parallel;
use parallel::try_parallel_for;

// Dependencies
use byteorder::{BigEndian, ByteOrder, LittleEndian, ReadBytesExt};
use memmap2::{Mmap, MmapOptions};

// Big-endian 2bit files are supported by this crate, but big-endian compile targets are not.
#[cfg(target_endian = "big")]
compile_error!("twobitreader is not yet implemented for big-endian targets.");

/// A reader for a single [2bit file](http://genome.ucsc.edu/FAQ/FAQformat.html#format7).
///
#[derive(Debug)]
pub struct TwobitReader {
    mmap: Mmap,                          // Memory map of entire 2bit file.
    seqs: Vec<TwobitSequence>,           // Sequence data, in same order as in the file.
    seq_by_name: HashMap<String, usize>, // Lookup sequence index (seqs[index]) by name.
}

/// A reader for a specific sequence record within the 2bit file.
///
#[derive(Debug)]
struct TwobitSequence {
    dna_offset: usize, // Offset (within mmap) to first byte of packed DNA
    dna_bytes: usize,  // Number of bytes (not nucleotides!) of packed DNA
    dna_len: usize,    // Number of nucleotides (not bytes!) in DNA sequence
    nblocks: Blocks,   // N block starts and ends [start0, end0, start1, end1, ...]
    masks: Blocks,     // Mask block starts and ends
    name: String,      // Sequence name ("chr3", etc.)
}

// Sorted parallel arrays where (starts[i], ends[i]) is the range of one N-block or
// lowercase-block in a sequence record.
#[derive(Debug)]
struct Blocks {
    starts: Vec<u32>,
    ends: Vec<u32>,
}

// Mark TwobitSequence as safe to send to a thread and safe to and share between threads.
// The compiler cannot assume that dna_ptr references immutable memory, so it will not allow
// allow a thread to borrow a TwobitSequence unless we explicitly indicate that it is safe.
unsafe impl Sync for TwobitSequence {}
unsafe impl Send for TwobitSequence {}

impl TwobitReader {
    /// Opens a 2bit file for reading.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;  // Human genome, build 38
    /// let seq = tbr.get("chr2", 10000, 10010);     // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::open_inner(path, false)
    }

    /// Opens a 2bit file for reading, with lowercase masks applied.
    ///
    /// The lowercase mask feature of 2bit files is mainly relevant for
    /// sequence search, for example with the
    /// [BLAT suite](https://genome.ucsc.edu/goldenpath/help/blatSpec.html) of tools.
    /// Specifically, lowercase letters typically indicate a region that should be ignored
    /// (not searched) during a sequence search.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open_masked("hg38.2bit")?;  // Human genome, build 38
    /// let seq = tbr.get("chr2", 10000, 10010);            // "CGTATcccac" (mixed case)
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn open_masked<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::open_inner(path, true)
    }

    // Implements opening the file and initializing a TwobitReader instance for it.
    fn open_inner<P: AsRef<Path>>(path: P, use_mask: bool) -> io::Result<Self> {
        // Map the file and read the header and mask arrays in the order they appear in the file
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let seqs = read_seqs(&mmap, use_mask)?;

        // Allow fast lookup of sequences by name
        let seq_by_name = HashMap::from_iter(seqs.iter().enumerate().map(|(i, seq)| (seq.name.clone(), i)));
        if seq_by_name.len() != seqs.len() {
            return Err(io::Error::new(InvalidData, "duplicate sequence name detected."));
        }

        Ok(Self { mmap, seqs, seq_by_name })
    }

    /// Returns the number of sequence records in the file.
    pub fn num_seqs(&self) -> usize {
        self.seqs.len()
    }

    /// Iterates over the sequence record names, in the order they appear in the file.
    pub fn iter_names(&self) -> impl Iterator<Item = &str> + '_ {
        self.seqs.iter().map(|seq| seq.name.as_str())
    }

    /// Returns a vector of sequence record names, in the order they appear in the file.
    pub fn names(&self) -> Vec<&str> {
        self.iter_names().collect()
    }

    /// Returns whether `name` matches a sequence record in the file.
    pub fn contains_name<N: AsRef<str>>(&self, name: N) -> bool {
        self.seq_by_name.contains_key(name.as_ref())
    }

    /// Returns the length of DNA for the named sequence record.
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found.
    ///
    pub fn seq_len<N: AsRef<str>>(&self, name: N) -> usize {
        self.get_seq(name).dna_len
    }

    /// Extracts range `start..end` (0-based, exclusive end) from the named sequence record.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let seq = tbr.get("chr2", 10000, 10010);     // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or the range was invalid.
    ///
    #[inline(never)]
    pub fn get<N: AsRef<str>>(&self, name: N, start: usize, end: usize) -> String {
        let mut dst = String::new();
        self.get_into(name, start, end, &mut dst);
        dst
    }

    /// A version of [`get`](Self::get) where the result is stored in `dst`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let mut dst = String::new();
    /// tbr.get_into("chr2", 10000, 10010, &mut dst);  // dst = "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or the range was invalid.
    ///
    pub fn get_into<N: AsRef<str>>(&self, name: N, start: usize, end: usize, dst: &mut String) {
        // Check range before trying to do any might-panic arithmetic with (start, end)
        let seq = self.get_seq(name);
        check_range(seq, start, end);

        let len = end - start;
        if len == 0 {
            // Empty range, so clear the output string.
            dst.clear();
        } else {
            // Non-empty range, so size dst and decode into its u8 buffer.
            // This is safe because the decoded ASCII is already valid utf8.
            if dst.capacity() < len {
                dst.reserve_exact(len - dst.len());
                assert!(dst.capacity() == len);
            }
            let dst = unsafe { dst.as_mut_vec() };
            unsafe { dst.set_len(len) };
            self.get_into_u8(seq, start, dst);
        }
    }

    /// A version of [`get`](Self::get) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let seq = tbr.get("chr2", 10001, 10010);  // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or the range was invalid.
    ///
    pub fn get_inclusive<N: AsRef<str>>(&self, name: N, start: usize, end: usize) -> String {
        // Check valid start and then convert range to 0-based exclusive.
        check_start_inclusive(start);
        self.get(name, start - 1, end)
    }

    /// A version of [`get_into`](Self::get_into) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let mut dst = String::new();
    /// tbr.get_inclusive_into("chr2", 10001, 10010, &mut dst);  // dst = "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or the range was invalid.
    ///
    pub fn get_inclusive_into<N: AsRef<str>>(&self, name: N, start: usize, end: usize, dst: &mut String) {
        // Check valid start and then convert range to 0-based exclusive.
        check_start_inclusive(start);
        self.get_into(name, start - 1, end, dst);
    }

    /// Concatenates a batch of sequence ranges into a single string.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let ranges = [(10000, 10006), (10006, 10010)];  // "CGTATC" "CCAC"
    /// let seq = tbr.concat("chr2", ranges);           // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or if any range was invalid.
    ///
    pub fn concat<N, R>(&self, name: N, ranges: R) -> String
    where
        N: AsRef<str>,
        R: IntoIterator<Item = (usize, usize)>,
        R::IntoIter: Clone,
    {
        // Compute total length by summing individual interval lengths. Also check the ranges here,
        // since get_into_u8 does not produce user-friendly panic messages for invalid ranges.
        let seq = self.get_seq(name);
        let ranges = ranges.into_iter();
        let total_len = ranges
            .clone()
            .map(|(start, end)| {
                check_range(seq, start, end);
                end - start
            })
            .sum();
        let mut dst = String::new();
        if total_len > 0 {
            // Pre-size dst to the length needed. Safe because all decoded nucleotides are valid utf8.
            let dst_vec = unsafe { dst.as_mut_vec() };
            dst_vec.reserve_exact(total_len);
            unsafe { dst_vec.set_len(total_len) };

            // Decode each interval into its respective slice of dst.
            let mut offset = 0;
            for (start, end) in ranges {
                let len = end - start;
                self.get_into_u8(seq, start, &mut dst_vec[offset..offset + len]);
                offset += len;
            }
        }
        dst
    }

    /// A version of [`concat`](Self::concat) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let ranges = [(10001, 10006), (10007, 10010)];  // "CGTATC" "CCAC"
    /// let seq = tbr.concat("chr2", ranges);           // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or if any range was invalid.
    ///
    pub fn concat_inclusive<N, R>(&self, name: N, ranges: R) -> String
    where
        N: AsRef<str>,
        R: IntoIterator<Item = (usize, usize)>,
        R::IntoIter: Clone,
    {
        // Check each start, and convert range to 0-based exclusive.
        self.concat(
            name,
            ranges.into_iter().map(|(start, end)| {
                check_start_inclusive(start);
                (start - 1, end)
            }),
        )
    }

    // Returns a reference to the [`TwobitSequence`] for the named sequence record.
    // Panics if no sequence record has that name.
    fn get_seq<N: AsRef<str>>(&self, name: N) -> &TwobitSequence {
        let index = *self.seq_by_name.get(name.as_ref()).expect("sequence name not found");
        &self.seqs[index]
    }

    // Decodes sequence [start..end] into dst, where end = start + dst.len().
    // All bytes in dst are replaced. The first byte in slice is assumed
    fn get_into_u8(&self, seq: &TwobitSequence, start: usize, dst: &mut [u8]) {
        // Slice spanning all packed 2bit DNA data for this sequence record.
        let dna = &self.mmap[seq.dna_offset..seq.dna_offset + seq.dna_bytes];

        // Decode packed 2-bit from dna into dst from the given start position.
        decode(start, dna, dst);

        // Apply N-block and lowercase masks to dst in-place.
        search_blocks(&seq.nblocks, start, dst, nblock_fill);
        search_blocks(&seq.masks, start, dst, mask_fill);
    }
}

// Panics if (start, end) is an invalid 0-based exclusive range for this TwobitSequence.
// (The error message does not format the actual start/end values that were passed in.
//  This is because the 0-based get functions are used to implement the 1-based get
//  functions, and so the start/end seen here may not be the actual ones users provided,
//  which would result in more confusion for the user than simply omitting the numbers.)
#[inline]
fn check_range(seq: &TwobitSequence, start: usize, end: usize) {
    if start > end {
        panic!("invalid range (start > end)");
    }
    if end > seq.dna_len {
        panic!("invalid end (end > dna_len)");
    }
}

// Panics ifs start is not valid for an inclusive range
#[inline]
fn check_start_inclusive(start: usize) {
    if start == 0 {
        panic!("invalid start (start = 0)");
    }
}

// Reads the header and sequence data from a .2bit file.
// The resulting vec of TwobitSequences is intended to be assigned to the TwobitReader::seqs field.
// Block intervals are accessed (paged in from disk), but the sequence data itself is not paged in.
fn read_seqs(mmap: &Mmap, use_mask: bool) -> io::Result<Vec<TwobitSequence>> {
    // Twobit file signature as little-endian (native) or big-endian (non-native)
    const SIGNATURE_LILEND: u32 = 0x1A412743;
    const SIGNATURE_BIGEND: u32 = 0x4327411A;

    // Start reading the file, using appropriate byte order
    let signature = Cursor::new(mmap).read_u32::<LittleEndian>()?;
    match signature {
        SIGNATURE_LILEND => read_seqs_from_endian::<LittleEndian>(mmap, use_mask),
        SIGNATURE_BIGEND => read_seqs_from_endian::<BigEndian>(mmap, use_mask),
        _ => Err(io::Error::new(InvalidData, "invalid file signature.")),
    }
}

// Implements read_seqs after the endian-ness of the file has been determined from the signature.
fn read_seqs_from_endian<B: ByteOrder>(mmap: &Mmap, use_mask: bool) -> io::Result<Vec<TwobitSequence>> {
    // Start cursor immediately after the signature field.
    let mut cursor = Cursor::new(mmap);
    cursor.set_position(size_of::<u32>() as u64);

    // Check the file version
    let version = cursor.read_u32::<B>()?;
    if version != 0 {
        return Err(io::Error::new(InvalidData, "file version not recognized."));
    }

    // Read number of sequences and read past the reserved dword.
    let num_seqs = cursor.read_u32::<B>()? as usize;
    let _reserved = cursor.read_u32::<B>()?;

    // Read each sequence header (name, data_offset).
    let mut seq_headers = Vec::with_capacity(num_seqs);
    for _ in 0..num_seqs {
        // Read the sequence name, first into a byte vec.
        let name_size = cursor.read_u8()? as usize;
        let mut name_bytes = vec![0u8; name_size];
        cursor.read_exact(&mut name_bytes)?;

        // Convert the name bytes into a String object
        let name = match String::from_utf8(name_bytes) {
            Ok(name) => name,
            Err(_) => return Err(io::Error::new(InvalidData, "failed reading sequence name as utf8.")),
        };

        // Get the offset of the data arrays for this sequence
        let data_offset = cursor.read_u32::<B>()? as usize;
        let curr_offset = cursor.position() as usize;
        if data_offset <= curr_offset {
            return Err(io::Error::new(InvalidData, "invalid data offset; file may be malformed."));
        }

        seq_headers.push((name, data_offset));
    }

    // Create a pre-sized Vec<Option<TwobitSequence>, initially filled with None.
    // (Use from_iter: cannot use vec![] or resize because TwobitSequence cannot clone.)
    let mut seqs = Vec::from_iter(repeat_with(|| None).take(num_seqs));

    // Read each sequence's data block, in parallel. Reading in parallel is very important
    // for the speed of opening a "cold" file that is not yet in the page cache.
    // seq_headers[i] is used to initialize TwobitSequence, which is then assigned
    // to seqs[i], replacing the initial None value.
    try_parallel_for(zip(seq_headers, &mut seqs), |((name, data_offset), seq)| -> io::Result<()> {
        let mut cursor = Cursor::new(mmap);
        cursor.set_position(data_offset as u64);

        // Read and process the block arrays, but the packed DNA can just stay mapped as-is.
        let dna_len = cursor.read_u32::<B>()? as usize;
        let nblocks = read_blocks::<B>(&mut cursor, true)?;
        let masks = read_blocks::<B>(&mut cursor, use_mask)?;
        let _reserved = cursor.read_u32::<B>()?;

        // Check blocks are valid.
        check_blocks(&nblocks, dna_len)?;
        check_blocks(&masks, dna_len)?;

        // Get offset and length of packed DNA, and make sure it's within the file's size.
        let dna_offset = cursor.position() as usize;
        let dna_bytes = (dna_len + NUCS_PER_U8 - 1) / NUCS_PER_U8;
        if dna_offset + dna_bytes > cursor.get_ref().len() {
            return Err(io::Error::new(UnexpectedEof, "dna bytes were truncated in file."));
        }

        // Initialize the sequence.
        *seq = Some(TwobitSequence { dna_offset, dna_bytes, dna_len, nblocks, masks, name });
        Ok(())
    })?;

    // If we got this far, all seqs entries are Some(_), so return as unwrapped.
    Ok(seqs.into_iter().map(|x| x.unwrap()).collect())
}

// Check if items are sorted. Pasted from nightly slice::is_sorted until in stable rust.
#[cfg(debug_assertions)]
fn is_sorted<T>(items: &[T]) -> bool
where
    T: PartialOrd<T>,
{
    items.is_empty() || items.windows(2).all(|x| x[0] < x[1])
}

// Reads a Twobit block section from the given cursor, returning a Blocks struct.
// If used=false, then the cursor is simply advanced, and an empty Blocks struct is returned.
fn read_blocks<B: ByteOrder>(cursor: &mut Cursor<&Mmap>, used: bool) -> io::Result<Blocks> {
    let num_blocks = cursor.read_u32::<B>()? as usize;
    let mut starts = Vec::new();
    let mut ends = Vec::new();

    if used {
        // Read starts.
        starts.resize(num_blocks, 0);
        cursor.read_u32_into::<B>(&mut starts)?;

        // Read sizes and convert them to ends.
        ends.resize(num_blocks, 0);
        cursor.read_u32_into::<B>(&mut ends)?;
        for (start, end) in zip(&starts, &mut ends) {
            *end += *start;
        }
    } else {
        // Advance the cursor, without reading the data itself.
        let new_offset = cursor.position() as usize + 2 * num_blocks * size_of::<u32>();
        if new_offset > cursor.get_ref().len() {
            return Err(io::Error::new(UnexpectedEof, "block indices truncated; file may be malformed."));
        }
        cursor.set_position(new_offset as u64);
    }

    Ok(Blocks { starts, ends })
}

#[cfg(debug_assertions)]
fn check_blocks(blocks: &Blocks, dna_len: usize) -> io::Result<()> {
    if !is_sorted(&blocks.starts) || !is_sorted(&blocks.ends) {
        return Err(io::Error::new(InvalidData, "block indices not sorted; file may be malformed."));
    }

    for (start, end) in zip(&blocks.starts, &blocks.ends) {
        if start >= end {
            return Err(io::Error::new(InvalidData, "block start >= end; file may be malformed."));
        }
        if *start as usize > dna_len {
            return Err(io::Error::new(InvalidData, "invalid block start; file may be malformed."));
        }
        if *end as usize > dna_len {
            return Err(io::Error::new(InvalidData, "invalid block end; file may be malformed."));
        }
    }

    Ok(())
}

#[cfg(not(debug_assertions))]
fn check_blocks(_blocks: &Blocks, _dna_len: usize) -> io::Result<()> {
    Ok(())
}

// Searches for block ranges that overlap query range [start..end] where end = start+dst.len().
// For each block range found, calls f() with the corresponding sub-slice of dst,
// so that its contents may be modified (e.g., replaced with N, or replaced with lowercase).
fn search_blocks<F>(blocks: &Blocks, start: usize, dst: &mut [u8], f: F)
where
    F: Fn(&mut [u8]),
{
    // Find index of first block having block_end > start.
    let i = blocks.ends.partition_point(|&x| x <= start as u32);

    // For each consecutive (block_start, block_end) pair, clamp to the query
    // (start, end) range and call f() to mutate the excerpted u8 slice.
    // Stop early if the clamped range is empty.
    let end = start + dst.len();
    for (block_start, block_end) in zip(&blocks.starts[i..], &blocks.ends[i..]) {
        let block_start = (*block_start as usize).max(start);
        let block_end = (*block_end as usize).min(end);
        if block_start >= block_end {
            break;
        }
        f(&mut dst[block_start - start..block_end - start]);
    }
}

// Replaces all bytes in dst with the letter 'N'
#[inline]
fn nblock_fill(dna: &mut [u8]) {
    dna.fill(b'N');
}

// Replaces all uppercase letters in dst with lowercase.
// The contents of dna are assumed to be uppercase letters already.
#[inline]
fn mask_fill(dna: &mut [u8]) {
    for nuc in dna {
        debug_assert!(*nuc >= b'A' && *nuc <= b'Z');
        *nuc += b'a' - b'A';
    }
}
