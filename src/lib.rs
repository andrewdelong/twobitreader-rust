#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
#![deny(unreachable_pub)]
#![warn(rust_2018_idioms)]
//!
//! This crate provides fast DNA sequence extraction from 2bit files, a
//! [standard format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) in bioinformatics.
//!
//! The motivation for this crate is speed.
//! Extracting sequences is consistently faster than the best alternative.
//! The focus is raw reading from 2bit, but fast concatenation and reverse-complement methods
//! are also provided to make higher-level use cases easier.
//!
//! [![CI](https://github.com/andrewdelong/twobitreader-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/andrewdelong/twobitreader-rust/actions/workflows/ci.yml)
//!
//! # Examples
//!
//! **Extracting sequences** is straightforward:
//! ```no_run
//! # use std::io;
//! # use twobitreader::TwobitReader;
//! let tbr = TwobitReader::open("hg38.2bit")?; // Human genome, build 38
//! let seq = tbr.get("chr1", 10000, 10005);    // -> String ("TAACC")
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
//! let transcript = tbr.concat("chr6", &exons); // -> String
//! # //
//! # // Same, but with zipped parallel arrays
//! # let starts = [1389575, 1394695];
//! # let ends = [1391118, 1395603];
//! # let transcript = tbr.concat_iter("chr6", zip(starts, ends));
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
//!     .collect::<Vec<_>>(); // -> Vec<String>
//! # Ok::<(), io::Error>(())
//! ```
//! Or, assembling a batch of spliced transcripts in parallel:
//! ```no_run
//! # use std::io;
//! # use std::collections::HashMap;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! use twobitreader::reverse_complement;
//! use rayon::prelude::*;
//!
//! fn stranded(seq: String, strand: char) -> String {
//!     if strand == '+' { seq } else { reverse_complement(seq) }
//! }
//!
//! let transcripts = [                   // (transcript_id, chromosome, exons)
//!     ("ENST00000407983.7", "chr2", '+', vec![(264899, 265007),     // Exon 1
//!                                             (271865, 271939),     // Exon 2
//!                                             (272036, 272557)]),   // Exon 3
//!     ("ENST00000319331.4", "chr3", '+', vec![(3799430, 3799919),   // Exon 1
//!                                             (3844363, 3849834)]), // Exon 2
//!     /* ... */
//! ];
//! // Concatenate exons and then reverse-complement if necessary.
//! // (Correct if exons listed in genome-coordinate order.)
//! let seqs = transcripts.into_par_iter()
//!     .map(|(id, chrom, strand, exons)| (id, stranded(tbr.concat(chrom, exons), strand)))
//!     .collect::<HashMap<_, _>>();      // HashMap<&str, String>
//! let seq = &seqs["ENST00000407983.7"]; // -> &String to transcript sequence
//! # Ok::<(), io::Error>(())
//! ```
//!
//! **Cold files** are an order of magnitude slower to access than files already in memory ("hot").
//! Use prefetching to dramatically improve single-threaded speed:
//! ```no_run
//! # use std::io;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! let exons = [("chr1", 10000, 10200),
//!              ("chr1", 10500, 10700), /* ... */ ];
//! tbr.prefetch(&exons);             // Ask the operating system to start paging this data from disk.
//! let seqs = tbr.get_batch(&exons); // Access the memory as it arrives.
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
//! - **hot** runs represent repeated or interactive dna extraction scenarios;
//! - **cold** runs represent a first run of a genomics pipeline, bound by disk speed;
//! - **prefetch** runs are cold but with a `prefetch` call preceding extraction.
//!
//! The table below shows running times in milliseconds. Experimental details are `BENCH.md`.
//!
//! | EXONS | 1-thread / hot | 1-thread / cold | 1-thread / prefetch | 16-thread / hot | 16-thread / cold |
//! |---|---:|---:|---:|---:|---:|
//! | **twobitreader** (rust)     | 90    | 1,600  | 180 | 11    | 190   |
//! | **py2bit** (C, python)      | 220   | 2,800  | n/a | n/a   | n/a   |
//! | **GenomeKit** (C++, python) | 340   | 2,400  | n/a | n/a   | n/a   |
//! | **twobit** (rust)           | 390   | 3,500  | n/a | n/a   | n/a   |
//! | **twobitToFa** (C)          | 1,000 | 4,500  | n/a | 340   | 850   |
//! | **twobitreader** (python)   | 7,200 | 13,000 | n/a | 1,900 | 2,300 |
//! | **Biopython** (python)      | 8,100 | 11,000 | n/a | n/a   | n/a   |
//!
//! | TRANSCRIPTS | 1-thread / hot | 1-thread / cold | 1-thread / prefetch | 16-thread / hot | 16-thread / cold |
//! |---|---:|---:|---:|---:|---:|
//! | **twobitreader** (rust)     | 160    | 2,500  | 240 | 18    | 210   |
//! | **py2bit** (C, python)      | 490    | 3,600  | n/a | n/a   | n/a   |
//! | **GenomeKit** (C++, python) | 720    | 2,900  | n/a | n/a   | n/a   |
//! | **twobit** (rust)           | 930    | 1,900  | n/a | n/a   | n/a   |
//! | **twobitToFa** (C)          | 2,100  | 6,900  | n/a | 840   | 1,200 |
//! | **twobitreader** (python)   | 13,000 | 21,000 | n/a | 2,800 | 3,500 |
//! | **Biopython** (python)      | 19,000 | 25,000 | n/a | n/a   | n/a   |
//!
//! # Dependencies
//!
//! * `byteorder` for handling endian-ness
//! * `memmap2` for memory mapping the 2bit file
//! * `seq-macro` for generating 2bit decoder lookup table
//! * `libc` for prefetching file ranges on Apple targets
//! * `windows-sys` for prefetching file ranges on Windows targets

// Standard library
use std::borrow::Borrow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Cursor, ErrorKind::InvalidData, Read};
use std::iter::zip;
use std::mem::size_of;
#[cfg(debug_assertions)]
use std::mem::MaybeUninit;
use std::ops::Range;
use std::path::Path;
use std::sync::OnceLock;

// Crate modules
mod decode;
use decode::{decode, NUCS_PER_U8};
mod prefetch;
use prefetch::PrefetchBatcher;

// Dependencies
use byteorder::{BigEndian, ByteOrder, LittleEndian, ReadBytesExt};
use memmap2::{Mmap, MmapOptions};
use seq_macro::seq;

// Big-endian 2bit files are supported by this crate, but big-endian compile targets are not.
#[cfg(target_endian = "big")]
compile_error!("twobitreader is not yet implemented for big-endian targets.");

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Endianness {
    Big,
    Little,
}

/// A reader for a single [2bit file](http://genome.ucsc.edu/FAQ/FAQformat.html#format7).
///
#[derive(Debug)]
pub struct TwobitReader {
    file: File,                          // Open 2bit file; kept for issuing prefetch hints.
    mmap: Mmap,                          // Memory map of entire 2bit file.
    masked: bool,                        // Whether to apply the lowercase mask.
    endianness: Endianness,              // Whether the 2bit file is big-endian.
    seqs: Vec<TwobitSequence>,           // Sequence data, in same order as in the file.
    seq_by_name: HashMap<String, usize>, // Lookup sequence index (seqs[index]) by name.
}

/// A specific sequence record within the 2bit file.
///
#[derive(Debug)]
struct TwobitSequence {
    data: OnceLock<TwobitSequenceData>, // Heavy data loaded lazily from the sequence record.
    name: String,                       // Sequence name ("chr3", etc.)
    data_offset: u64,                   // Offset to data block
}

/// Details of a sequence record located deeper in the file, such as block indices.
///
#[derive(Debug)]
struct TwobitSequenceData {
    dna_offset: usize, // Offset (within mmap) to first byte of packed DNA
    dna_bytes: usize,  // Number of bytes (not nucleotides!) of packed DNA
    dna_len: usize,    // Number of nucleotides (not bytes!) in DNA sequence
    nblocks: Blocks,   // N block starts and ends [start0, end0, start1, end1, ...]
    masks: Blocks,     // Mask block starts and ends
}

// Sorted parallel arrays where (starts[i], ends[i]) is the range of one N-block or
// mask-block (lowercase block) in a sequence record.
#[derive(Debug)]
struct Blocks {
    starts: Vec<u32>,
    ends: Vec<u32>,
}

impl TwobitReader {
    /// Opens a 2bit file for reading.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?; // Human genome, build 38
    /// let seq = tbr.get("chr2", 10000, 10010);    // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::open_impl(path, false)
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
    /// let tbr = TwobitReader::open_masked("hg38.2bit")?; // Human genome, build 38
    /// let seq = tbr.get("chr2", 10000, 10010);           // -> "CGTATcccac" (mixed case)
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn open_masked<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::open_impl(path, true)
    }

    // Implements opening the file and initializing a TwobitReader instance for it.
    fn open_impl<P: AsRef<Path>>(path: P, masked: bool) -> io::Result<Self> {
        // Map the file and read the header and mask arrays in the order they appear in the file
        let file = File::open(path)?;

        // SAFETY: Mapping a file is inherently unsafe, because another process could modify the
        // the file while it is mapped. There is nothing we can do about that.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let endianness = read_endianness(&mmap)?;
        let seqs = read_seqs(&mmap, endianness)?;

        // Allow fast lookup of sequences by name
        let seq_by_name = HashMap::from_iter(seqs.iter().enumerate().map(|(i, seq)| (seq.name.clone(), i)));
        if seq_by_name.len() != seqs.len() {
            return Err(io::Error::new(InvalidData, "duplicate sequence name detected."));
        }

        Ok(Self { file, mmap, masked, endianness, seqs, seq_by_name })
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
    /// Panics if an IO error occurs due to a malformed file.
    ///
    pub fn seq_len<N: AsRef<str>>(&self, name: N) -> usize {
        self.get_seq_data_by_name(name).dna_len
    }

    /// Extracts range `start..end` (0-based, exclusive end) from the named sequence record.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let seq = tbr.get("chr2", 10000, 10010); // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or the range was invalid.
    /// Panics if an IO error occurs due to a malformed file.
    ///
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
    /// tbr.get_into("chr2", 10000, 10010, &mut dst); // dst = "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`get`](Self::get).
    ///
    pub fn get_into<N: AsRef<str>>(&self, name: N, start: usize, end: usize, dst: &mut String) {
        // Check range before trying to do any might-panic arithmetic with (start, end)
        let seq = self.get_seq_data_by_name(name);
        check_range(seq, start, end);

        // Prepare an empty buffer with sufficient capacity.
        // SAFETY: clearing dst ensures that its visible portion remains valid utf8
        // until decode_and_append lengthens it.
        dst.clear();
        let buf = unsafe { dst.as_mut_vec() };
        buf.reserve_exact(end - start);
        self.decode_and_append(seq, start, end, buf);
    }

    /// Returns an iterator that calls [`get`](Self::get) for each query in the batch, for convenience.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// # let tbr = TwobitReader::open("hg38.2bit")?;
    /// let exons = [("chr1", 10000, 15000),
    ///              ("chr2", 30000, 35000), /* ... */ ];
    /// let seq = tbr.get_batch(&exons);
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`get`](Self::get).
    ///
    pub fn get_batch<'a, N, T, I>(&'a self, batch: I) -> impl Iterator<Item = String> + 'a
    where
        // The Borrow is for flexibility. We borrow a reference even if batch was passed by value.
        // The lifetime bounds express that the returned iterator must not outlive self (captured below).
        N: AsRef<str>,
        T: Borrow<(N, usize, usize)>,
        I: IntoIterator<Item = T>,
        I::IntoIter: 'a,
    {
        batch.into_iter().map(|item| {
            let (chrom, start, end) = item.borrow();
            self.get(chrom, *start, *end)
        })
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
    /// let seq = tbr.get_inclusive("chr2", 10001, 10010); // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`get`](Self::get).
    ///
    pub fn get_inclusive<N: AsRef<str>>(&self, name: N, start: usize, end: usize) -> String {
        // Check valid start and then convert range to 0-based exclusive.
        check_start_inclusive(start, 1);
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
    /// tbr.get_inclusive_into("chr2", 10001, 10010, &mut dst); // dst = "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`get`](Self::get).
    ///
    pub fn get_inclusive_into<N: AsRef<str>>(&self, name: N, start: usize, end: usize, dst: &mut String) {
        // Check valid start and then convert range to 0-based exclusive.
        check_start_inclusive(start, 1);
        self.get_into(name, start - 1, end, dst);
    }

    /// A version of [`get_batch`](Self::get_batch) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// # let tbr = TwobitReader::open("hg38.2bit")?;
    /// let exons = [("chr1", 10001, 15000),
    ///              ("chr2", 30001, 35000), /* ... */ ];
    /// let seq = tbr.get_batch_inclusive(&exons);
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`get`](Self::get).
    ///
    pub fn get_batch_inclusive<'a, N, T, I>(&'a self, batch: I) -> impl Iterator<Item = String> + 'a
    where
        N: AsRef<str>,
        T: Borrow<(N, usize, usize)>,
        I: IntoIterator<Item = T>,
        I::IntoIter: 'a,
    {
        batch.into_iter().map(|item| {
            let (chrom, start, end) = item.borrow();
            self.get_inclusive(chrom, *start, *end)
        })
    }

    /// Concatenates a batch of sequence ranges into a single string.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let ranges = [(10000, 10006), (10006, 10010)]; // "CGTATC" "CCAC"
    /// let seq = tbr.concat("chr2", &ranges);         // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// This version performs one allocation only, of the required total length.
    /// However, it requires `ranges` to be iterable twice.
    /// See [`concat_iter`](Self::concat_iter) for a version that grows rather than pre-allocates.
    ///
    /// # Panics
    ///
    /// Panics if the sequence name was not found or if any range was invalid.
    /// Panics if an IO error occurs due to a malformed file.
    ///
    pub fn concat<N, R>(&self, name: N, ranges: R) -> String
    where
        N: AsRef<str>,
        R: AsRef<[(usize, usize)]>,
    {
        self.concat_impl(name, ranges, 0)
    }

    // Implements concat for 0-based-exclusive and 1-based-inclusive ranges.
    fn concat_impl<N, R>(&self, name: N, ranges: R, base: usize) -> String
    where
        N: AsRef<str>,
        R: AsRef<[(usize, usize)]>,
    {
        // Compute total length by summing range lengths; also check each range.
        let seq = self.get_seq_data_by_name(name);
        let ranges = ranges.as_ref();
        let total_len = ranges
            .iter()
            .map(|&(start, end)| {
                check_start_inclusive(start, base); // Check start >= base before subtracting
                check_range(seq, start - base, end);
                end - (start - base)
            })
            .sum();

        // Pre-size a byte buffer to the total length needed.
        let mut buf = Vec::<u8>::new();
        buf.reserve_exact(total_len);

        // Decode each interval into its respective slice of the buffer.
        for &(start, end) in ranges.iter() {
            self.decode_and_append(seq, start - base, end, &mut buf);
        }

        // SAFETY: buf is valid utf8 here because decode_and_append writes through whatever spare
        // capacity it adds before setting the length.
        debug_assert!(buf.is_ascii(), "decoded nucleotides were not valid ascii");
        unsafe { String::from_utf8_unchecked(buf) }
    }

    /// A version of [`concat`](Self::concat) that iterates through `ranges`, consuming it.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use std::iter::zip;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let starts = [10000, 10006];
    /// let ends = [10006, 10010];
    /// let seq = tbr.concat_iter("chr2", zip(starts, ends)); // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// This version iterates through `ranges` exactly once, growing the output String as needed.
    ///
    /// # Panics
    ///
    /// See [`concat`](Self::concat).
    ///
    pub fn concat_iter<N, I>(&self, name: N, ranges: I) -> String
    where
        N: AsRef<str>,
        I: IntoIterator<Item = (usize, usize)>,
    {
        let seq = self.get_seq_data_by_name(name);
        let mut buf = Vec::<u8>::new();

        // Decode each interval into its respective slice of the buffer.
        for (start, end) in ranges.into_iter() {
            check_range(seq, start, end);
            buf.reserve(end - start);
            self.decode_and_append(seq, start, end, &mut buf);
        }

        // SAFETY: see concat()
        debug_assert!(buf.is_ascii(), "decoded nucleotides were not valid ascii");
        unsafe { String::from_utf8_unchecked(buf) }
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
    /// let seq = tbr.concat_inclusive("chr2", ranges); // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`concat`](Self::concat).
    ///
    pub fn concat_inclusive<N, R>(&self, name: N, ranges: R) -> String
    where
        N: AsRef<str>,
        R: AsRef<[(usize, usize)]>,
    {
        self.concat_impl(name, ranges, 1)
    }

    /// A version of [`concat_iter`](Self::concat_iter) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use std::iter::zip;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let starts = [10001, 10007];
    /// let ends = [10006, 10010];
    /// let seq = tbr.concat_iter_inclusive("chr2", zip(starts, ends)); // -> "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// See [`concat`](Self::concat).
    ///
    pub fn concat_iter_inclusive<N, I>(&self, name: N, ranges: I) -> String
    where
        N: AsRef<str>,
        I: IntoIterator<Item = (usize, usize)>,
    {
        // Check each start, and convert range to 0-based exclusive.
        self.concat_iter(
            name,
            ranges.into_iter().map(|(start, end)| {
                check_start_inclusive(start, 1);
                (start - 1, end)
            }),
        )
    }

    /// Asks the operating system ensure that the data for all the ranges gets paged into memory,
    /// and returns without waiting for the data to arrive.
    ///
    /// Use this method only if are reading a "cold" file that is not already paged in from disk.
    /// Accessing a cold file is almost entirely IO-bound. This function provides a hint to the
    /// operating system, specifying what data you'll soon be asking for. It does not affect the
    /// output of subsequent calls to [`get`](Self::get) and [`concat`](Self::concat), only their
    /// speed. For maximum benefit, prefetch the largest batch of ranges that you can in one call.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// # let tbr = TwobitReader::open("hg38.2bit")?;
    /// let exons = [("chr1", 10000, 15000),
    ///              ("chr2", 30000, 35000), /* ... */ ];
    /// tbr.prefetch(&exons);             // Ask the operating system to start paging this data from disk.
    /// let seqs = tbr.get_batch(&exons); // Access the memory as it arrives.
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if any sequence name was not found or if `start > end``.
    /// Panics if an IO error occurs while processing any newly-accessed sequence records.
    ///
    pub fn prefetch<N, T, I>(&self, args: I)
    where
        N: AsRef<str>,
        T: Borrow<(N, usize, usize)>,
        I: IntoIterator<Item = T>,
    {
        // Bucket the ranges by sequence name. Most of the new vecs will remain empty and unallocated.
        let mut ranges_by_seq: Vec<Vec<Range<usize>>> = self.seqs.iter().map(|_| Vec::new()).collect();
        for item in args.into_iter() {
            let (name, start, end) = item.borrow();
            let seq_index = self.seq_by_name.get(name.as_ref()).expect("sequence name not found");
            assert!(start <= end, "invalid range (start > end");
            ranges_by_seq[*seq_index].push(*start..*end);
        }

        let mut prefetch = PrefetchBatcher::new(&self.file, &self.mmap);

        // Prefetch the leading bytes of each sequence's data block, since those will be read for sure.
        // Makes prefetch on a cold file about 5% faster. For *masked* files, 4kb may be too conservative.
        const BLOCKS_PREFETCH_BYTES: usize = 4 * 1024;
        for (seq, ranges) in zip(self.seqs.iter(), ranges_by_seq.iter()) {
            if !ranges.is_empty() {
                let data_offset = seq.data_offset as usize;
                prefetch.push(data_offset..data_offset + BLOCKS_PREFETCH_BYTES);
            }
        }
        prefetch.flush();

        // Prefetch the ranges in each bucket.
        for (seq_index, mut ranges) in ranges_by_seq.into_iter().enumerate() {
            // If empty, do NOT call get_seq_data; that would process the blocks of an unused sequence.
            if ranges.is_empty() {
                continue;
            }

            // Before sending hints to the operating system, aggregate the ranges to reduce
            // overhead of repeated hints, which can be a few microseconds each.
            // The aggregated ranges will be over-complete, but tuning the threshold below
            // gives a significant performance improvement. The threshold below means that
            // two consecutive ranges need a gap > 16 kilobytes to avoid being aggregated.
            const AGGREGATION_GAP_THRESHOLD: usize = NUCS_PER_U8 * 16 * 1024;
            ranges.sort_unstable_by_key(|range| range.start);
            let mut i = 0;
            for j in 1..ranges.len() {
                if ranges[i].end.saturating_add(AGGREGATION_GAP_THRESHOLD) >= ranges[j].start {
                    ranges[i].end = ranges[i].end.max(ranges[j].end);
                } else {
                    i += 1;
                    ranges[i] = ranges[j].clone();
                }
            }
            ranges.truncate(i + 1);

            // Send the aggregated hints to the operating system. This also pages in and processes
            // the blocks arrays for the corresponding sequence record (by calling get_seq_data).
            let seq = self.get_seq_data_by_index(seq_index);
            for range in ranges.iter_mut() {
                if range.start < range.end {
                    let first_byte = seq.dna_offset + range.start / NUCS_PER_U8;
                    let last_byte = seq.dna_offset + range.end.div_ceil(NUCS_PER_U8);
                    prefetch.push(first_byte..last_byte);
                }
            }
        }
        prefetch.flush();
    }

    // Returns a reference to the [`TwobitSequence`] for the named sequence record.
    // Panics if no sequence record has that name.
    fn get_seq_data_by_name<N: AsRef<str>>(&self, name: N) -> &TwobitSequenceData {
        let index = *self.seq_by_name.get(name.as_ref()).expect("sequence name not found");
        self.get_seq_data_by_index(index)
    }

    // Returns a reference to the [`TwobitSequenceData`] for the named sequence record.
    // Panics if no sequence record has that name.
    fn get_seq_data_by_index(&self, index: usize) -> &TwobitSequenceData {
        let seq = &self.seqs[index];
        seq.data.get_or_init(|| read_seq_data(&self.mmap, self.masked, self.endianness, seq.data_offset))
    }

    // Decodes sequence [start..end] and appends it to dst.
    // Requires dst to already have the capacity to hold the decoded bytes.
    fn decode_and_append(&self, seq: &TwobitSequenceData, start: usize, end: usize, dst: &mut Vec<u8>) {
        if start >= end {
            return;
        }

        // Check that capacity is already sufficient.
        debug_assert!(dst.capacity() - dst.len() >= end - start);

        // Slice spanning all packed 2bit DNA data for this sequence record.
        let dna = &self.mmap[seq.dna_offset..seq.dna_offset + seq.dna_bytes];

        // Decode into the leading bytes of dst's unused capacity. Writing to unused capacity
        // ensures that the used portion never contains uninitialized bytes or broken utf8.
        let range_len = end - start;
        let buf_uninit = &mut dst.spare_capacity_mut()[..range_len];

        // Poison the buffer in debug builds just in case the logic of decode() violates the
        // assumption that all bytes of buf get written to.
        #[cfg(debug_assertions)]
        buf_uninit.fill(MaybeUninit::new(0xfe));

        // Decode packed 2-bit dna from the given start position.
        decode(start, dna, buf_uninit);

        // SAFETY: decode() wrote every byte of its destination slice. If we arrived
        // here without a panic, then range_len additional bytes are now valid ASCII
        // and can be safely appended (via set_len) to dst.
        let old_len = dst.len();
        unsafe {
            dst.set_len(old_len + range_len);
        }
        let buf_init = &mut dst[old_len..];

        // Apply N-block and lowercase masks to the newly-written portion of dst, in-place.
        search_blocks(&seq.nblocks, start, buf_init, nblock_fill);
        search_blocks(&seq.masks, start, buf_init, mask_fill);

        // Catch any decoding errors or unwritten bytes.
        debug_assert!(buf_init.is_ascii());
    }
}

// Panics if (start, end) is an invalid 0-based exclusive range for this TwobitSequence.
// (The error message does not format the actual start/end values that were passed in.
//  This is because the 0-based get functions are used to implement the 1-based get
//  functions, and so the start/end seen here may not be the actual ones users provided,
//  which would result in more confusion for the user than simply omitting the numbers.)
#[inline]
fn check_range(seq: &TwobitSequenceData, start: usize, end: usize) {
    if start > end {
        panic!("invalid range (start > end)");
    }
    if end > seq.dna_len {
        panic!("invalid end (end > dna_len)");
    }
}

// Panics if start is not valid for an inclusive range
#[inline]
fn check_start_inclusive(start: usize, base: usize) {
    if start < base {
        debug_assert_eq!(base, 1, "expected base=1 but found {base}");
        panic!("invalid start (0) for 1-based range");
    }
}

fn read_endianness(mmap: &Mmap) -> io::Result<Endianness> {
    // Twobit file signature as little-endian (native) or big-endian (non-native)
    const SIGNATURE_LILEND: u32 = 0x1A412743;
    const SIGNATURE_BIGEND: u32 = 0x4327411A;

    // Start reading the file, using appropriate byte order
    let signature = Cursor::new(mmap).read_u32::<LittleEndian>()?;
    match signature {
        SIGNATURE_LILEND => Ok(Endianness::Little),
        SIGNATURE_BIGEND => Ok(Endianness::Big),
        _ => Err(io::Error::new(InvalidData, "failed to read a valid 2bit file signature.")),
    }
}

// Reads the sequence records from a .2bit file.
// The resulting vec of TwobitSequences is intended to be assigned to the TwobitReader::seqs field.
// Block intervals and sequence data are not paged in at this point, only the header.
fn read_seqs(mmap: &Mmap, endianness: Endianness) -> io::Result<Vec<TwobitSequence>> {
    match endianness {
        Endianness::Little => read_seqs_endian::<LittleEndian>(mmap),
        Endianness::Big => read_seqs_endian::<BigEndian>(mmap),
    }
}

// Implements read_seqs after the endian-ness of the file has been determined from the signature.
fn read_seqs_endian<B: ByteOrder>(mmap: &Mmap) -> io::Result<Vec<TwobitSequence>> {
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
    // (Deliberate choice not to prefetch; a typical 2bit file's header fits in a single page.)
    let mut seqs = Vec::with_capacity(num_seqs);
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
        let data_offset = cursor.read_u32::<B>()? as u64;
        let curr_offset = cursor.position();
        if data_offset <= curr_offset {
            return Err(io::Error::new(InvalidData, "invalid data offset; file may be malformed."));
        }

        // The block indices and dna_offsets are loaded lazily for each sequence record when first accessed
        // through the OnceLock below. This avoids cost of paging in block data that may never be used,
        // which tends to be scattered around the file; important for opening cold files quickly.
        let data = OnceLock::new();
        seqs.push(TwobitSequence { data, name, data_offset });
    }
    Ok(seqs)
}

// Reads the block data from a single sequence record at the given data_offset.
// This must be done prior to decoding any DNA from the given sequence record.
fn read_seq_data(mmap: &Mmap, masked: bool, endianness: Endianness, data_offset: u64) -> TwobitSequenceData {
    match endianness {
        Endianness::Little => read_seq_data_endian::<LittleEndian>(mmap, masked, data_offset),
        Endianness::Big => read_seq_data_endian::<BigEndian>(mmap, masked, data_offset),
    }
}

// Implements read_seq_data for a specific endian-ness.
fn read_seq_data_endian<B: ByteOrder>(mmap: &Mmap, masked: bool, data_offset: u64) -> TwobitSequenceData {
    // Prepare to read the sequence's data block.
    let mut cursor = Cursor::new(mmap);
    cursor.set_position(data_offset);

    // Read and process the block arrays. Do not touch the packed DNA data.
    let dna_len = cursor.read_u32::<B>().expect("Failed to read DNA length from 2bit file.") as usize;
    let nblocks = read_blocks::<B>(&mut cursor, true);
    let masks = read_blocks::<B>(&mut cursor, masked);

    // Skip reserved bytes
    cursor.set_position(cursor.position() + size_of::<u32>() as u64);

    // Check blocks are valid.
    check_blocks(&nblocks, dna_len);
    check_blocks(&masks, dna_len);

    // Get offset and length of packed DNA, and make sure it's within the file's size.
    let dna_offset = cursor.position() as usize;
    let dna_bytes = dna_len.div_ceil(NUCS_PER_U8);

    assert!(
        dna_offset + dna_bytes <= cursor.get_ref().len(),
        "Failed to read DNA from 2bit file. DNA data was truncated."
    );

    TwobitSequenceData { dna_offset, dna_bytes, dna_len, nblocks, masks }
}

// Reads a Twobit block section from the given cursor, returning a Blocks struct.
// If used=false, then the cursor is simply advanced, and an empty Blocks struct is returned.
fn read_blocks<B: ByteOrder>(cursor: &mut Cursor<&Mmap>, used: bool) -> Blocks {
    let num_blocks = cursor.read_u32::<B>().expect("Failed to read number of blocks from 2bit file.") as usize;
    let mut starts = Vec::new();
    let mut ends = Vec::new();

    // Check to avoid huge allocation on corrupt or malicious block count.
    assert!(2 * num_blocks * size_of::<u32>() <= cursor.get_ref().len(), "num_blocks too large, may be corrupt");

    if used {
        // Read starts.
        starts.resize(num_blocks, 0);
        cursor.read_u32_into::<B>(&mut starts).expect("Failed to read block starts from 2bit file.");

        // Read sizes and convert them to ends.
        ends.resize(num_blocks, 0);
        cursor.read_u32_into::<B>(&mut ends).expect("Failed to read block ends from 2bit file.");
        for (start, end) in zip(&starts, &mut ends) {
            *end += *start;
        }
    } else {
        // Advance the cursor, without reading the data itself.
        let new_offset = cursor.position() as usize + 2 * num_blocks * size_of::<u32>();
        assert!(new_offset <= cursor.get_ref().len(), "Failed to read blocks from 2bit file. Block indices truncated.");
        cursor.set_position(new_offset as u64);
    }

    Blocks { starts, ends }
}

#[cfg(debug_assertions)]
fn check_blocks(blocks: &Blocks, dna_len: usize) {
    assert!(blocks.starts.is_sorted(), "block start indices are not sorted; file may be malformed.");
    assert!(blocks.ends.is_sorted(), "block end indices are not sorted; file may be malformed.");
    for (start, end) in zip(&blocks.starts, &blocks.ends) {
        assert!(start < end, "block start >= end; file may be malformed.");
        assert!(*start as usize <= dna_len, "block start beyond DNA extents; file may be malformed.");
        assert!(*end as usize <= dna_len, "block end beyond DNA extents; file may be malformed.");
    }
}

#[cfg(not(debug_assertions))]
fn check_blocks(_blocks: &Blocks, _dna_len: usize) {}

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

// Compile-time generated lookup table of nucleotide complements for all u8 byte values.
const NUC_COMPLEMENT_U8: [u8; 256] = seq!(i in 0..256 {[#(
    match i {
       b'A' => b'T', b'a' => b't',
       b'C' => b'G', b'c' => b'g',
       b'G' => b'C', b'g' => b'c',
       b'T' => b'A', b't' => b'a',
       b'N' => b'N', b'n' => b'n',
       _    => b'?',
    },
)*]});

/// Returns a reverse-complemented version of the input sequence.
///
/// Useful for higher-level code that wants to assemble strand-sensitive transcripts.
///
/// The string is modified in-place and returned, so no allocation takes place.
///
/// Note that `seq` must contain only characters from `ACGTNacgtn`. Otherwise,
/// the invalid character will assert (debug) or be replaced with '?' (release).
///
pub fn reverse_complement(mut seq: String) -> String {
    // SAFETY: this is safe if dna contains ACGTN bytes, as the buffer will remain
    // valid utf8 at every step.
    if !seq.is_empty() {
        let buf = unsafe { seq.as_mut_vec() };
        buf.reverse();
        for byte in buf {
            *byte = NUC_COMPLEMENT_U8[*byte as usize];
        }
    }
    debug_assert!(!seq.contains('?'), "invalid character in DNA string.");
    seq
}
