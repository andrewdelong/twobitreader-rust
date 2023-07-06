//! This crate provides fast DNA sequence extraction from 2bit files, a
//! [standard format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) in bioinformatics.
//! 
//! The motivation for this crate is speed.
//! Extracting sequences is 1.5-30x faster than the best alternative, depending on use case.
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
//! **Batch extraction** is fast.
//! It works by iterating over (chrom, start, end) triplets in parallel:
//! ```no_run
//! # use std::io;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! let args = [("chr1", 10000, 15000),
//!             ("chr1", 30000, 35000), /* ... */ ];
//! let seqs = tbr.get_batch(args);     // Vec<String>
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
//! For example, assembling a batch of spliced transcripts in parallel:
//! ```no_run
//! # use std::io;
//! # use std::collections::HashMap;
//! # use twobitreader::TwobitReader;
//! # let tbr = TwobitReader::open("hg38.2bit")?;
//! use rayon::prelude::*;
//! let transcripts = [     // (transcript_id, chromosome, exons)
//!     ("ENST00000407983.7", "chr2", vec![(264899, 265007),      // Exon 1 (start, end)
//!                                        (271865, 271939),      // Exon 2 (start, end)
//!                                        (272036, 272557)]),    // Exon 3 (start, end)
//!     ("ENST00000319331.4", "chr3", vec![(3799430, 3799919),    // Exon 1 (start, end)
//!                                        (3844363, 3849834)]),  // Exon 2 (start, end)
//!     /* ... */
//! ];
//! let seqs = transcripts.into_par_iter()
//!                       .map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons)))
//!                       .collect::<HashMap<_, _>>();  // HashMap<&str, String>
//! let seq = &seqs["ENST00000407983.7"];               // Look up transcript sequence 
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
//! * `num_cpus` for spawning parallel workers
//! * `seq-macro` for generating 2bit decoder lookup table

// Standard library
use std::fs::File;
use std::io::{self, Cursor, Read};
use std::ops::Index;
use std::mem::size_of;
use std::ptr;
use std::vec;
use std::iter::{zip, repeat_with};
use std::path::Path;
use std::panic;
use std::slice;
use std::collections::HashMap;
use core::borrow::Borrow;
use core::hash::Hash;

// Crate modules
mod parallel;
use parallel::{parallel_for, try_parallel_for};
mod decode;
use decode::{decode, NUCS_PER_U8};

// Dependencies
use byteorder::{BigEndian, LittleEndian, ReadBytesExt, ByteOrder};
use memmap2::{Mmap, MmapOptions};

// Big-endian 2bit files are supported by this crate, but big-endian compile targets are not.
#[cfg(target_endian = "big")]
compile_error!("twobitreader is not yet implemented for big-endian targets.");

/// A reader for a single [2bit file](http://genome.ucsc.edu/FAQ/FAQformat.html#format7).
/// 
#[allow(dead_code)]                       // Suppress "mmap unused" warning; only referenced via *const u8s.
pub struct TwobitReader {    
    mmap: Mmap,                           // Memory map of entire 2bit file
    seqs: Vec<TwobitSequence>,            // Sequence data, in same order as in the file
    seq_by_name: HashMap<String, usize>,  // Lookup sequence index (seqs[index]) by name.
}

/// A reader for a specific sequence record within the 2bit file.
/// 
pub struct TwobitSequence {
    dna_ptr: *const u8,     // Pointer to first byte of packed DNA sequence data
    dna_bytes: usize,       // Number of bytes (not nucleotides!) of packed DNA
    dna_len: usize,         // Number of nucleotides (not bytes!) in DNA sequence
    nblocks: Blocks,        // N block starts and ends [start0, end0, start1, end1, ...]
    masks: Blocks,          // Mask block starts and ends
    name: String,           // Sequence name ("chr3", etc.)
}

// Sorted parallel arrays where (starts[i], ends[i]) is the range of one N-block or
// lowercase-block in a sequence record.
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
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;  // Human genome, build 38
    /// let seq = tbr.get("chr2", 10000, 10010);     // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<TwobitReader> {
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
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open_masked("hg38.2bit")?;  // Human genome, build 38
    /// let seq = tbr.get("chr2", 10000, 10010);            // "CGTATcccac" (mixed case)
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn open_masked<P: AsRef<Path>>(path: P) -> io::Result<TwobitReader> {
        Self::open_inner(path, true)
    }

    // Implements opening the file and initializing a TwobitReader instance for it.
    fn open_inner<P: AsRef<Path>>(path: P, use_mask: bool) -> io::Result<TwobitReader> {
        // Map the file and read the header and mask arrays in the order they appear in the file
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let seqs = read_seqs(&mmap, use_mask)?;

        // Allow fast lookup of sequences by name
        let seq_by_name = HashMap::from_iter(seqs.iter().enumerate().map(|(i, seq)| (seq.name.clone(), i)));
        if seq_by_name.len() != seqs.len() {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "TwobitReader: duplicate sequence name detected."));
        }

        Ok(TwobitReader{mmap, seqs, seq_by_name})
    }
    
    /// Returns the number of sequence records in the file.
    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    /// Iterates over the [`TwobitSequence`] records, in the order they appear in the file.
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// for tbs in &tbr {
    ///     println!("{} {}", tbs.name(), tbs.len());
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    /// 
    pub fn iter(&self) -> impl Iterator<Item=&TwobitSequence> + '_ {
        self.seqs.iter()
    }

    /// Iterates over the sequence record names, in the order they appear in the file.
    pub fn iter_names(&self) -> impl Iterator<Item=&str> + '_ {
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

    /// Extracts range `start..end` (0-based, exclusive end) from the named sequence record.
    /// 
    /// # Examples
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
    pub fn get<N: AsRef<str>>(&self, name: N, start: usize, end: usize) -> String {
        self[name.as_ref()].get(start, end)
    }

    /// A version of [`get`](Self::get) where the result is stored in `dst`.
    /// 
    /// # Examples
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
        self[name.as_ref()].get_into(start, end, dst);
    }

    /// Extracts a batch of sequences by calling [`get`](Self::get) on every 
    /// `(name, start, end)` item returned by iterating over `iter`.
    /// 
    /// Iteration is parallel and potentially much faster than calling
    /// [`get`](Self::get) sequentially.
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let args = [("chr2", 10000, 10006),
    ///             ("chr2", 10006, 10010), /* ... */];
    /// let seqs = tbr.get_batch(args);  // ["CGTATC", "CCAC", ...]
    /// # Ok::<(), io::Error>(())
    /// ```
    /// 
    /// Internally this method is very similar to using `rayon` (see below), except a
    /// global thread pool is not maintained.
    /// ```no_run
    /// # use rayon::prelude::*;
    /// # use twobitreader::TwobitReader;
    /// # let tbr = TwobitReader::open("hg38.2bit").unwrap();
    /// # let args = [("chr2", 10000, 10006),
    /// #             ("chr2", 10006, 10010)];
    /// let seqs = args.into_par_iter().map(
    ///     |(name, start, end)| tbr.get(name, start, end)
    /// ).collect::<Vec<_>>();  // Vec<String>
    /// ```
    /// 
    /// # Panics
    ///
    /// Panics if any sequence name was not found or if any range was invalid.
    /// 
    pub fn get_batch<N, I>(&self, iter: I) -> Vec<String>
    where
        N: AsRef<str> + Sync + Send,   // N = String, &String, &str, etc...
        I: IntoIterator<Item = (N, usize, usize)>,
    {
        // Collect iter into a Vec. This makes sizing, preallocating, and chunking easy.
        // Comment on the overhead of collecting iter into a Vec:
        // - If iter is vec::IntoIter<(N, usize, usize)>, collecting iter back into a Vec is free.
        // - Specifically, Vec::from_iter is specialized for vec::IntoIter, and will
        //   simply claim the capacity for itself in O(1), with no overhead.
        // - For details, see SpecFromIter<T, IntoIter<T>>::from_iter for Vec<T>.
        // - However, if collect is forced to iterate, then this is the only step
        //   that iterates sequentially, and might be slower than rayon as a result.
        let mut args = Vec::from_iter(iter.into_iter());
        let mut seqs = Vec::with_capacity(args.len());
        
        // Chunk the input/output vecs.
        // Initially, the output seqs vec has size zero, and all its memory is uninitialized capacity.
        // The seq_chunks provides a view into that capacity, in the form of MaybeUninit<String> items.
        // If all chunks extract without panicking, seqs.set_len() will claim the now-initialized Strings.
        // Comment on speed of using fixed-sized chunks:
        // - Using a granular fixed chunk size tends to be faster than one-chunk-per-thread.
        //   Why? Because there tends to be large variation in the time to process each chunk,
        //   and the one-chunk-per-thread approach tends to leave many threads idle while the
        //   slowest remaining chunks are processed. This is especially true when the mapped
        //   file is 'cold', because there can be huge variability in the delays incurred for
        //   paging sequence data in from disk.
        // - If chunk size is too big, then small batches won't benefit from parallelism, and
        //   this is extra important for "cold" files that need to be paged in.
        // - If chunk size is very small, the synchronization overhead becomes significant.
        const CHUNK_SIZE: usize = 100;
        let arg_chunks = args.chunks_mut(CHUNK_SIZE);
        let seq_chunks = seqs.spare_capacity_mut().chunks_mut(CHUNK_SIZE);

        // Main loop. Call get() on each input, and store the result.
        // - Use MaybeUninit::write to emplace the String into the output slot represented by seq.
        // - Use drop_in_place to consume the name and dealloc if applicable.
        //   For example, if N = String, then we own the string and drop_in_place will release its
        //   u8 buffer. If N = &String or &str, we don't own the string and drop_in_place is a noop.
        // Comment on speed of using MaybeUninit and drop_in_place:
        // - For "hot" data where threads are not blocked, it actually becomes important that
        //   args/seqs be consumed/initialized in parallel, and NOT sequentially before/after
        //   the worker loop. For example, we could initialize the seqs array at the start:
        //      let mut seqs = Vec::new();
        //      seqs.resize(args.len(), String::new());
        //   But cloning all those strings sequentially does affect performance, surprisingly.
        // - Similarly, we could just let args and all its contents be automatically dropped
        //   outside the worker loop, as get_batch returns, but the sequential String deallocs
        //   turn out to be slow, so we do  drop_in_place inside the parallel workers instead.
        parallel_for(zip(arg_chunks, seq_chunks), |(args, seqs)| {
            for ((name, start, end), seq) in zip(args, seqs) {
                seq.write(self.get(name.as_ref(), *start, *end));
                unsafe { ptr::drop_in_place(name) };  // Drop here emulates an "Into" version of ChunksMut
            }
        });

        // Mark outputs as initialized. Mark inputs as consumed (already dropped in-place).
        unsafe {
            seqs.set_len(args.len());
            args.set_len(0);
        }

        seqs
    }
    
    /// A version of [`get`](Self::get) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let seq = tbr.get("chr2", 10001, 10010);     // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    /// 
    /// # Panics
    ///
    /// Panics if the sequence name was not found or the range was invalid.
    ///
    pub fn get_inclusive<N: AsRef<str>>(&self, name: N, start: usize, end: usize) -> String {
        self[name.as_ref()].get_inclusive(start, end)
    }

    /// A version of [`get_into`](Self::get_into) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    /// 
    /// # Examples
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
        self[name.as_ref()].get_inclusive_into(start, end, dst);
    }

    /// A version of [`get_batch`](Self::get_batch) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let args = [("chr2", 10001, 10006),
    ///             ("chr2", 10007, 10010), /* ... */];
    /// let seqs = tbr.get_inclusive_batch(args);  // ["CGTATC", "CCAC", ...]
    /// # Ok::<(), io::Error>(())
    /// ```
    /// 
    /// Internally this method is very similar to using `rayon` (see below), except a
    /// global thread pool is not maintained.
    /// ```no_run
    /// # use rayon::prelude::*;
    /// # use twobitreader::TwobitReader;
    /// # let tbr = TwobitReader::open("hg38.2bit").unwrap();
    /// # let args = [("chr2", 10001, 10006),
    /// #             ("chr2", 10007, 10010)];
    /// let seqs = args.into_par_iter().map(
    ///     |(name, start, end)| tbr.get_inclusive(name, start, end)
    /// ).collect::<Vec<_>>();  // Vec<String>
    /// ```
    /// 
    /// # Panics
    ///
    /// Panics if any sequence name was not found or if any range was invalid.
    /// 
    pub fn get_inclusive_batch<N, I>(&self, iter: I) -> Vec<String>
    where
        N: AsRef<str> + Sync + Send,   // S = String, &String, &str, etc...
        I: IntoIterator<Item = (N, usize, usize)>,
    {
        self.get_batch(iter.into_iter().map(|(chrom, start, end)| {
            check_start_inclusive(start);
            (chrom, start-1, end)
        }))
    }

    /// Concatenates a batch of sequence ranges into a single string.
    /// 
    /// # Examples
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
        self[name.as_ref()].concat(ranges)
    }

    /// A version of [`concat`](Self::concat) using 1-based inclusive ranges;
    /// see [genomic interval notations]((https://standage.github.io/on-genomic-interval-notation.html)).
    /// 
    /// # Examples
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
        self[name.as_ref()].concat_inclusive(ranges)
    }
}

impl<Q: ?Sized> Index<&Q> for TwobitReader
where
    Q: Eq + Hash,
    String: Borrow<Q>,
{
    type Output = TwobitSequence;

    /// Returns a reference to the [`TwobitSequence`] for the named sequence record.
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let tbs = &tbr["chr2"];           // &TwobitSequence
    /// let seq = tbs.get(10000, 10010);  // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    /// 
    /// # Panics
    ///
    /// Panics if no sequence record has the given name.
    /// 
    fn index(&self, name: &Q) -> &TwobitSequence {
        let index = *self.seq_by_name.get(name).expect("TwobitReader: sequence name not found");
        &self.seqs[index]
    }
}

impl Index<String> for TwobitReader
 {
     type Output = TwobitSequence;

    /// Returns a reference to the [`TwobitSequence`] for the named sequence record.
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// let tbs = &tbr["chr2".to_string()];   // &TwobitSequence
    /// let seq = tbs.get(10000, 10010);      // "CGTATCCCAC"
    /// # Ok::<(), io::Error>(())
    /// ```
    /// 
    /// # Panics
    ///
    /// Panics if no sequence record has the given name.
    /// 
    fn index(&self, name: String) -> &TwobitSequence {
        let index = *self.seq_by_name.get(&name).expect("TwobitReader: sequence name not found");
         &self.seqs[index]
     }
}

impl<'a> IntoIterator for &'a TwobitReader {
    type Item = &'a TwobitSequence;
    type IntoIter = <&'a Vec<TwobitSequence> as IntoIterator>::IntoIter;

    /// Allows a TwobitReader reference to be turned into an iterator over its [`TwobitSequence`] records.
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use std::io;
    /// # use twobitreader::TwobitReader;
    /// let tbr = TwobitReader::open("hg38.2bit")?;
    /// for tbs in &tbr {
    ///     println!("{} {}", tbs.name(), tbs.len());
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    fn into_iter(self) -> Self::IntoIter {
        self.seqs.iter()
    }
}

impl TwobitSequence {

    /// Returns the name of this sequence record.
    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    /// Returns the length of this sequence record, in nucleotides.
    pub fn len(&self) -> usize {
        self.dna_len
    }

    /// See [`TwobitReader::get`].
    pub fn get(&self, start: usize, end: usize) -> String {
        let mut dst = String::new();
        self.get_into(start, end, &mut dst);
        dst
    }

    /// See [`TwobitReader::get_into`].
    pub fn get_into(&self, start: usize, end: usize, dst: &mut String) {
        // Check range before trying to do any might-panic arithmetic with (start, end)
        self.check_range(start, end);
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
            self.get_into_u8(start, dst);
        }
    }

    /// See [`TwobitReader::get_inclusive`].
    pub fn get_inclusive(&self, start: usize, end: usize) -> String {
        // Check valid start and then convert range to 0-based exclusive.
        check_start_inclusive(start);
        self.get(start-1, end)
    }

    /// See [`TwobitReader::get_inclusive_into`].
    pub fn get_inclusive_into(&self, start: usize, end: usize, dst: &mut String) {
        // Check valid start and then convert range to 0-based exclusive.
        check_start_inclusive(start);
        self.get_into(start-1, end, dst);
    }

    /// See [`TwobitReader::concat`].
    pub fn concat<R>(&self, ranges: R) -> String
    where
        R: IntoIterator<Item = (usize, usize)>,
        R::IntoIter: Clone,
    {
        // Compute total length by summing individual interval lengths. Also check the ranges here,
        // since get_into_u8 does not produce user-friendly panic messages for invalid ranges.
        let ranges = ranges.into_iter();
        let total_len = ranges.clone().map(|(start, end)| { self.check_range(start, end); end - start } ).sum();
        
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
                self.get_into_u8(start, &mut dst_vec[offset..offset+len]);
                offset += len;
            }
        }
        dst
    }

    /// See [`TwobitReader::concat_inclusive`].
    pub fn concat_inclusive<R>(&self, ranges: R) -> String
    where
        R: IntoIterator<Item = (usize, usize)>,
        R::IntoIter: Clone,
    {
        // Check each start, and convert range to 0-based exclusive.
        self.concat(ranges.into_iter().map(|(start, end)| {
            check_start_inclusive(start);
            (start-1, end)
        }))
    }

    // Decodes sequence [start..end] into dst, where end = start + dst.len().
    // All bytes in dst are replaced. The first byte in slice is assumed
    fn get_into_u8(&self, start: usize, dst: &mut [u8]) {
        // Slice spanning all packed 2bit DNA data for this sequence record.
        let dna = unsafe { slice::from_raw_parts(self.dna_ptr, self.dna_bytes) };

        // Decode packed 2-bit from dna into dst from the given start position.
        decode(start, dna, dst);

        // Apply N-block and lowercase masks to dst in-place.
        search_blocks(&self.nblocks, start, dst, nblock_fill);
        search_blocks(&self.masks, start, dst, mask_fill);
    }

    // Panics if (start, end) is an invalid 0-based exclusive range for this TwobitSequence.
    // (The error message does not format the actual start/end values that were passed in.
    //  This is because the 0-based get functions are used to implement the 1-based get
    //  functions, and so the start/end seen here may not be the actual ones users provided,
    //  which would result in more confusion for the user than simply omitting the numbers.)
    #[inline]
    fn check_range(&self, start: usize, end: usize) {
        if start > end { panic!("TwobitReader: invalid range (start > end)"); }
        if end > self.dna_len { panic!("TwobitReader: invalid end (end > dna_len)"); }
    }
}

// Panics ifs start is not valid for an inclusive range
#[inline]
fn check_start_inclusive(start: usize) {
    if start == 0 { panic!("TwobitReader: invalid start (start = 0)"); }
}

// Reads the header and sequence data from a .2bit file.
// The resulting vec of TwobitSequences is intended to be assigned to the TwobitReader::seqs field.
// Block intervals are accessed (paged in from disk), but the sequence data itself is not paged in.
fn read_seqs(mmap: &Mmap, use_mask: bool) -> io::Result<Vec::<TwobitSequence>> {
    // Twobit file signature
    const SIGNATURE_LILEND: u32 = 0x1A412743;  // Twobit signature as little-endian u32 (native)
    const SIGNATURE_BIGEND: u32 = 0x4327411A;  // Twobit signature as big-endian u32 (non-native)

    // Start reading the file, using appropriate byte order
    let signature = Cursor::new(mmap).read_u32::<LittleEndian>()?;
    return match signature {
        SIGNATURE_LILEND => read_seqs_from_endian::<LittleEndian>(mmap, use_mask),
        SIGNATURE_BIGEND => read_seqs_from_endian::<BigEndian>(mmap, use_mask),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "Twobit file signature not found.")),
    }
}

// Implements read_seqs after the endian-ness of the file has been determined from the signature.
fn read_seqs_from_endian<B: ByteOrder>(mmap: &Mmap, use_mask: bool) -> io::Result<Vec::<TwobitSequence>> {
    // Start cursor immediately after the signature field.
    let mut cursor = Cursor::new(mmap);
    cursor.set_position(size_of::<u32>() as u64);

    // Check the file version
    let version = cursor.read_u32::<B>()?;
    if version != 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
                    "TwobitReader: file version not recognized."));
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
            Err(_) => return Err(io::Error::new(io::ErrorKind::InvalidData,
                                 "TwobitReader: failed reading sequence name as utf8.")),
        };

        // Get the offset of the data arrays for this sequence 
        let data_offset = cursor.read_u32::<B>()? as usize;
        let curr_offset = cursor.position() as usize;
        if data_offset <= curr_offset {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                       "TwobitReader: invalid data offset; file may be malformed."));
        }

        seq_headers.push((name, data_offset));
    }

    // Create a pre-sized result vec, with all slots initially filled with None.
    // (Can't use vec![] or resize; TwobitSequence doesn't implement Clone.)
    let mut seqs = Vec::from_iter(repeat_with(|| None).take(num_seqs));

    // Read each sequence's data block, in parallel. Reading in parallel is very important
    // for the speed of opening a "cold" file that is not yet in the page cache.
    // sequence_header[i] will be used to read and create a TwobitSequence instance, which
    // will then be assigned to seqs[i] (replacing the None).
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

        // Create (unsafe) pointer to dna bytes, if entire span was mapped.
        let dna_offset = cursor.position() as usize;
        let file_len = cursor.get_ref().len();
        if dna_offset + (dna_len + NUCS_PER_U8 - 1) / NUCS_PER_U8 > file_len {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof,
                "TwobitReader: dna bytes were truncated in file."));
        }
        let dna_ptr = unsafe { cursor.get_ref().as_ptr().add(dna_offset) };
        let dna_bytes = (dna_len + NUCS_PER_U8 - 1) / NUCS_PER_U8;

        // Initialize the sequence.
        *seq = Some(TwobitSequence{dna_ptr, dna_bytes, dna_len, nblocks, masks, name});
        Ok(())
    })?;

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
            *end = *start + *end;
        }

    } else {
        // Advance the cursor, without reading the data itself.
        let new_offset = cursor.position() as usize + 2 * num_blocks * size_of::<u32>();
        if new_offset > cursor.get_ref().len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof,
                "TwobitReader: block indices were truncated in file."));
        }
        cursor.set_position(new_offset as u64);
    }

    Ok(Blocks{starts, ends})
}

#[cfg(debug_assertions)]
fn check_blocks(blocks: &Blocks, dna_len: usize) -> io::Result<()> {
    if !is_sorted(&blocks.starts) || !is_sorted(&blocks.ends) {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            "TwobitReader: block starts / ends were not strictly ascending; file may be malformed."));
    }

    for (start, end) in zip(&blocks.starts, &blocks.ends) {
        if start >= end {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "TwobitReader: expected all block start < block end; file may be malformed."));
        }
        if *start as usize > dna_len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "TwobitReader: block start was beyond sequence length; file may be malformed."));
        }
        if *end as usize > dna_len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                "TwobitReader: block end was beyond sequence length; file may be malformed."));
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
        let block_end   = (*block_end as usize).min(end);
        if block_start >= block_end { break; }
        f(&mut dst[block_start-start..block_end-start]);
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
