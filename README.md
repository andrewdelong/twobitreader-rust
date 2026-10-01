# twobitreader


This crate provides fast DNA sequence extraction from 2bit files, a
[standard format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) in bioinformatics.

The motivation for `twobitreader` is speed; see benchmarks below.
It is also available as a Python package `twobitreader_rs`.

The focus is raw reading from 2bit, but fast concatenation and reverse-complement methods
are also provided to make higher-level use cases easier.

[![CI](https://github.com/andrewdelong/twobitreader-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/andrewdelong/twobitreader-rust/actions/workflows/ci.yml)
&nbsp;
![Windows](https://img.shields.io/badge/Windows-supported-blue?logo=windows)
![macOS](https://img.shields.io/badge/macOS-supported-blue?logo=apple)
![Linux](https://img.shields.io/badge/Linux-supported-blue?logo=linux)

## Examples

**Extracting sequences** is straightforward:
```rust
let tbr = TwobitReader::open("hg38.2bit")?; // Human genome, build 38
let seq = tbr.get("chr1", 10000, 10005);    // -> String ("TAACC")
```

**Concatenation** works by iterating over (start, end) pairs.
For example, assembling a spliced transcript:
```rust
// Exon ranges for human FKHL6 gene transcript (Gencode v43)
let exons = [(1389575, 1391118),             // Exon 1 (start, end)
             (1394695, 1395603)];            // Exon 2 (start, end)
let transcript = tbr.concat("chr6", &exons); // -> String
```

**Parallelism** is easy with crates like [`rayon`](https://docs.rs/rayon/latest/rayon/).
For example, batch extraction of sequences:
```rust
use rayon::prelude::*;
let args = [("chr1", 10000, 15000),
            ("chr1", 30000, 35000), /* ... */ ];
let seqs = args.into_par_iter()
    .map(|(chrom, start, end)| tbr.get(chrom, start, end))
    .collect::<Vec<_>>(); // -> Vec<String>
```
Or, assembling a batch of spliced transcripts in parallel:
```rust
use twobitreader::reverse_complement;
use rayon::prelude::*;

fn stranded(seq: String, strand: char) -> String {
    if strand == '+' { seq } else { reverse_complement(seq) }
}

let transcripts = [                   // (transcript_id, chromosome, exons)
    ("ENST00000407983.7", "chr2", '+', vec![(264899, 265007),     // Exon 1
                                            (271865, 271939),     // Exon 2
                                            (272036, 272557)]),   // Exon 3
    ("ENST00000319331.4", "chr3", '+', vec![(3799430, 3799919),   // Exon 1
                                            (3844363, 3849834)]), // Exon 2
    /* ... */
];
// Concatenate exons and reverse-complement if negative strand.
// (Correct when exons are listed in genome-coordinate order.)
let seqs = transcripts.into_par_iter()
    .map(|(id, chrom, strand, exons)| (id, stranded(tbr.concat(chrom, exons), strand)))
    .collect::<HashMap<_, _>>();      // HashMap<&str, String>
let seq = &seqs["ENST00000407983.7"]; // -> &String to transcript sequence
```

**Cold files** are an order of magnitude slower to access than files already in memory ("hot").
Use prefetching to dramatically improve single-threaded speed:
```rust
let exons = [("chr1", 10000, 10200),
             ("chr1", 10500, 10700), /* ... */ ];
tbr.prefetch(&exons);             // Ask the operating system to start paging this data from disk.
let seqs = tbr.get_batch(&exons); // Access the memory as it arrives.
```

## Benchmarks

Two tasks were benchmarked:
- **exons**: extract 133,388 distinct human exon sequences;
- **transcripts**: concatenate 319,468 exons into 29,211 human spliced transcript sequences.

Speed depends on parallelism and page cache (hot vs cold):
- **hot** runs represent repeated or interactive dna extraction scenarios;
- **cold** runs represent a first run of a genomics pipeline, bound by disk speed;
- **prefetch** runs are cold but with a `prefetch` call preceding extraction.

The table below shows running times in milliseconds. Experimental details are `BENCH.md`.
This crate provides *twobitreader* (pure rust) and *twobitreader_rs* (python wrapper).

| EXONS | 1-thread / hot | 1-thread / cold | 1-thread / prefetch | 16-thread / hot | 16-thread / cold |
|---|---:|---:|---:|---:|---:|
| **twobitreader** (rs)        | 90    | 1,700  | 180 | 10    | 190   |
| **twobitreader_rs** (rs, py) | 100   | 1,900  | 190 | *27   | *210  |
| **py2bit** (c, py)           | 220   | 2,800  | n/a | n/a   | n/a   |
| **GenomeKit** (cpp, py)      | 340   | 2,400  | n/a | n/a   | n/a   |
| **twobit** (rs)              | 390   | 3,500  | n/a | n/a   | n/a   |
| **twobitToFa** (c)           | 1,000 | 4,500  | n/a | 340   | 850   |
| **twobitreader** (py)        | 7,200 | 13,000 | n/a | 1,900 | 2,300 |
| **Biopython** (py)           | 8,100 | 11,000 | n/a | n/a   | n/a   |

| TRANSCRIPTS | 1-thread / hot | 1-thread / cold | 1-thread / prefetch | 16-thread / hot | 16-thread / cold |
|---|---:|---:|---:|---:|---:|
| **twobitreader** (rs)        | 140    | 1,700  | 230 | 13    | 190   |
| **twobitreader_rs** (rs, py) | 180    | 2,100  | 310 | *27   | *220  |
| **py2bit** (c, py)           | 490    | 3,600  | n/a | n/a   | n/a   |
| **GenomeKit** (cpp, py)      | 720    | 2,900  | n/a | n/a   | n/a   |
| **twobit** (rs)              | 930    | 1,900  | n/a | n/a   | n/a   |
| **twobitToFa** (c)           | 2,100  | 6,900  | n/a | 840   | 1,200 |
| **twobitreader** (py)        | 13,000 | 21,000 | n/a | 2,800 | 3,500 |
| **Biopython** (py)           | 19,000 | 25,000 | n/a | n/a   | n/a   |

Entries marked * were run in free-threaded Python.

## Dependencies

* `byteorder` for handling endian-ness
* `memmap2` for memory mapping the 2bit file
* `seq-macro` for generating 2bit decoder lookup table
* `libc` for prefetching file ranges on Apple targets
* `windows-sys` for prefetching file ranges on Windows targets

## Use of AI

Claude Code: generated the OS-specific prefetch loops; improved handling of corrupt or malicious files;
documented the Python bindings and mirrored their tests and benchmark; improved error checking and
propagation to Python more broadly; and generated the CI configurations.

License: MIT OR Apache-2.0
