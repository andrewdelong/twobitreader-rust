# Details of speed comparison

The table of timings shown in the README is described below.

### Libraries
* [`twobitreader 0.1`](https://github.com/andrewdelong/twobitreader-rust) - Rust crate for reading 2bit files
* [`twobit 0.2`](https://github.com/jbethune/rust-twobit) - Rust crate for reading and writing 2bit files
* [`py2bit 0.3`](https://github.com/deeptools/py2bit) - Python package written in C
* [`GenomeKit 7.6.1`](https://github.com/deepgenomics/GenomeKit) - Python package written in C++
* [`twobitreader 3.1`](https://github.com/benjschiller/twobitreader) - Python package
* [`twobitToFa`](https://genome.ucsc.edu/goldenPath/help/twoBit.html) - Command line utility written in C

### Data
* All experiments extracted from `hg38.p13.2bit` (850MB, [link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/)).
* All experiments read intervals from this crate's `gencode-exons.bed` (133,388 lines) or 
  `gencode-transcripts.bed` (319,468 lines) located in `benches/assets`.
  See also [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
* Exon and transcript intervals were extracted from a subset of [Gencode v43](https://www.gencodegenes.org/human/)
  basic gene annotations; specifically, those transcripts having `transcript_type=protein_coding`
  and `transcript_support_level=1`.

### System
* MacBook Pro running macOS 26.7
* 2.3 GHz 8-core Intel i9
* SSD hard drive
* Rust 1.98

### Details
* Lowercase masks were not used during extraction. Only N-blocks were enabled.
* Ranges were treated as positive-strand ("+"), for consistency across libraries.
* Times are average of 10 runs, rounded to two digits.
* Time for opening the `.2bit` file was included for all methods.
* Time for reading the `.bed` file was not included, except for command-line `twobitToFa` (`bed.gz` was pre-unzipped).
* All parallel experiments used 16 threads, except `twobitToFa` which was faster with 8.
* Missing timings (`n/a`):
  * Parallel timings for `py2bit` were not collected because Python's `multiprocessing` requires pickling
    (`TypeError: cannot pickle 'py2bit.pyTwoBit' object`).
  * Parallel timings for `twobit` were not collected because it requires a mutable reference
    that cannot be shared across threads by parallel libraries such as `rayon`.
* Hot vs cold:
  * In "hot" experiments, the entire 2bit file was read into page cache beforehand.
  * In "cold" experiments, cache was disabled via `fcntl(F_NOCACHE)` (twobitreader) or wiped via `sudo purge` (others).


### Extra details for Rust benchmarking

See `benches/bench_twobitreader.rs` in this crate. Nearly identical code was used to benchmark `twobit`.

### Extra details for Python benchmarking

A simplified version of the benchmarking for `py2bit` and `twobit` is shown below.
```python
def bench_exons_sequential_py2bit():
    exons = read_exons()
    start_time = time()
    
    tb = py2bit.open("hg38.p13.2bit")
    seqs = [tb.sequence(*exon) for exon in exons]
    
    return time() - start_time
```

Parallel experiments were the same, but used `multiprocessing.Pool` to implement a parallel for loop.

### Extra details for `twobitToFa`
* The command was `twobitToFa -bed={input}.bed -noMask hg38.p13.2bit /dev/null`.
* Parallel runs used `xargs -P8`, where 8 processes was faster than defaults.
* Parallel inputs were pre-split `.bed` files (10,000 line chunks). Split time was not included.
  * Strictly speaking, splitting `gencode-transcripts.bed` this way is incorrect, 
    because it stores one exon per line (does not use the "block" columns), but 
    naive splitting is still representative of speed.

### Extra details for `GenomeKit`

* Timings include constructing the `Interval(chrom, "+", start, end, "hg38.p13")` 
  object for each dna query, since that's a GenomeKit-specific overhead.
* The strand was kept as '+' even though GenomeKit is capable of 
  strand-sensitive extraction. Reverse complementing has a small overhead, 
  so this was for comparability with the more low-level libraries.
* The parallel timings say "n/a" because GenomeKit is not yet compatible
  with free-threaded Python, and the pickling overhead of `multiprocessing.Pool` negates any performance benefit of parallelism this way.