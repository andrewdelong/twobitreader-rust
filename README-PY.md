# twobitreader-rs

Fast DNA sequence extraction from [2bit files](http://genome.ucsc.edu/FAQ/FAQformat.html#format7),
the compact genome format used by the UCSC Genome Browser.

The motivation for this package is speed. See the underlying [`twobitreader`](https://crates.io/crates/twobitreader) Rust crate for benchmarks.

```console
pip install twobitreader-rs
```

## Usage

```python
from twobitreader_rs import TwobitReader

tbr = TwobitReader.open("hg38.2bit")
tbr.get("chr1", 10000, 10005)        # 'TAACC'
tbr.seq_len("chr1")                  # 248956422
tbr.names()                          # ['chr1', 'chr2', ...]
```

Ranges are 0-based with an exclusive end, like a Python `start:end` slice.
Every method also has an `_inclusive` counterpart taking 1-based inclusive ranges,
as used by GFF/GTF and genome browsers.

**Concatenation** helps to assemble a spliced transcript from its exons:

```python
from twobitreader_rs import reverse_complement

tx = tbr.concat("chr6", [(1389575, 1391118), (1394695, 1395603)])
tx = reverse_complement(tx)          # if the transcript is on the minus strand
```

**Prefetching** dramatically speeds up access to "cold" files that are not yet in
the page cache:

```python
tbr.prefetch(exons)                  # tell the OS what is coming
seqs = tbr.get_batch(exons)          # reads land as the data arrives
```


## Documentation

The Python methods closely mirror their counterparts in the [Rust crate's documentation](https://docs.rs/twobitreader).
Every Python method carries a full docstring with arguments, return values and exceptions, visible
through `help(TwobitReader)`.

## License

MIT OR Apache-2.0
