# twobitreader

This crate provides fast DNA sequence extraction from 2bit files, a
[standard format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) in bioinformatics.

The motivation for this crate is speed.
Extracting sequences is 1.5-30x faster than the best alternative, depending on use case.

## Examples

**Extracting sequences** is straightforward with [`TwobitReader`]:
```rust
let tbr = TwobitReader::open("hg38.2bit")?;  // Human genome, build 38
let seq = tbr.get("chr1", 10000, 10005);     // String ("TAACC")
```

**Batch extraction** is fast.
It works by iterating over (chrom, start, end) triplets in parallel:
```rust
let args = [("chr1", 10000, 15000),
            ("chr1", 30000, 35000), /* ... */ ];
let seqs = tbr.get_batch(args);     // Vec<String>
```

**Concatenation** works by iterating over (start, end) pairs.
For example, assembling a spliced transcript:
```rust
// Exon ranges for human FKHL6 gene transcript (Gencode v43)
let exons = [(1389575, 1391118),             // Exon 1 (start, end)
             (1394695, 1395603)];            // Exon 2 (start, end)
let transcript = tbr.concat("chr6", exons);  // String
```

**Parallelism** is easy with crates like [`rayon`](https://docs.rs/rayon/latest/rayon/).
For example, assembling a batch of spliced transcripts in parallel:
```rust
use rayon::prelude::*;
let transcripts = [     // (transcript_id, chromosome, exons)
    ("ENST00000407983.7", "chr2", vec![(264899, 265007),      // Exon 1 (start, end)
                                       (271865, 271939),      // Exon 2 (start, end)
                                       (272036, 272557)]),    // Exon 3 (start, end)
    ("ENST00000319331.4", "chr3", vec![(3799430, 3799919),    // Exon 1 (start, end)
                                       (3844363, 3849834)]),  // Exon 2 (start, end)
    /* ... */
];
let seqs = transcripts.into_par_iter()
                      .map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons)))
                      .collect::<HashMap<_, _>>();  // HashMap<&str, String>
let seq = &seqs["ENST00000407983.7"];               // Look up transcript sequence
```

## Speed

Two tasks were benchmarked:
- **exons**: extract 133,388 distinct human exon sequences;
- **transcripts**: concatenate 319,468 exons into 29,180 human spliced transcript sequences.

Speed depends on parallelism and page cache (hot vs cold):
- **hot** represents computations that are repeated, or are run interactively, on the same server;
- **cold** represents a first run of a computational pipeline, and is bound by disk speed.

The table below shows running times in milliseconds. Experimental details are `doc/comparison-details.md`.
<table style="font-size:.85em;">
<thead>
<tr>
    <th rowspan=2 colspan=2 style="border:none; vertical-align:bottom;"></th>
    <th colspan=5>file in page cache (hot)</th>
    <th style="border:none; padding:4pt;"></th>
    <th colspan=5>file <i>not</i> in page cache (cold)</th>
</tr>
<tr>
    <th colspan=2>sequential</th>
    <th style="border:none; padding:4pt;"></th>
    <th colspan=2>parallel</th>
    <th style="border:none; padding:4pt;"></th>
    <th colspan=2>sequential</th>
    <th style="border:none; padding:4pt;"></th>
    <th colspan=2>parallel</th>
</tr>
<tr>
    <th style="border:none;"></th>
    <th style="border:none; padding:4pt;"></th>
    <th>exons</th>
    <th>transcripts</th>
    <th style="border:none; padding:4pt;"></th>
    <th>exons</th>
    <th>transcripts</th>
    <th style="border:none; padding:4pt;"></th>
    <th>exons</th>
    <th>transcripts</th>
    <th style="border:none; padding:4pt;"></th>
    <th>exons</th>
    <th>transcripts</th>
</tr>
</thead>
<tbody>
<tr>
    <td style="font-weight:bold;">twobitreader (rust)</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">92</td>
    <td style="text-align:right">160</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right">12</td>
    <td style="text-align:right">26</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">2,200</td>
    <td style="text-align:right">2,600</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right">210</td>
    <td style="text-align:right">220</td>
</tr>
<tr>
    <td style="font-weight:bold">py2bit (C, python)</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">180</td>
    <td style="text-align:right">400</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right;color:#888">n/a</td>
    <td style="text-align:right;color:#888">n/a</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">3,100</td>
    <td style="text-align:right">4,000</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right;color:#888">n/a</td>
    <td style="text-align:right;color:#888">n/a</td>
</tr>
<tr>
    <td style="font-weight:bold">twobit (rust)</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">390</td>
    <td style="text-align:right">890</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right;color:#888">n/a</td>
    <td style="text-align:right;color:#888">n/a</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">3,500</td>
    <td style="text-align:right">7,700</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right;color:#888">n/a</td>
    <td style="text-align:right;color:#888">n/a</td>
</tr>
<tr>
    <td style="font-weight:bold">twobitToFa (C)</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">1,000</td>
    <td style="text-align:right">2,100</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right">340</td>
    <td style="text-align:right">840</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">4,500</td>
    <td style="text-align:right">6,900</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right">850</td>
    <td style="text-align:right">1,200</td>
</tr>
<tr>
    <td style="font-weight:bold">twobitreader (python)</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">7,200</td>
    <td style="text-align:right">13,000</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right">1,900</td>
    <td style="text-align:right">2,800</td>
    <td style="border:none; padding:4pt;"></td>
    <td style="text-align:right">13,000</td>
    <td style="text-align:right">21,000</td>
    <th style="border:none; padding:4pt;"></th>
    <td style="text-align:right">2,300</td>
    <td style="text-align:right">3,500</td>
</tr>
</tbody>
</table>

## Dependencies

* `byteorder` for handling endian-ness
* `memmap2` for memory mapping the 2bit file
* `num_cpus` for spawning parallel workers
* `seq-macro` for generating 2bit decoder lookup table
