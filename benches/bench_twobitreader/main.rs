use twobitreader::TwobitReader;

// Standard library
use std::collections::HashMap;
use std::error::Error;
use std::hint::black_box;
use std::time::{Duration, Instant};

// Modules
mod util;
use util::*;

// Dependencies
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

// Load hg38 exon intervals and time how long it takes to extract their DNA.
fn bench_hg38_exons(use_rayon: bool) -> Result<Duration, Box<dyn Error>> {
    let path = download_hg38()?;
    let exons = read_exons()?;
    let tic = Instant::now();

    // Open file and collect strings into a vec, either in parallel or sequentially.
    let tbr = TwobitReader::open(path)?;
    let dst: Vec<_> = if use_rayon {
        exons.par_iter().map(|(chrom, start, end)| tbr.get(chrom, *start, *end)).collect()
    } else {
        exons.iter().map(|(chrom, start, end)| tbr.get(chrom, *start, *end)).collect()
    };

    let duration = tic.elapsed();
    black_box(dst);
    Ok(duration)
}

// Load hg38 transcript intervals (groups of exons) from a BED file and then time how
// long it takes to extract DNA sequences for all those transcripts.
fn bench_hg38_transcripts(use_rayon: bool) -> Result<Duration, Box<dyn Error>> {
    let path = download_hg38()?;
    let transcripts = read_transcripts()?;
    let tic = Instant::now();

    // Open file and collect transcript strings into a hashmap, either in parallel or sequentially.
    let tbr = TwobitReader::open(path)?;
    let dst: HashMap<_, _> = if use_rayon {
        transcripts.par_iter().map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons.iter().copied()))).collect()
    } else {
        transcripts.iter().map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons.iter().copied()))).collect()
    };

    let duration = tic.elapsed();
    black_box(dst);
    Ok(duration)
}

fn main() -> Result<(), Box<dyn Error>> {
    ThreadPoolBuilder::new().build_global().unwrap();
    load_into_page_cache(download_hg38()?.as_path())?;
    bench("hg38_exons_basic", || bench_hg38_exons(false), 5, "ms");
    bench("hg38_exons_rayon", || bench_hg38_exons(true), 5, "ms");
    bench("hg38_transcripts_basic", || bench_hg38_transcripts(false), 5, "ms");
    bench("hg38_transcripts_rayon", || bench_hg38_transcripts(true), 5, "ms");
    Ok(())
}
