use twobitreader::TwobitReader;

// Standard library
use std::collections::HashMap;
use std::error::Error;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::{Duration, Instant};

// Modules
mod util;
use util::*;

// Dependencies
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

// Which cache scenario to simulate in the benchmark run.
#[derive(Clone, Copy, PartialEq)]
enum Cache {
    Hot,          // File already resides in the page cache.
    Cold,         // File is not yet in the page cache at all.
    ColdPrefetch, // Cold, but with a prefetch hint at upcoming data access.
}

// How and whether to use parallelism in the benchmark run.
#[derive(Clone, Copy, PartialEq)]
enum Parallel {
    None,          // No parallelism.
    Rayon,         // Use Rayon for parallelism.
}

// Returns the path to benchmark against, preparing the page cache as the scenario requires.
// Setup happens outside the timed region, so its cost is not part of any reported duration.
fn setup(cache: Cache) -> Result<PathBuf, Box<dyn Error>> {
    let path = download_hg38()?;
    match cache {
        Cache::Hot => {
            load_into_page_cache(&path)?;
            Ok(path)
        }
        Cache::Cold | Cache::ColdPrefetch => Ok(make_cold_copy(&path)?),
    }
}

// Load hg38 and time processing the header
fn bench_hg38_open(cache: Cache) -> Result<Duration, Box<dyn Error>> {
    let path = setup(cache)?;
    let tic = Instant::now();

    // Open file, either the hot version or the cold version
    let tbr = TwobitReader::open(path)?;

    let duration = tic.elapsed();
    black_box(tbr);
    Ok(duration)
}

// Load hg38 exon intervals and time how long prefetch takes, without waiting for data to finish paging in.
// The goal is to estimate how much time in the ColdPrefetch cases is taken by the prefetch phase itself.
fn bench_hg38_exons_prefetch_only() -> Result<Duration, Box<dyn Error>> {
    let path = setup(Cache::Cold)?;
    let exons = read_exons()?;
    let tic = Instant::now();

    // Open file and collect strings into a vec, either in parallel or sequentially.
    let tbr = TwobitReader::open(path)?;
    tbr.prefetch(&exons);

    let duration = tic.elapsed();
    Ok(duration)
}

// Load hg38 exon intervals and time how long it takes to extract their DNA.
fn bench_hg38_exons(cache: Cache, parallel: Parallel) -> Result<Duration, Box<dyn Error>> {
    let path = setup(cache)?;
    let exons = read_exons()?;
    let tic = Instant::now();

    // Open file and collect strings into a vec, either in parallel or sequentially.
    let tbr = TwobitReader::open(path)?;
    if cache == Cache::ColdPrefetch {
        tbr.prefetch(&exons);
    }
    let dst: Vec<_> = if parallel == Parallel::Rayon {
        exons.par_iter().map(|(chrom, start, end)| tbr.get(chrom, *start, *end)).collect()
    } else {
        tbr.get_batch(&exons).collect()
    };

    let duration = tic.elapsed();
    black_box(dst);
    Ok(duration)
}

// Load hg38 transcript intervals (groups of exons) from a BED file and then time how
// long it takes to extract DNA sequences for all those transcripts.
fn bench_hg38_transcripts(cache: Cache, parallel: Parallel) -> Result<Duration, Box<dyn Error>> {
    let path = setup(cache)?;
    let transcripts = read_transcripts()?;
    let tic = Instant::now();

    // Open file and collect transcript strings into a hashmap, either in parallel or sequentially.
    let tbr = TwobitReader::open(path)?;
    if cache == Cache::ColdPrefetch {
        // Iterator across the exons of all transcripts.
        let exons = transcripts.iter()
            .flat_map(|(_id, chrom, exons)| {
                exons.iter().map(move |&(start, end)| (chrom, start, end))
            }
        );
        tbr.prefetch(exons);
    }
    let dst: HashMap<_, _> = if parallel == Parallel::Rayon {
        transcripts.par_iter().map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons))).collect()
    } else {
        transcripts.iter().map(|(id, chrom, exons)| (id, tbr.concat(chrom, exons))).collect()
    };

    let duration = tic.elapsed();
    black_box(dst);
    Ok(duration)
}

fn main() -> Result<(), Box<dyn Error>> {
    ThreadPoolBuilder::new().build_global().unwrap();

    // Hot: the file is already in the page cache, so these measure decoding speed.
    bench("hg38_open_hot", || bench_hg38_open(Cache::Hot), 3, "ms");
    bench("hg38_exons_hot_basic", || bench_hg38_exons(Cache::Hot, Parallel::None), 10, "ms");
    bench("hg38_exons_hot_rayon", || bench_hg38_exons(Cache::Hot, Parallel::Rayon), 10, "ms");
    bench("hg38_transcripts_hot_basic", || bench_hg38_transcripts(Cache::Hot, Parallel::None), 10, "ms");
    bench("hg38_transcripts_hot_rayon", || bench_hg38_transcripts(Cache::Hot, Parallel::Rayon), 10, "ms");

    // Cold: the file must be read from disk, so these measure how well the reads are overlapped.
    if make_cold_copy(download_hg38()?).is_ok() {
        bench("hg38_open_cold", || bench_hg38_open(Cache::Cold), 3, "ms");
        bench("hg38_exons_cold_prefetch_only", || bench_hg38_exons_prefetch_only(), 10, "ms");
        bench("hg38_exons_cold_prefetch", || bench_hg38_exons(Cache::ColdPrefetch, Parallel::None), 10, "ms");
        bench("hg38_exons_cold_basic", || bench_hg38_exons(Cache::Cold, Parallel::None), 10, "ms");
        bench("hg38_exons_cold_rayon", || bench_hg38_exons(Cache::Cold, Parallel::Rayon), 10, "ms");
        bench("hg38_transcripts_cold_prefetch", || bench_hg38_transcripts(Cache::ColdPrefetch, Parallel::None), 10, "ms");
        bench("hg38_transcripts_cold_basic", || bench_hg38_transcripts(Cache::Cold, Parallel::None), 10, "ms");
        bench("hg38_transcripts_cold_rayon", || bench_hg38_transcripts(Cache::Cold, Parallel::Rayon), 10, "ms");
        std::fs::remove_file(COLD_PATH)?;
    } else {
        println!("     skipping cold benchmarks: cannot write an uncached file on this target");
    }
    Ok(())
}
