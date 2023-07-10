// Standard library
use std::error::Error;
use std::fs::File;
use std::hint::black_box;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

// Dependencies
use flate2::read::GzDecoder;
use itertools::Itertools;
use memmap2::MmapOptions;

// When a bench callback is invoked several times, only this many will
// be printed to console describing the result (for legibility).
const NUM_BENCH_TIMES_PRINTED: usize = 5;

// Returns mean of a vector of elapsed time values.
fn mean_time(x: &Vec<u128>) -> u128 {
    assert!(!x.is_empty());
    let n = x.iter().copied().sum::<u128>();
    n / (x.len() as u128)
}

// Returns median of a vector of elapsed time values.
fn median_time(x: &Vec<u128>) -> u128 {
    assert!(!x.is_empty());
    let mut y = x.to_owned();
    y.sort();
    y[y.len() / 2].clone()
}

// Invokes f() num_iter times and prints the mean/median/min/max of the reported durations.
// The time durations are reported as integers in the specified unit ("s", "ms", "us", "ns").
pub fn bench(name: &str, f: fn() -> Result<Duration, Box<dyn Error>>, num_iter: usize, unit: &str) {
    // Convert duration to particular time unit as u128.
    let as_unit = match unit {
        "ns" => |d: Duration| d.as_nanos(),
        "us" => |d: Duration| d.as_micros(),
        "ms" => |d: Duration| d.as_millis(),
        "s" => |d: Duration| d.as_secs().into(),
        _ => panic!("invalid unit '{unit}'"),
    };

    // Invoke the callback num_iter times, recording the reported duration each time.
    let mut times = Vec::with_capacity(num_iter);
    for _ in 0..num_iter {
        let duration = f().unwrap();
        times.push(as_unit(duration));
    }

    // Compute summary stats from the times.
    let avg_time = mean_time(&times);
    let med_time = median_time(&times);
    let min_time = times.iter().min().unwrap();
    let max_time = times.iter().max().unwrap();

    // Format the output string, in a way that also reports the first few measurements.
    let mut avg_str = format!("{avg_time}");
    let mut times_str = String::default();
    let mut stats_str = String::default();
    if times.len() > 1 {
        avg_str.insert_str(0, "avg=");
        stats_str = format!("(med={med_time} min={min_time} max={max_time})");
        times_str = format!("{:?}", &times[..NUM_BENCH_TIMES_PRINTED.min(times.len())]);
        if times.len() > NUM_BENCH_TIMES_PRINTED {
            times_str = times_str.replace("]", ", ...]");
        }
    }

    // Print bench result indented, and format the test name in green
    println!("     \x1b[32m{name}\x1b[0m: {avg_str} {unit} {stats_str} {times_str}");
}

// Reads the first 4 columns of a bed.gz file. The returned Vec has structure:
//   Vec<(chrom, start, end, name)>
pub fn read_bed<P: AsRef<Path>>(path: P) -> Result<Vec<(String, usize, usize, String)>, Box<dyn Error>> {
    let reader = io::BufReader::new(GzDecoder::new(File::open(path)?));
    let mut result = Vec::new();
    // Using read_until('\n') faster than lines(), but negligible benefit for .gz files
    for line in reader.lines() {
        let line = line?;
        let mut cols = line.split('\t');
        let chrom = cols.next().ok_or("expected >=4 tab-separated columns")?.to_string();
        let start = cols.next().ok_or("expected >=4 tab-separated columns")?.parse::<usize>()?;
        let end = cols.next().ok_or("expected >=4 tab-separated columns")?.parse::<usize>()?;
        let name = cols.next().ok_or("expected >=4 tab-separated columns")?.to_string();
        result.push((chrom, start, end, name));
    }
    Ok(result)
}

// Reads gencode-exons.bed.gz intervals. The returned Vec has structure:
//   Vec<(chrom, start, end)>
pub fn read_exons() -> Result<Vec<(String, usize, usize)>, Box<dyn Error>> {
    let bed = read_bed("benches/assets/gencode-exons.bed.gz")?;
    let exons = bed.into_iter().map(|(chrom, start, end, _)| (chrom, start, end)).collect();
    Ok(exons)
}

// Reads gencode-transcripts.bed.gz and groups the intervals by transcript ID (BED name field).
// The returned Vec has structure:
//   Vec<(transcript_id, chrom, Vec<(start, end)>)>
pub fn read_transcripts() -> Result<Vec<(String, String, Vec<(usize, usize)>)>, Box<dyn Error>> {
    let bed = read_bed("benches/assets/gencode-transcripts.bed.gz")?;
    let transcripts = bed.into_iter().map(|(chrom, start, end, id)| ((id, chrom), (start, end))).into_group_map();
    let transcripts = transcripts.into_iter().map(|((id, chrom), exons)| (id, chrom, exons)).collect();
    Ok(transcripts)
}

// Maps the given file and accesses all pages, so that it resides in the page cache.
// If there is sufficient unused physical memory, the operating system will keep
// the data in memory (hot) for subsequent reuse.
pub fn load_into_page_cache<P: AsRef<Path>>(path: P) -> io::Result<()> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    // Access all pages by computing an arbitrary value.
    // (Avoid mmap madvise because it's not portable.)
    black_box(mmap.iter().max());
    Ok(())
}

// Returns a local path to the hg38 2bit file used for tests on human reference genome.
// If the local file does not yet exist, the function will attempt to download it from UCSC.
pub fn download_hg38() -> Result<PathBuf, Box<dyn Error>> {
    let path = Path::new("tests/assets/hg38.p13.2bit");
    if !path.exists() {
        let url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.2bit";
        let status_code = Command::new("curl")
            .args(["-o", path.to_str().unwrap(), url])
            .status()?
            .code()
            .expect("Curl failed to terminate successfully");
        if status_code != 0 {
            panic!("Curl failed with exit code {status_code}");
        }
    }
    Ok(path.to_path_buf())
}
