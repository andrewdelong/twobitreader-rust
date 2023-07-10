// Standard library
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::process::Command;

// Reads a FASTA file where description lines all adhere to ">chrom:start-end" format.
// The returned Vec has structure:
//   Vec<(chrom, start, end, seq)>
pub fn read_test_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<(String, usize, usize, String)>, Box<dyn Error>> {
    let lines = io::BufReader::new(File::open(path)?).lines();
    let mut r = Vec::new();
    // TODO: Lines<B> creates new string each iteration; instead call read_line
    // directly or read_until('\n') which is even faster for ASCII
    for line in lines {
        let line = line?;
        if line.starts_with('>') {
            // Parse ">chrom:start-end"
            let i = line.find(':').ok_or("expected description line chrom:start-end")?;
            let j = line.find('-').ok_or("expected description line chrom:start-end")?;
            let chrom = line[1..i].to_string();
            let start = line[i + 1..j].parse::<usize>()?;
            let end = line[j + 1..].parse::<usize>()?;
            r.push((chrom, start, end, String::new()));
        } else if !line.is_empty() && !line.starts_with(';') {
            let seq = &mut r.last_mut().ok_or("expected description line before sequence")?.3; // .3 = String of dna
            seq.push_str(&line);
        }
    }
    Ok(r)
}

// Returns a local path to the hg38 2bit file used for tests on human reference genome.
// If the local file does not yet exist, the function will attempt to download it from UCSC.
// WARNING: Should not be called by multiple tests. Tests are run in parallel, and this code
//          has an obvious race condition on the output file if used in multiple tests.
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
