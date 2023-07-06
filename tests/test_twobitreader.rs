use twobitreader::TwobitReader;

// Standard library
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};
use std::iter::{zip, once, repeat};
use std::error::Error;
use std::time::Duration;
use std::collections::HashMap;

// Dependencies
use reqwest::blocking::Client;

fn open_tiny() -> TwobitReader {
    let r = TwobitReader::open("tests/assets/tiny-lilend.2bit");
    assert!(r.is_ok());
    r.unwrap()    
}

#[test]
fn test_malformed_file_error() {
    // File has bad Twobit signature
    let r = TwobitReader::open("tests/assets/malformed-bad-signature.2bit");
    assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::InvalidData);

    // File has a duplicate sequence name
    let r = TwobitReader::open("tests/assets/malformed-duplicate-name.2bit");
    assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::InvalidData);

    // File block indices were truncated by end of file
    let r = TwobitReader::open("tests/assets/malformed-truncated-blocks.2bit");
    assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::UnexpectedEof);

    // File DNA sequence was truncated by end of file
    let r = TwobitReader::open("tests/assets/malformed-truncated-dna.2bit");
    assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::UnexpectedEof);
}

fn test_tiny<P: AsRef<Path>>(path: P) -> Result<(), Box<dyn Error>> {
    let tbr = TwobitReader::open_masked(path)?;

    // Reproduce the expected contents of tiny.2bit file here, as strings.
    let expect_names = vec!["seq_xy", "seq_abc"];
    let expect_seqs = vec!["gTCCTGCTccaGAAGCAATAACTGATAACNNNNnnnngatcaGCAAGACAATTGAAGAAt",
                           "NNGCTgtgcanNNNNNGACTCCTAcctcnn"];
    let expect_twobit = zip(expect_names.iter().copied(),
                            expect_seqs.iter().copied()).collect::<HashMap<_,_>>();

    // Test names and iter_names
    assert_eq!(expect_names, tbr.names());
    assert_eq!(expect_names, tbr.iter_names().collect::<Vec<_>>());

    // Test name lookups and sequence len
    for name in expect_names.iter().copied() {
        assert!(tbr.contains_name(name));
        assert_eq!(tbr[name].name(), name);               // name: &str 
        assert_eq!(tbr[name.to_string()].name(), name);   // name: String
        assert_eq!(tbr[&name.to_string()].name(), name);  // name: &String
    }

    // Test sequence len
    for name in expect_names.iter().copied() {
        assert_eq!(tbr[name].len(), expect_twobit[name].len());
    }

    // Test contains_name negatives
    assert!(!tbr.contains_name(""));
    assert!(!tbr.contains_name("no_such_name"));               // name: &str
    assert!(!tbr.contains_name("no_such_name".to_string()));   // name: String
    assert!(!tbr.contains_name(&"no_such_name".to_string()));  // name: &String

    // Test iteration over &TwobitReader
    assert_eq!(tbr.len(), expect_twobit.len());
    for tbs in &tbr {
        assert!(expect_twobit.contains_key(tbs.name()));        
    }

    // Iterator over all possible subsequences
    let windows = expect_twobit.iter().flat_map(|(name, seq)| {
        (0..seq.len()).flat_map(move |window_size| {
            (0..seq.len()-window_size+1).map(move |start| {
                (*name, start, start+window_size)
            })
        })
    });

    // Each window's expected sequence
    let expect_batch = windows.clone().map(|(name, start, end)| &expect_twobit[name][start..end]).collect::<Vec<_>>();

    // Test get for every sub-window
    for ((name, start, end), expect_seq) in zip(windows.clone(), expect_batch.iter().copied()) {
        // TwobitReader get methods
        assert_eq!(expect_seq, tbr.get(name, start, end));
        assert_eq!(expect_seq, tbr.get_inclusive(name, start+1, end));
        assert_eq!(expect_seq, { let mut dst = String::new(); tbr.get_into(name, start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::new(); tbr.get_inclusive_into(name, start+1, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbr.get_into(name, start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbr.get_inclusive_into(name, start+1, end, &mut dst); dst });

        // TwobitSequence get methods
        let tbs = &tbr[name];
        assert_eq!(expect_seq, tbs.get(start, end));
        assert_eq!(expect_seq, tbs.get_inclusive(start+1, end));
        assert_eq!(expect_seq, { let mut dst = String::new(); tbs.get_into(start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::new(); tbs.get_inclusive_into(start+1, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbs.get_into(start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbs.get_inclusive_into(start+1, end, &mut dst); dst });
    }

    // Test get_batch on all sub-windows at once
    let to_inclusive = |(name, start, end)| (name, start+1, end);
    assert_eq!(expect_batch, tbr.get_batch(windows.clone()));                      // as iterator
    assert_eq!(expect_batch, tbr.get_batch(windows.clone().collect::<Vec<_>>()));  // as collection
    assert_eq!(expect_batch, tbr.get_inclusive_batch(windows.clone().map(to_inclusive)));
    assert_eq!(expect_batch, tbr.get_inclusive_batch(windows.clone().map(to_inclusive).collect::<Vec<_>>()));

    // Test get_batch on empty
    assert_eq!(tbr.get_batch(Vec::<(String, usize, usize)>::new().into_iter()), Vec::<String>::new());
    assert_eq!(tbr.get_inclusive_batch(Vec::<(String, usize, usize)>::new().into_iter()), Vec::<String>::new());

    // Test concat with zip
    let expect = "aACNNNcaGCATTGA";
    let starts = [0, 10, 20, 30, 40, 50];
    let ends   = [0, 11, 22, 33, 44, 55];
    let inclusive_starts = starts.iter().map(|start| start+1);
    assert_eq!(expect, tbr.concat("seq_xy", zip(starts, ends)));                      // ranges: Iterator
    assert_eq!(expect, tbr.concat("seq_xy", zip(starts, ends).collect::<Vec<_>>()));  // ranges: Vec
    assert_eq!(expect, tbr["seq_xy"].concat(zip(starts, ends)));
    assert_eq!(expect, tbr["seq_xy"].concat(zip(starts, ends).collect::<Vec<_>>()));
    assert_eq!(expect, tbr.concat_inclusive("seq_xy", zip(inclusive_starts.clone(), ends)));
    assert_eq!(expect, tbr.concat_inclusive("seq_xy", zip(inclusive_starts.clone(), ends).collect::<Vec<_>>()));
    assert_eq!(expect, tbr["seq_xy"].concat_inclusive(zip(inclusive_starts.clone(), ends)));
    assert_eq!(expect, tbr["seq_xy"].concat_inclusive(zip(inclusive_starts.clone(), ends).collect::<Vec<_>>()));

    // Test concat with array of (start, end) pairs
    let ranges = [(0, 0), (10, 11), (20, 22), (30, 33), (40, 44), (50, 55)];
    assert_eq!(expect, tbr.concat("seq_xy", ranges));              // ranges: Array
    assert_eq!(expect, tbr.concat("seq_xy", ranges.into_iter()));  // ranges: Iterator
    assert_eq!(expect, tbr.concat("seq_xy", Vec::from(ranges)));   // ranges: Vec
    assert_eq!(expect, tbr["seq_xy"].concat(ranges));
    assert_eq!(expect, tbr["seq_xy"].concat(ranges.into_iter()));
    assert_eq!(expect, tbr["seq_xy"].concat(Vec::from(ranges)));

    // Test concat with total_len=0, and different name types
    assert_eq!(tbr["seq_xy"].concat(zip(starts, starts)), "");
    assert_eq!(tbr.concat("seq_xy", zip(starts, starts)), "");                // name: &str
    assert_eq!(tbr.concat("seq_xy".to_string(), zip(starts, starts)), "");    // name: String
    assert_eq!(tbr.concat(&"seq_xy".to_string(), zip(starts, starts)), "");   // name: &String

    Ok(())
}

// Run test_tiny on the little-endian version of test.2bit
#[test]
fn test_tiny_lilend() -> Result<(), Box<dyn Error>> {
    test_tiny("tests/assets/tiny-lilend.2bit")
}

// Run test_tiny on the big-endian version of test.2bit
#[test]
fn test_tiny_bigend() -> Result<(), Box<dyn Error>> {
    test_tiny("tests/assets/tiny-bigend.2bit")
}

#[test]
#[should_panic(expected="invalid start")]
fn test_get_inclusive_invalid_start_panic() {
    open_tiny().get_inclusive("seq_xy", 0, 10);
}

#[test]
#[should_panic(expected="invalid range")]
fn test_get_invalid_range_panic() {
    open_tiny().get("seq_xy", 11, 10);
}

#[test]
#[should_panic(expected="invalid end")]
fn test_get_invalid_end_panic() {
    open_tiny().get("seq_xy", 0, 61);
}

#[test]
#[should_panic(expected="sequence name not found")]
fn test_index_invalid_name_panic() {
    let _tbs = &open_tiny()["seq_XY"];
}

#[test]
#[should_panic(expected="sequence name not found")]
fn test_get_invalid_name_panic() {
    open_tiny().get("seq_XY", 0, 10);
}

#[test]
#[should_panic(expected="sequence name not found")]
fn test_get_batch_invalid_name_panic() {
    // Create large batch to spawn additional threads, to test that panic is propagated.
    let args = repeat(("seq_xy", 0, 10)).take(1000).chain(once(("seq_XY", 0, 10)));
    open_tiny().get_batch(args);
}

#[test]
#[should_panic(expected="invalid end")]
fn test_get_batch_invalid_end_panic() {
    // Create large batch to spawn additional threads, to test that panic is propagated.
    let args = repeat(("seq_xy", 0, 10)).take(1000).chain(once(("seq_xy", 0, 61)));
    open_tiny().get_batch(args);
}

#[test]
#[should_panic(expected="invalid start")]
fn test_get_inclusive_batch_invalid_start_panic() {
    // Create large batch to spawn additional threads, to test that panic is propagated.
    let args = repeat(("seq_xy", 1, 10)).take(1000).chain(once(("seq_xy", 0, 10)));
    open_tiny().get_inclusive_batch(args);
}

#[test]
#[should_panic(expected="sequence name not found")]
fn test_concat_invalid_name_panic() {
    open_tiny().concat("seq_XY", [(0, 1), (2, 3)]);
}

#[test]
#[should_panic(expected="invalid end")]
fn test_concat_invalid_end_panic() {
    open_tiny().concat("seq_xy", [(0, 1), (2, 61)]);
}

// Reads a FASTA file where description lines all adhere to ">chrom:start-end" format.
// The returned Vec has structure:
//   Vec<(chrom, start, end, seq)>
fn read_test_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<(String, usize, usize, String)>, Box<dyn Error>> {
    let lines = io::BufReader::new(File::open(path)?).lines();
    let mut r = Vec::new();
    for line in lines {  // TODO: Lines<B> creates new string each iteration; instead call read_line directly or read_until('\n') which is even faster for ASCII
        let line = line?;
        if line.starts_with('>') {
            // Parse ">chrom:start-end"
            let i = line.find(':').ok_or("expected description line format chrom:start-end")?;
            let j = line.find('-').ok_or("expected description line format chrom:start-end")?;
            let chrom = line[1..i].to_string();
            let start = line[i+1..j].parse::<usize>()?;
            let end = line[j+1..].parse::<usize>()?;
            r.push((chrom, start, end, String::new()));
        } else if !line.is_empty() && !line.starts_with(';') {
            let seq = &mut r.last_mut().ok_or("expected description line before sequence")?.3;  // .3 = String of dna
            seq.push_str(&line);
        }
    }
    Ok(r)
}

// Returns a local path to the hg38 2bit file used for tests on human reference genome.
// If the local file does not yet exist, the function will attempt to download it from UCSC.
// WARNING: Should not be called by multiple tests. Tests are run in parallel, and this code
//          has an obvious race condition on the output file if used in multiple tests.
fn download_hg38() -> Result<PathBuf, Box<dyn Error>> {
    let path = Path::new("tests/assets/hg38.p13.2bit");
    if !path.exists() {
        // Download binary contents of 2bit and write to a local file.
        let url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.2bit";
        let client = Client::builder()
                    .timeout(Duration::new(240, 0))   // Give 4 minutes to download ~800MB
                    .build()?;
        let data = client.get(url).send()?.error_for_status()?.bytes()?;
        let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
        file.write_all(&data)?;
    }
    Ok(path.to_path_buf())
}

#[test]
#[ignore]  // hg38.2bit is not committed to git; run only if requested by user.
fn test_hg38() -> Result<(), Box<dyn Error>> {
    // Check that all sample intervals match reference sequence
    let r = TwobitReader::open_masked(download_hg38()?.as_path())?;
    let fasta = read_test_fasta("tests/assets/hg38-expected-seqs.fasta")?;
    for (chrom, start, end, seq) in fasta {
        assert_eq!(r.get(&chrom, start, end), seq);
    }
    Ok(())
}
