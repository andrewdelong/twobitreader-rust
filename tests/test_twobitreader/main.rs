use twobitreader::TwobitReader;

// Standard library
use std::collections::HashMap;
use std::error::Error;
use std::io;
use std::iter::zip;
use std::path::Path;

// Modules
mod util;
use util::*;

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

#[rustfmt::skip]
fn test_tiny<P: AsRef<Path>>(path: P) -> Result<(), Box<dyn Error>> {
    let tbr = TwobitReader::open_masked(path)?;

    // Reproduce the expected contents of tiny.2bit file here, as strings.
    let expect_names = vec!["seq_xy", "seq_abc"];
    let expect_seqs = vec![
        "gTCCTGCTccaGAAGCAATAACTGATAACNNNNnnnngatcaGCAAGACAATTGAAGAAt",
        "NNGCTgtgcanNNNNNGACTCCTAcctcnn",
    ];
    let expect_twobit = zip(expect_names.iter().copied(), expect_seqs.iter().copied()).collect::<HashMap<_, _>>();

    // Test names and iter_names
    assert_eq!(expect_names, tbr.names());
    assert_eq!(expect_names, tbr.iter_names().collect::<Vec<_>>());

    // Test name lookups (&str, String, &String) and sequence len
    for name in expect_names.iter().copied() {
        assert!(tbr.contains_name(name));
    }

    // Test sequence len
    for name in expect_names.iter().copied() {
        assert_eq!(tbr.seq_len(name), expect_twobit[name].len());
    }

    // Test contains_name negatives (&str, String, &String)
    assert!(!tbr.contains_name(""));
    assert!(!tbr.contains_name("no_such_name"));
    assert!(!tbr.contains_name("no_such_name".to_string()));
    assert!(!tbr.contains_name(&"no_such_name".to_string()));

    // Test iteration over &TwobitReader
    assert_eq!(tbr.num_seqs(), expect_twobit.len());

    // Iterator over all possible subsequences
    let windows = expect_twobit.iter().flat_map(|(name, seq)| {
        (0..seq.len()).flat_map(move |window_size| {
            (0..seq.len() - window_size + 1).map(move |start| {
                (*name, start, start + window_size)
            })
        })
    });

    // Each window's expected sequence
    let expect_batch = windows.clone().map(|(name, start, end)| &expect_twobit[name][start..end]).collect::<Vec<_>>();

    // Test get for every sub-window
    for ((name, start, end), expect_seq) in zip(windows.clone(), expect_batch.iter().copied()) {
        // TwobitReader get methods
        assert_eq!(expect_seq, tbr.get(name, start, end));
        assert_eq!(expect_seq, tbr.get_inclusive(name, start + 1, end));
        assert_eq!(expect_seq, { let mut dst = String::new(); tbr.get_into(name, start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::new(); tbr.get_inclusive_into(name, start + 1, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbr.get_into(name, start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbr.get_inclusive_into(name, start + 1, end, &mut dst); dst });
    }

    // Test concat with zip (from iterator, from collection)
    let expect = "aACNNNcaGCATTGA";
    let starts = [0, 10, 20, 30, 40, 50];
    let ends   = [0, 11, 22, 33, 44, 55];
    let inclusive_starts = starts.iter().map(|start| start + 1);
    assert_eq!(expect, tbr.concat("seq_xy", zip(starts, ends)));
    assert_eq!(expect, tbr.concat("seq_xy", zip(starts, ends).collect::<Vec<_>>()));
    assert_eq!(expect, tbr.concat_inclusive("seq_xy", zip(inclusive_starts.clone(), ends)));
    assert_eq!(expect, tbr.concat_inclusive("seq_xy", zip(inclusive_starts.clone(), ends).collect::<Vec<_>>()));

    // Test concat with array of (start, end) pairs (from array, from iterator, from collection)
    let ranges = [(0, 0), (10, 11), (20, 22), (30, 33), (40, 44), (50, 55)];
    assert_eq!(expect, tbr.concat("seq_xy", ranges));
    assert_eq!(expect, tbr.concat("seq_xy", ranges.into_iter()));
    assert_eq!(expect, tbr.concat("seq_xy", Vec::from(ranges)));

    // Test concat with total_len=0, and different name types (&str, String, &String)
    assert_eq!(tbr.concat("seq_xy", zip(starts, starts)), "");
    assert_eq!(tbr.concat("seq_xy".to_string(), zip(starts, starts)), "");
    assert_eq!(tbr.concat(&"seq_xy".to_string(), zip(starts, starts)), "");

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
#[should_panic(expected = "invalid start")]
fn test_get_inclusive_invalid_start_panic() {
    open_tiny().get_inclusive("seq_xy", 0, 10);
}

#[test]
#[should_panic(expected = "invalid range")]
fn test_get_invalid_range_panic() {
    open_tiny().get("seq_xy", 11, 10);
}

#[test]
#[should_panic(expected = "invalid end")]
fn test_get_invalid_end_panic() {
    open_tiny().get("seq_xy", 0, 61);
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_seq_len_invalid_name_panic() {
    open_tiny().seq_len("seq_XY");
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_get_invalid_name_panic() {
    open_tiny().get("seq_XY", 0, 10);
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_concat_invalid_name_panic() {
    open_tiny().concat("seq_XY", [(0, 1), (2, 3)]);
}

#[test]
#[should_panic(expected = "invalid end")]
fn test_concat_invalid_end_panic() {
    open_tiny().concat("seq_xy", [(0, 1), (2, 61)]);
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_concat_inclusive_invalid_name_panic() {
    open_tiny().concat_inclusive("seq_XY", [(1, 1), (3, 3)]);
}

#[test]
#[should_panic(expected = "invalid end")]
fn test_concat_inclusive_invalid_end_panic() {
    open_tiny().concat("seq_xy", [(1, 1), (3, 61)]);
}

#[test]
#[ignore] // hg38.2bit is not committed to git; run only if requested by user.
fn test_hg38() -> Result<(), Box<dyn Error>> {
    // Check that all sample intervals match reference sequence
    let r = TwobitReader::open_masked(download_hg38()?.as_path())?;
    let fasta = read_test_fasta("tests/assets/hg38-expected-seqs.fasta")?;
    for (chrom, start, end, seq) in fasta {
        assert_eq!(r.get(&chrom, start, end), seq);
    }
    Ok(())
}
