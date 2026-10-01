"""Tests for the twobitreader_rs Python bindings.

These mirror the Rust integration tests in tests/test_twobitreader/main.rs, so that the
bindings are checked against the same expectations as the crate they wrap. Where the Rust
tests assert a panic, the equivalent here asserts the exception PyO3 maps that panic to.
"""

from contextlib import contextmanager
from pathlib import Path

import pytest  # type: ignore

from twobitreader_rs import TwobitReader, reverse_complement


# PyO3 maps a Rust panic to pyo3_runtime.PanicException, which derives from BaseException
# rather than Exception. That module is synthetic and cannot be imported, so the exception
# is matched by name instead of by class.
@contextmanager
def raises_panic(match):
    with pytest.raises(BaseException, match=match) as excinfo:
        yield
    assert type(excinfo.value).__name__ == "PanicException", excinfo.value


# Assets are located relative to this file, so the tests do not depend on the working directory.
ASSETS = Path(__file__).resolve().parents[1] / "assets"

# Expected contents of the tiny .2bit files, in the order the records appear in the file.
EXPECT_NAMES = ["seq_xy", "seq_abc"]
EXPECT_SEQS_MASKED = [
    "gTCCTGCTccaGAAGCAATAACTGATAACNNNNnnnngatcaGCAAGACAATTGAAGAAt",
    "NNGCTgtgcanNNNNNGACTCCTAcctcnn",
]
EXPECT_SEQS_UNMASKED = [
    "GTCCTGCTCCAGAAGCAATAACTGATAACNNNNNNNNGATCAGCAAGACAATTGAAGAAT",
    "NNGCTGTGCANNNNNNGACTCCTACCTCNN",
]

# Every (name, start, end) sub-window of both sequences, matching the Rust test's generator.
# This covers all four start alignments and every quad/tail split in the decoder.
def all_windows(expect):
    for name, seq in expect.items():
        for size in range(len(seq) + 1):
            for start in range(len(seq) - size + 1):
                yield (name, start, start + size)


def open_tiny():
    return TwobitReader.open(str(ASSETS / "tiny-lilend.2bit"))


def test_malformed_header_error():
    # File has bad Twobit signature.
    with pytest.raises(OSError, match="signature"):
        TwobitReader.open(str(ASSETS / "malformed-bad-signature.2bit"))

    # File has a duplicate sequence name.
    with pytest.raises(OSError, match="duplicate sequence name"):
        TwobitReader.open(str(ASSETS / "malformed-duplicate-name.2bit"))


def test_missing_file_error():
    with pytest.raises(FileNotFoundError):
        TwobitReader.open(str(ASSETS / "no-such-file.2bit"))


def test_truncated_blocks_error():
    # File block indices were truncated by end of file.
    tbr = TwobitReader.open_masked(str(ASSETS / "malformed-truncated-blocks.2bit"))
    with raises_panic(match="num_blocks too large"):
        for name in tbr.names():
            tbr.get(name, 0, tbr.seq_len(name))


def test_truncated_dna_error():
    # File DNA sequence was truncated by end of file.
    tbr = TwobitReader.open_masked(str(ASSETS / "malformed-truncated-dna.2bit"))
    with raises_panic(match="Failed to read DNA"):
        for name in tbr.names():
            tbr.get(name, 0, tbr.seq_len(name))


# Each tiny file is checked in both endiannesses and with the mask applied or not,
# mirroring test_tiny_lilend / test_tiny_lilend_masked / test_tiny_bigend / test_tiny_bigend_masked.
@pytest.mark.parametrize("endian", ["lilend", "bigend"])
@pytest.mark.parametrize("masked", [False, True])
def test_tiny(endian, masked):
    path = str(ASSETS / f"tiny-{endian}.2bit")
    tbr = TwobitReader.open_masked(path) if masked else TwobitReader.open(path)
    expect_seqs = EXPECT_SEQS_MASKED if masked else EXPECT_SEQS_UNMASKED
    expect = dict(zip(EXPECT_NAMES, expect_seqs))

    # Test names.
    assert tbr.names() == EXPECT_NAMES

    # Test name lookups and sequence length.
    for name in EXPECT_NAMES:
        assert tbr.contains_name(name)
        assert tbr.seq_len(name) == len(expect[name])

    # Test contains_name negatives.
    assert not tbr.contains_name("")
    assert not tbr.contains_name("no_such_name")

    windows = list(all_windows(expect))
    expect_batch = [expect[name][start:end] for name, start, end in windows]

    # Test get and get_inclusive for every sub-window.
    for (name, start, end), expect_seq in zip(windows, expect_batch):
        assert tbr.get(name, start, end) == expect_seq
        assert tbr.get_inclusive(name, start + 1, end) == expect_seq

    # Test prefetch with all valid ranges, as a list and as a generator. The windows include
    # empty ranges, which in 1-based form are (start, start - 1).
    assert tbr.prefetch(windows) is None
    assert tbr.prefetch(w for w in windows) is None
    windows_inclusive = [(n, s + 1, e) for n, s, e in windows]
    assert tbr.prefetch_inclusive(windows_inclusive) is None
    assert tbr.prefetch_inclusive(w for w in windows_inclusive) is None

    # Test concat over ranges of differing length, including an empty one.
    expect_concat = "aACNNNcaGCATTGA" if masked else "AACNNNCAGCATTGA"
    ranges = [(0, 0), (10, 11), (20, 22), (30, 33), (40, 44), (50, 55)]
    ranges_inclusive = [(start + 1, end) for start, end in ranges]

    # Test concat with a list, a tuple, and a generator (the latter is iterable only once,
    # which selects the grow-the-buffer branch of the binding rather than the pre-sized one).
    assert tbr.concat("seq_xy", ranges) == expect_concat
    assert tbr.concat("seq_xy", tuple(ranges)) == expect_concat
    assert tbr.concat("seq_xy", (r for r in ranges)) == expect_concat
    assert tbr.concat_inclusive("seq_xy", ranges_inclusive) == expect_concat
    assert tbr.concat_inclusive("seq_xy", (r for r in ranges_inclusive)) == expect_concat

    # Test concat with total length zero.
    empty = [(start, start) for start, _ in ranges]
    assert tbr.concat("seq_xy", empty) == ""
    assert tbr.concat("seq_xy", (r for r in empty)) == ""
    assert tbr.concat("seq_xy", []) == ""

    # Test get_batch and get_batch_inclusive, as a list (sized) and a generator (unsized).
    assert tbr.get_batch(windows) == expect_batch
    assert tbr.get_batch(w for w in windows) == expect_batch
    assert tbr.get_batch_inclusive(windows_inclusive) == expect_batch
    assert tbr.get_batch_inclusive(w for w in windows_inclusive) == expect_batch


def test_get_inclusive_invalid_start_error():
    with pytest.raises(ValueError, match="start"):
        open_tiny().get_inclusive("seq_xy", 0, 10)


def test_get_invalid_range_error():
    with pytest.raises(ValueError, match="invalid range"):
        open_tiny().get("seq_xy", 11, 10)


def test_get_invalid_end_error():
    with pytest.raises(ValueError, match="invalid end"):
        open_tiny().get("seq_xy", 0, 61)


# An unknown sequence name is reported as KeyError by every method that takes one, the way a
# dict reports a missing key. The Rust crate panics in the same situation; the bindings resolve
# the name through a fallible lookup so that the panic never reaches Python.
@pytest.mark.parametrize(
    "call",
    [
        lambda tbr: tbr.seq_len("seq_XY"),
        lambda tbr: tbr.get("seq_XY", 0, 10),
        lambda tbr: tbr.get_inclusive("seq_XY", 1, 10),
        lambda tbr: tbr.get_batch([("seq_XY", 0, 10)]),
        lambda tbr: tbr.get_batch_inclusive([("seq_XY", 1, 10)]),
        lambda tbr: tbr.concat("seq_XY", [(0, 1), (2, 3)]),
        lambda tbr: tbr.concat_inclusive("seq_XY", [(1, 1), (3, 3)]),
        lambda tbr: tbr.prefetch([("seq_XY", 0, 10)]),
        lambda tbr: tbr.prefetch_inclusive([("seq_XY", 1, 10)]),
    ],
)
def test_invalid_name_raises_key_error(call):
    with pytest.raises(KeyError, match="seq_XY"):
        call(open_tiny())


def test_prefetch_invalid_range_error():
    with pytest.raises(ValueError, match="invalid range"):
        open_tiny().prefetch([("seq_xy", 11, 10)])


def test_prefetch_inclusive_invalid_range_error():
    with pytest.raises(ValueError, match="invalid range"):
        open_tiny().prefetch_inclusive([("seq_xy", 12, 10)])


def test_prefetch_inclusive_invalid_start_error():
    with pytest.raises(ValueError, match="invalid start"):
        open_tiny().prefetch_inclusive([("seq_xy", 0, 10)])


# A negative index cannot convert to the unsigned type the bindings take, so PyO3 rejects it
# before any of our own checks run. OverflowError is what Python raises elsewhere for the same
# reason, so it is the contract rather than a quirk; note that it is not a ValueError.
@pytest.mark.parametrize(
    "call",
    [
        lambda tbr: tbr.get("seq_xy", -1, 10),
        lambda tbr: tbr.get("seq_xy", 0, -1),
        lambda tbr: tbr.get_inclusive("seq_xy", -1, 10),
        lambda tbr: tbr.get_batch([("seq_xy", -1, 10)]),
        lambda tbr: tbr.get_batch_inclusive([("seq_xy", -1, 10)]),
        lambda tbr: tbr.concat("seq_xy", [(-1, 10)]),
        lambda tbr: tbr.concat_inclusive("seq_xy", [(-1, 10)]),
        lambda tbr: tbr.prefetch([("seq_xy", -1, 10)]),
        lambda tbr: tbr.prefetch_inclusive([("seq_xy", -1, 10)]),
    ],
)
def test_negative_index_raises_overflow_error(call):
    with pytest.raises(OverflowError, match="negative"):
        call(open_tiny())


# A non-integer index is likewise rejected during argument conversion.
@pytest.mark.parametrize(
    "call",
    [
        lambda tbr: tbr.get("seq_xy", 1.5, 10),
        lambda tbr: tbr.get("seq_xy", "0", 10),
        lambda tbr: tbr.get_batch([("seq_xy", 1.5, 10)]),
        lambda tbr: tbr.concat("seq_xy", [(1.5, 10)]),
        lambda tbr: tbr.get_batch(["not a tuple"]),
    ],
)
def test_non_integer_index_raises_type_error(call):
    with pytest.raises(TypeError):
        call(open_tiny())


# bool is a subclass of int in Python, so True is accepted as the index 1. Not a defect,
# but worth pinning so the behaviour is not changed by accident.
def test_bool_index_is_accepted_as_int():
    tbr = open_tiny()
    assert tbr.get("seq_xy", True, 10) == tbr.get("seq_xy", 1, 10)


def test_concat_invalid_end_error():
    with pytest.raises(ValueError, match="invalid end"):
        open_tiny().concat("seq_xy", [(0, 1), (2, 61)])


def test_concat_inclusive_invalid_start_error():
    with pytest.raises(ValueError, match="invalid start"):
        open_tiny().concat_inclusive("seq_xy", [(0, 5)])


# Malformed arguments should raise rather than panic or silently misbehave.
def test_batch_item_type_errors():
    tbr = open_tiny()
    with pytest.raises(TypeError):
        tbr.prefetch([(1, 2, 3)])
    with pytest.raises(ValueError):
        tbr.get_batch([("seq_xy", 0, 1, 9)])


def test_reverse_complement():
    def check_pair(x, y):
        assert x == reverse_complement(y)
        assert y == reverse_complement(x)

    check_pair("", "")
    check_pair("A", "T")
    check_pair("C", "G")
    check_pair("N", "N")
    check_pair("a", "t")
    check_pair("c", "g")
    check_pair("n", "n")
    check_pair("AC", "GT")
    check_pair("ac", "gt")
    check_pair("ACG", "CGT")
    check_pair("acg", "cgt")
    check_pair("ACGT", "ACGT")
    check_pair("acgt", "acgt")
    check_pair("ACGTN", "NACGT")
    check_pair("acgtn", "nacgt")
    check_pair(
        "ACGTNaaccggttnnAAACCCGGGTTTNNNaaaaccccggggttttnnnn",
        "nnnnaaaaccccggggttttNNNAAACCCGGGTTTnnaaccggttNACGT",
    )


def test_reverse_complement_non_string():
    with pytest.raises(TypeError):
        reverse_complement(1)
