//! PyO3 bindings for the twobitreader crate, built as the `twobitreader_rs` extension module.

// Standard library
use std::cell::Cell;
use std::io;
use std::iter::repeat_n;
use std::mem::MaybeUninit;
use std::slice;

// Crate modules
use crate::TwobitSequenceData;

// Dependencies
use memmap2::Mmap;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::ffi::{Py_DECREF, PyObject, PyUnicode_1BYTE_DATA, PyUnicode_CheckExact, PyUnicode_New};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::{PyList, PyNone, PyString};
use scopeguard::ScopeGuard;

/// A reader for a 2bit file.
///
/// The constructors are `TwobitReader.open` and `TwobitReader.open_masked`.
/// Instances are immutable and safe to share between threads.
///
///     >>> tbr = TwobitReader.open("hg38.2bit")
///     >>> tbr.get("chr2", 10000, 10010)
///     'CGTATCCCAC'
#[pyclass(frozen)]
#[derive(Debug)]
pub struct TwobitReader {
    inner: crate::TwobitReader,
}

#[pymethods]
impl TwobitReader {
    /// Open a 2bit file for reading.
    ///
    /// Args:
    ///     path: Path to the 2bit file.
    ///
    /// Returns:
    ///     TwobitReader: A reader for the file.
    ///
    /// Raises:
    ///     FileNotFoundError: The path does not exist.
    ///     OSError: The file is not a 2bit file, or its header is malformed.
    ///
    ///     >>> tbr = TwobitReader.open("hg38.2bit")  # Human genome, build 38
    ///     >>> tbr.get("chr2", 10000, 10010)
    ///     'CGTATCCCAC'
    #[staticmethod]
    pub fn open(path: &str) -> io::Result<Self> {
        crate::TwobitReader::open(path).map(|inner| TwobitReader { inner })
    }

    /// Open a 2bit file for reading, with the lowercase mask applied.
    ///
    /// The lowercase mask feature of 2bit files is mainly relevant for
    /// sequence search, for example with the
    /// [BLAT suite](https://genome.ucsc.edu/goldenpath/help/blatSpec.html) of tools.
    /// Specifically, lowercase letters typically indicate a region that should be ignored
    /// (not searched) during a sequence search.
    ///
    /// Args:
    ///     path: Path to the 2bit file.
    ///
    /// Returns:
    ///     TwobitReader: A reader that returns mixed-case sequences.
    ///
    /// Raises:
    ///     FileNotFoundError: The path does not exist.
    ///     OSError: The file is not a 2bit file, or its header is malformed.
    ///
    ///     >>> tbr = TwobitReader.open_masked("hg38.2bit")
    ///     >>> tbr.get("chr2", 10000, 10010)
    ///     'CGTATcccac'
    #[staticmethod]
    pub fn open_masked(path: &str) -> io::Result<Self> {
        crate::TwobitReader::open_masked(path).map(|inner| TwobitReader { inner })
    }

    /// Return the sequence record names, in the order they appear in the file.
    ///
    /// Returns:
    ///     list[str]: Names such as `["chr1", "chr2", ...]`.
    pub fn names(&self) -> Vec<String> {
        self.inner.iter_names().map(|s| s.to_string()).collect()
    }

    /// Return whether `name` matches a sequence record in the file.
    ///
    /// Args:
    ///     name: Sequence record name to look for.
    ///
    /// Returns:
    ///     bool: True if the record exists.
    pub fn contains_name(&self, name: &str) -> bool {
        self.inner.contains_name(name)
    }

    /// Return the number of nucleotides in the named sequence record.
    ///
    /// Args:
    ///     name: Sequence record name.
    ///
    /// Returns:
    ///     int: Number of nucleotides in the record.
    ///
    /// Raises:
    ///     KeyError: No record has that name.
    ///     BaseException: `pyo3_runtime.PanicException` if the record is malformed.
    pub fn seq_len(&self, name: &str) -> PyResult<usize> {
        Ok(self.seq_data(name)?.dna_len)
    }

    /// Extract the range `start:end` (0-based, exclusive) from the named sequence record.
    ///
    /// Args:
    ///     name: Sequence record name.
    ///     start: 0-based index of the first nucleotide.
    ///     end: 0-based index one past the last nucleotide.
    ///
    /// Returns:
    ///     str: Nucleotides as uppercase `ACGTN`, or mixed case if the reader was opened
    ///         with `open_masked`. Unknown bases are `N`.
    ///
    /// Raises:
    ///     KeyError: No record has that name.
    ///     OverflowError: `start` or `end` is negative.
    ///     TypeError: `start` or `end` is not an integer.
    ///     ValueError: `start > end`, or `end` is past the end of the record.
    ///     BaseException: `pyo3_runtime.PanicException` if the record is malformed.
    ///
    ///     >>> tbr.get("chr2", 10000, 10010)
    ///     'CGTATCCCAC'
    pub fn get<'py>(&self, py: Python<'py>, name: &str, start: usize, end: usize) -> PyResult<Bound<'py, PyString>> {
        // The implementation used for get decodes directly into a PyUnicode buffer and returns, eliding any intermediate
        // allocation or copying. However, it is only ~20% faster than the one-line version that fully leverages PyO3:
        //
        //    pub fn get(&self, name: &str, start: usize, end: usize) -> String {
        //        self.inner.get(name, start, end)
        //    }
        //
        // Since this library is low-level, for a very stable file format, and is unlikely to change, the
        // faster-but-harder-to-maintain version below is still used, and likewise for several other methods.

        // Check range before subtracting. If invalid, raise a Python error here.
        let seq = self.seq_data(name)?;
        check_range(seq, start, end)?;

        // Pre-allocate a PyUnicode instance and decode directly into it.
        // SAFETY: the unicode is initialized with kind ASCII with length matching
        // the number of bytes to be written by decode_from_mmap. Then, fill_blocks
        // can safely use the buffer as initialized (assume_init to drop MaybeUninit).
        // The PyUnicode is then valid ASCII and is safe to cast as PyString.
        let (unicode_guard, buf) = new_unicode_guarded(py, end - start)?;
        decode_and_fill_blocks(&self.inner.mmap, seq, start, buf);
        Ok(as_bound_unicode(py, unicode_guard))
    }

    /// Extract a batch of sequences in one call.
    ///
    /// Equivalent to `[tbr.get(*arg) for arg in batch]`, so mainly provided for concision
    /// and for consistency with `prefetch`.
    ///
    /// Args:
    ///     batch: Iterable of `(name, start, end)` tuples, 0-based with exclusive end.
    ///
    /// Returns:
    ///     list[str]: One sequence per input range, in input order.
    ///
    /// Raises:
    ///     KeyError: An item names a record that does not exist.
    ///     OverflowError: An item has a negative index.
    ///     TypeError: An item is not a 3-tuple of `(str, int, int)`.
    ///     ValueError: An item has an invalid range.
    ///
    ///     >>> tbr.get_batch([("chr1", 10000, 15000), ("chr2", 30000, 35000)])
    ///     ['TAACC...', 'CGTAT...']
    pub fn get_batch<'py>(&self, py: Python<'py>, batch: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyList>> {
        self.get_batch_impl(py, batch, 0)
    }

    /// A version of `get` using 1-based inclusive ranges.
    ///
    /// See https://standage.github.io/on-genomic-interval-notation.html.
    ///
    /// Args:
    ///     name: Sequence record name.
    ///     start: 1-based index of the first nucleotide.
    ///     end: 1-based index of the last nucleotide, included.
    ///
    /// Returns:
    ///     str: The nucleotides in the range.
    ///
    /// Raises:
    ///     KeyError: No record has that name.
    ///     OverflowError: `start` or `end` is negative.
    ///     TypeError: `start` or `end` is not an integer.
    ///     ValueError: `start` is 0, the range is inverted, or `end` is past the record.
    ///
    ///     >>> tbr.get_inclusive("chr2", 10001, 10010)
    ///     'CGTATCCCAC'
    #[rustfmt::skip]
    pub fn get_inclusive<'py>(&self, py: Python<'py>, name: &str, start: usize, end: usize) -> PyResult<Bound<'py, PyString>> {
        if start < 1 {
            return Err(PyValueError::new_err("Invalid start for inclusive range (start == 0)"));
        }
        self.get(py, name, start - 1, end)
    }

    /// A version of `get_batch` using 1-based inclusive ranges.
    ///
    /// See https://standage.github.io/on-genomic-interval-notation.html.
    ///
    /// Args:
    ///     batch: Iterable of `(name, start, end)` tuples, 1-based with inclusive end.
    ///
    /// Returns:
    ///     list[str]: One sequence per input range, in input order.
    ///
    /// Raises:
    ///     KeyError: An item names a record that does not exist.
    ///     OverflowError: An item has a negative index.
    ///     TypeError: An item is not a 3-tuple of `(str, int, int)`.
    ///     ValueError: An item has an invalid range, including a `start` of 0.
    #[rustfmt::skip]
    pub fn get_batch_inclusive<'py>(&self, py: Python<'py>, batch: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyList>> {
        self.get_batch_impl(py, batch, 1)
    }

    /// Concatenate several ranges of one sequence record into a single string.
    ///
    /// Args:
    ///     name: Sequence record name.
    ///     ranges: Iterable of `(start, end)` pairs, 0-based with exclusive end.
    ///
    /// Returns:
    ///     str: The concatenated nucleotides.
    ///
    /// Raises:
    ///     KeyError: No record has that name.
    ///     OverflowError: An item has a negative index.
    ///     TypeError: An item is not a 2-tuple of `(int, int)`.
    ///     ValueError: An item has an invalid range.
    ///
    ///     >>> tbr.concat("chr2", [(10000, 10006), (10006, 10010)])
    ///     'CGTATCCCAC'
    #[rustfmt::skip]
    pub fn concat<'py>(&self, py: Python<'py>, name: &str, ranges: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyString>> {
        self.concat_impl(py, name, ranges, 0)
    }

    /// A version of `concat` using 1-based inclusive ranges.
    ///
    /// Args:
    ///     name: Sequence record name.
    ///     ranges: Iterable of `(start, end)` pairs, 1-based with inclusive end.
    ///
    /// Returns:
    ///     str: The concatenated nucleotides.
    ///
    /// Raises:
    ///     KeyError: No record has that name.
    ///     OverflowError: An item has a negative index.
    ///     TypeError: An item is not a 2-tuple of `(int, int)`.
    ///     ValueError: An item has an invalid range, including a `start` of 0.
    ///
    ///     >>> tbr.concat_inclusive("chr2", [(10001, 10006), (10007, 10010)])
    ///     'CGTATCCCAC'
    pub fn concat_inclusive<'py>(
        &self,
        py: Python<'py>,
        name: &str,
        ranges: &Bound<'_, PyAny>,
    ) -> PyResult<Bound<'py, PyString>> {
        self.concat_impl(py, name, ranges, 1)
    }

    /// Ask the operating system to begin paging in the data for a batch of ranges.
    ///
    /// Use this method only if are reading a "cold" file that is not already paged in from disk.
    /// Accessing a cold file is almost entirely IO-bound. This function provides a hint to the
    /// operating system, specifying what data you'll soon be asking for. It does not affect the
    /// output of subsequent calls to `get` and `concat`, only their speed.
    /// For maximum benefit, prefetch the largest batch of ranges that you can in one call.
    ///
    /// Args:
    ///     batch: Iterable of `(name, start, end)` tuples, 0-based with exclusive end.
    ///
    /// Raises:
    ///     KeyError: An item names a record that does not exist.
    ///     OverflowError: An item has a negative index.
    ///     TypeError: An item is not a 3-tuple of `(str, int, int)`.
    ///     ValueError: An item has `start > end`.
    ///
    ///     >>> tbr.prefetch(args)          # Returns immediately.
    ///     >>> seqs = tbr.get_batch(args)  # Reads land as the data arrives.
    pub fn prefetch(&self, batch: &Bound<'_, PyAny>) -> PyResult<()> {
        self.prefetch_impl(batch, 0)
    }

    /// A version of `prefetch` using 1-based inclusive ranges.
    ///
    /// Args:
    ///     batch: Iterable of `(name, start, end)` tuples, 1-based with inclusive end.
    ///
    /// Raises:
    ///     KeyError: An item names a record that does not exist.
    ///     OverflowError: An item has a negative index.
    ///     TypeError: An item is not a 3-tuple of `(str, int, int)`.
    ///     ValueError: An item has a `start` of 0, or an inverted range.
    pub fn prefetch_inclusive(&self, batch: &Bound<'_, PyAny>) -> PyResult<()> {
        self.prefetch_impl(batch, 1)
    }
}

// Helpers shared by the 0-based and 1-based methods above. These live outside the
// #[pymethods] block so that PyO3 does not export them as undocumented Python methods.
impl TwobitReader {
    // Looks up the TwobitSequenceData record, or returns a suitable Python error type, rather
    // than letting the crate's implementation raise a PanicException.
    #[inline]
    fn seq_data(&self, name: &str) -> PyResult<&TwobitSequenceData> {
        self.inner.try_get_seq_data_by_name(name).ok_or_else(|| PyKeyError::new_err(name.to_string()))
    }

    // Implements prefetch for 0-based-exclusive and 1-based-inclusive ranges.
    //
    // Items are validated here and raised as a suitable Python error type, rather than letting
    // the crate's implementation raise a PanicException.
    fn prefetch_impl(&self, batch: &Bound<'_, PyAny>, base: usize) -> PyResult<()> {
        // Items are mapped to (PyBackedStr, usize, usize), since PyBackedStr implements the
        // AsRef<str> required by prefetch, yet only bumps the refcount on the string rather
        // than copying it.
        let ec = PyErrCapture::default();
        let ranges = batch.try_iter()?.map_while(|item| {
            let (name, start, end) = ec.ok(ec.ok(item)?.extract::<(PyBackedStr, usize, usize)>())?;
            ec.ok(self.check_prefetch_range(&name, start, end, base))?;
            Some((name, start, end))
        });
        if base == 1 {
            self.inner.prefetch_inclusive(ranges);
        } else {
            self.inner.prefetch(ranges);
        }
        ec.into_result()
    }

    // Returns a Python error for anything the crate's prefetch would panic on.
    #[inline]
    fn check_prefetch_range(&self, name: &str, start: usize, end: usize, base: usize) -> PyResult<()> {
        self.seq_data(name)?;
        check_start_inclusive(start, base)?;
        if start - base > end {
            return Err(PyValueError::new_err("invalid range (start > end)"));
        }
        Ok(())
    }

    #[rustfmt::skip]
    fn get_batch_impl<'py>(&self, py: Python<'py>, batch: &Bound<'_, PyAny>, base: usize) -> PyResult<Bound<'py, PyList>> {
        // Helper function to unpack each item and return the corresponding string from get()
        let get = |item: Result<Bound<'_, PyAny>, PyErr>| -> Result<Bound<'py, PyString>, PyErr> {
            let item = item?;
            let (name, start, end): (&str, usize, usize) = item.extract()?;
            if base == 1 {
                self.get_inclusive(py, name, start, end)
            } else {
                self.get(py, name, start, end)
            }
        };

        if let Ok(len) = batch.len() {
            // Batch is sized, so create a sized PyList to match, initialized with None.
            let list = PyList::new(py, repeat_n(PyNone::get(py), len))?;
            for (i, item) in batch.try_iter()?.enumerate() {
                list.set_item(i, get(item)?)?;
            }
            Ok(list)
        } else {
            // Batch is unsized, so create an empty PyList and append each item as we go.
            let list = PyList::empty(py);
            for item in batch.try_iter()? {
                list.append(get(item)?)?;
            }
            Ok(list)
        }
    }

    #[rustfmt::skip]
    fn concat_impl<'py>(&self, py: Python<'py>, name: &str, ranges: &Bound<'_, PyAny>, base: usize) -> PyResult<Bound<'py, PyString>> {
        // The implementation used for concat decodes, when possible, directly into a pre-sized PyUnicode buffer and returns,
        // eliding any conversion or copying from String to PyString. However, the extra complexity of doing so is only
        // about 10% faster than the simple implementation in the comment below, where PyO3 is leveraged for simplicity:
        //
        //    pub fn concat<'py>(&self, name: &str, ranges: &Bound<'_, PyAny>) -> String  {
        //        let range_iter = ranges.try_iter().expect("Expected ranges to be iterable.");
        //        self.inner.concat_iter(name, range_iter.map(|item| -> (usize, usize) {
        //            let item = item.expect("Could not get item.");
        //            item.extract().expect("Could not extract (usize, usize) item.")
        //        }))
        //    }
        //
        // Since this library is low-level, for a very stable file format, and is unlikely to change, the
        // faster-but-harder-to-maintain version below is still used, and likewise for several other methods.

        // Check if ranges is re-iterable or not, and branch implementations based on that.
        let mut range_iter = ranges.try_iter()?;
        if ranges.is(&range_iter) {
            // Ranges returned itself as the iterator. Assume ranges is only iterable
            // once (e.g., a generator). Grow the output buffer rather than pre-sizing it.
            let ec = PyErrCapture::default();
            let result = self.inner.concat_iter(
                name,
                ranges.try_iter()?.map_while(|item| {
                    // Convert to a 0-based exclusive range, checking start first so it cannot underflow.
                    let (start, end) = ec.ok(ec.ok(item)?.extract::<(usize, usize)>())?;
                    ec.ok(check_start_inclusive(start, base))?;
                    Some((start - base, end))
                }),
            );
            ec.check()?; // Return error if any.

            // Now that the size is known and the result is stored as a String, copy it into
            // a fresh PyUnicode object and return.
            let (unicode_guard, buf) = new_unicode_guarded(py, result.len())?;
            buf.write_copy_of_slice(result.as_bytes());
            Ok(as_bound_unicode(py, unicode_guard))
        } else {
            // Ranges returned a separate iterator object. Assume ranges can be re-iterated.
            // Use this fact to pre-size the output buffer and fill it progressively.
            let seq = self.seq_data(name)?;
            let total_len = range_iter.try_fold(0, |total, item| -> PyResult<usize> {
                let (start, end) = item?.extract()?;
                check_start_inclusive(start, base)?; // Check start >= base before subtracting
                check_range(seq, start - base, end)?;
                Ok(total + end - (start - base))
            })?;

            // Allocate new PyUnicode object and immediately guard it against leaking on return or panic.
            let (unicode_guard, mut buf) = new_unicode_guarded(py, total_len)?;

            // Decode each range into its respective slice of the output buffer. The ranges are
            // re-validated rather than trusted, because a second iteration of a user-supplied
            // object is not guaranteed to yield the same lengths as the total_len pass.
            for item in ranges.try_iter()? {
                let (start, end): (usize, usize) = item?.extract()?;
                check_start_inclusive(start, base)?; // Check start >= base before subtracting
                check_range(seq, start - base, end)?;
                let start = start - base;
                if end - start > buf.len() {
                    return Err(PyValueError::new_err("ranges grew between iterations"));
                }
                let (head, tail) = buf.split_at_mut(end - start);
                decode_and_fill_blocks(&self.inner.mmap, seq, start, head);
                buf = tail;
            }

            // Every byte must have been written, or the tail of the buffer would reach Python
            // uninitialized. This can only happen if re-iterating yielded a shorter sequence.
            if !buf.is_empty() {
                return Err(PyValueError::new_err("ranges shrank between iterations"));
            }

            // Release the guard now that buf is completely initialized, and return as PyString.
            Ok(as_bound_unicode(py, unicode_guard))
        }
    }
}

/// Fast DNA sequence extraction from 2bit files.
///
/// 2bit is a standard binary format in bioinformatics, documented at
/// http://genome.ucsc.edu/FAQ/FAQformat.html#format7
///
/// The motivation for this package is speed; see the Python timing results in the
/// Github repository for this package.
///
/// Extracting sequences is straightforward:
///
///     >>> from twobitreader_rs import TwobitReader
///     >>> tbr = TwobitReader.open("hg38.2bit")   # Human genome, build 38
///     >>> tbr.get("chr1", 10000, 10005)
///     'TAACC'
///
/// Concatenation assembles a spliced transcript from its exons:
///
///     >>> exons = [(1389575, 1391118),           # Exon 1 (start, end)
///     ...          (1394695, 1395603)]           # Exon 2 (start, end)
///     >>> transcript = tbr.concat("chr6", exons)
///
/// Cold files are an order of magnitude slower to read than files already in memory.
/// Prefetching recovers most of that without resorting to threads:
///
///     >>> tbr.prefetch(args)                     # Returns immediately; results unaffected.
///     >>> seqs = tbr.get_batch(args)             # Reads land as the data arrives.
///
/// Ranges are 0-based with an exclusive end, matching Python slicing and the BED format.
/// Every method has an `_inclusive` counterpart taking 1-based inclusive ranges, as used by
/// GFF/GTF and by genome browsers. See
/// https://standage.github.io/on-genomic-interval-notation.html
///
/// Threads: the extension is marked free-threading compatible, so on a free-threaded build
/// of CPython using `ThreadPoolExecutor` achieves significant parallelism.
//
// This module is safe for free-threaded Python distributions.
// - TwobitReader is Send + Sync and its interior mutability is synchronized via OnceLock.
// - Global state is const, never modified.
// - The GIL itself is never directly used.
// - The Cell used by PyErrCapture only appears as a local non-shared variable.
#[pymodule(gil_used = false)]
fn twobitreader_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<TwobitReader>()?;
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    Ok(())
}

/// Return the reverse complement of a DNA sequence.
///
/// Useful when assembling strand-sensitive transcripts, in combination with `concat`.
///
/// Args:
///     dna: A string containing only the characters `ACGTNacgtn`.
///
/// Returns:
///     str: The reverse complement.
///
/// Raises:
///     TypeError: `dna` is not a string.
///
///     >>> reverse_complement("ACGTN")
///     'NACGT'
///     >>> reverse_complement("acgtn")
///     'nacgt'
///
/// Note:
///     Characters outside `ACGTNacgtn` are replaced with `?` rather than reported, so
///     validate the input if needed for your application.
#[pyfunction]
fn reverse_complement(dna: String) -> PyResult<String> {
    Ok(crate::reverse_complement(dna))
}

#[inline]
#[allow(clippy::type_complexity)]
fn new_unicode_guarded<'py>(
    py: Python<'py>,
    len: usize,
) -> PyResult<(ScopeGuard<*mut PyObject, impl FnOnce(*mut PyObject)>, &'py mut [MaybeUninit<u8>])> {
    if len >= isize::MAX as usize {
        return Err(PyValueError::new_err("range length exceeded string size limit"));
    }

    // Allocate a new PyASCIIObject, and check that allocation succeeded.
    let unicode = unsafe { PyUnicode_New(len as isize, 127) };
    if unicode.is_null() {
        return Err(PyErr::fetch(py));
    }

    // Guard the pointer immediately, to ensure that Py_DECREF is called on early return
    // or on panic, unless the guard is explicitly released.
    // SAFETY: for PyASCIIObject, dealloc with uninitialized buffer is safe.
    let unicode_guard = scopeguard::guard(unicode, |o| {
        unsafe { Py_DECREF(o) };
    });

    // SAFETY: PyUnicode_1BYTE_DATA on a PyASCIIObject is guaranteed never to return NULL
    let buf = unsafe { slice::from_raw_parts_mut(PyUnicode_1BYTE_DATA(unicode).cast(), len) };
    Ok((unicode_guard, buf))
}

#[inline]
#[allow(clippy::type_complexity)]
fn as_bound_unicode<'py>(
    py: Python<'py>,
    unicode_guard: ScopeGuard<*mut PyObject, impl FnOnce(*mut PyObject)>,
) -> Bound<'py, PyString> {
    debug_assert!(!unicode_guard.is_null());
    debug_assert!(unsafe { PyUnicode_CheckExact(*unicode_guard) != 0 });

    // Release the guard and return the PyUnicode as bound to py.
    // SAFETY: the unchecked cast is safe, because this function is only called from within this crate
    // and every call passes exactly a PyUnicode instance.
    let unicode = ScopeGuard::into_inner(unicode_guard);
    unsafe { Bound::from_owned_ptr(py, unicode).cast_into_unchecked() }
}

// Returns a ValueError if the range is invalid or outside seq's range.
#[inline]
fn check_range(seq: &crate::TwobitSequenceData, start: usize, end: usize) -> PyResult<()> {
    if start > end {
        return Err(PyValueError::new_err("invalid range (start > end)"));
    }
    if end > seq.dna_len {
        return Err(PyValueError::new_err("invalid end (end > dna_len)"));
    }
    Ok(())
}

// Returns a ValueError if start is not valid for an inclusive range
#[inline]
fn check_start_inclusive(start: usize, base: usize) -> PyResult<()> {
    if start < base {
        debug_assert_eq!(base, 1, "expected base=1 but found {base}");
        return Err(PyValueError::new_err("invalid start (0) for 1-based range"));
    }
    Ok(())
}

// Decodes the DNA from seq at rang start..start+buf.len() into buf, including any block fills.
#[inline]
fn decode_and_fill_blocks(mmap: &Mmap, seq: &TwobitSequenceData, start: usize, buf: &mut [MaybeUninit<u8>]) {
    crate::decode_from_mmap(mmap, seq, start, buf);
    crate::fill_blocks(seq, start, unsafe { buf.assume_init_mut() });
}

// Captures the first error raised inside an iterator adaptor, so that a fallible closure can
// end the iteration instead of panicking, and the error can be returned once iteration is over.
//
// fn pure_rust_func<I, S>(items: I)
// where
//     I: IntoIterator<Item = S>,
//     S: AsRef<str>,
// {
//     for item in items.into_iter() {
//         println!("{}", item.as_ref());  // Print the str.
//     }
// }
//
// fn py_wrapper_func(items: &Bound<'_, PyAny>) -> PyResult<()> {
//     let ec = PyErrCapture::default();
//     pure_rust_func(items.try_iter()?.map_while(|item| {
//         let item = ec.ok(item)?;                          // Return None if error.
//         let item = ec.ok(item.extract::<PyBackedStr>());  // Capture error here, if any.
//         item
//     }));
//     ec.into_result()
// }
//
// Implemented as a Cell to allow nested ok() calls to pass the borrow check.
#[derive(Default)]
struct PyErrCapture(Cell<Option<PyErr>>);

impl PyErrCapture {
    // Converts Err into None, retaining the first error seen. Intended for map_while, where
    // returning None ends the iteration at the point of failure.
    fn ok<T>(&self, result: PyResult<T>) -> Option<T> {
        match result {
            Ok(value) => Some(value),
            Err(err) => {
                // Take-then-set, rather than set, so that a later error cannot displace the first.
                let first = self.0.take().or(Some(err));
                self.0.set(first);
                None
            }
        }
    }

    // Returns and clears the captured error. Intended as `pec.check()?` part way through a
    // function, so that code after the iteration is skipped once an error has occurred.
    fn check(&self) -> PyResult<()> {
        self.0.take().map_or(Ok(()), Err)
    }

    // Consuming form of check, for the end of a function.
    fn into_result(self) -> PyResult<()> {
        self.0.into_inner().map_or(Ok(()), Err)
    }
}
