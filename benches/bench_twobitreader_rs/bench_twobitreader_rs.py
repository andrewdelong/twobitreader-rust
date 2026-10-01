"""Benchmarks for the twobitreader_rs Python bindings.

This mirrors the Rust benchmark in benches/bench_twobitreader/main.rs -- same workloads, same
cache scenarios, same timing harness and same output format -- so the two sets of numbers can be
compared directly. The _threads variants correspond to the Rust _rayon variants.

The _threads variants only scale on a free-threaded interpreter, since a stock build holds the
GIL for the duration of every binding call. Whether the GIL is active is reported on startup.

Run from the repository root:

    python benches/bench_twobitreader_rs/bench_twobitreader_rs.py [name filter ...]
"""

# Standard library
import gzip
import os
import mmap
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from enum import Enum
from functools import partial
from pathlib import Path
from time import perf_counter, sleep

# Dependencies
from twobitreader_rs import TwobitReader, reverse_complement

# When a bench callback is invoked several times, only this many will
# be printed to console describing the result (for legibility).
NUM_BENCH_TIMES_PRINTED = 5

HG38_PATH = Path("tests/assets/hg38.p13.2bit")
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.2bit"

# Path of the scratch copy used by cold benchmarks.
COLD_PATH = Path("tests/assets/hg38.p13.cold.2bit")

# fcntl.F_NOCACHE is not exposed by Python's fcntl module; 48 is its value in <sys/fcntl.h>.
F_NOCACHE = 48


class Cache(Enum):
    """Which cache scenario to simulate in the benchmark run."""

    HOT = "hot"  # File already resides in the page cache.
    COLD = "cold"  # File is not yet in the page cache at all.
    PREFETCH = "prefetch"  # Cold, but with a prefetch hint at upcoming data access.


class Parallel(Enum):
    """How and whether to use parallelism in the benchmark run."""

    NONE = "none"  # No parallelism.
    THREADS = "threads"  # Use a thread pool, the counterpart of the Rust rayon runs.


# One worker per core, mirroring rayon's default pool size.
NUM_THREADS = os.cpu_count() or 1

# Whether this interpreter still has a GIL, which serializes every binding call and so makes
# the _threads runs pointless. sys._is_gil_enabled exists only on 3.13+; older builds always have one.
GIL_ENABLED = getattr(sys, "_is_gil_enabled", lambda: True)()

# Several chunks per worker, so the pool has some room to balance uneven work in the
# way rayon's work stealing does. Chunking also amortises the per-call binding overhead.
CHUNKS_PER_THREAD = 4

_EXECUTOR = None


def executor():
    """Return the shared thread pool, creating it on first use.

    The pool is built once and reused, so that thread creation never falls inside a timed
    region. This mirrors the Rust harness, which builds rayon's global pool once in main.
    """
    global _EXECUTOR
    if _EXECUTOR is None:
        _EXECUTOR = ThreadPoolExecutor(max_workers=NUM_THREADS)
    return _EXECUTOR


def chunked(items):
    """Split items into contiguous chunks, one work item each for the thread pool."""
    size = max(1, len(items) // (NUM_THREADS * CHUNKS_PER_THREAD))
    return [items[i : i + size] for i in range(0, len(items), size)]


def bench(name, f, num_iter, unit="ms"):
    """Invoke f() num_iter times and print the mean/median/min/max of the reported durations.

    The durations are reported as integers in the specified unit, matching the Rust harness.
    """
    # A threaded run on a GIL-enabled interpreter only reproduces the serial time, several
    # minutes at a time, so skip it rather than report a number that invites misreading.
    if GIL_ENABLED and name.endswith("_threads"):
        return

    # Any command line arguments act as filters, so that a subset can be re-run quickly.
    # Flags are skipped, to match the Rust harness where cargo passes --bench through.
    filters = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    if filters and not any(filter in name for filter in filters):
        return

    # Convert a duration in seconds to the requested unit, truncating as Duration::as_* does.
    scale = {"ns": 1e9, "us": 1e6, "ms": 1e3, "s": 1}
    if unit not in scale:
        raise ValueError(f"invalid unit '{unit}'")

    # Invoke the callback num_iter times, recording the reported duration each time.
    times = []
    for _ in range(num_iter):
        # Sleep 0.1 seconds prior to each run to reduce variance between invocations
        # (Allows OS to quiet down any console output, etc.)
        sleep(0.1)

        # Run the function and record its reported wall-clock time.
        times.append(int(f() * scale[unit]))

    # Compute summary stats from the times.
    avg_time = sum(times) // len(times)
    med_time = sorted(times)[len(times) // 2]
    min_time = min(times)
    max_time = max(times)

    # Format the output string, in a way that also reports the first few measurements.
    avg_str = f"{avg_time}"
    times_str = ""
    stats_str = ""
    if len(times) > 1:
        avg_str = f"avg={avg_str}"
        stats_str = f"(med={med_time} min={min_time} max={max_time})"
        times_str = f"{times[:NUM_BENCH_TIMES_PRINTED]}"
        if len(times) > NUM_BENCH_TIMES_PRINTED:
            times_str = times_str.replace("]", ", ...]")

    # Print bench result indented, and format the test name in green.
    print(f"     \x1b[32m{name}\x1b[0m: {avg_str} {unit} {stats_str} {times_str}")


def read_bed(path):
    """Read the first 6 columns of a bed.gz file, as a list of (chrom, start, end, strand, name)."""
    with gzip.open(path, "rt", encoding="ascii") as f:
        rows = []
        for line in f:
            cols = line.split("\t", 6)
            chrom, start, end, name, _score, strand = cols[:6]
            rows.append((chrom, int(start), int(end), strand.rstrip()[0], name))
    return rows


def read_exons():
    """Read gencode-exons.bed.gz intervals, as a list of (chrom, start, end)."""
    bed = read_bed("benches/assets/gencode-exons.bed.gz")
    return [(chrom, start, end) for (chrom, start, end, _strand, _name) in bed]


def read_transcripts():
    """Read gencode-transcripts.bed.gz, grouped by transcript ID (the BED name field).

    Returns a list of (transcript_id, chrom, strand, [(start, end), ...]).
    """
    bed = read_bed("benches/assets/gencode-transcripts.bed.gz")
    lookup = defaultdict(list)
    for chrom, start, end, strand, id in bed:
        lookup[(id, chrom, strand)].append((start, end))
    return [(id, chrom, strand, exons) for ((id, chrom, strand), exons) in lookup.items()]


def load_into_page_cache(path):
    """Map the given file and access all pages, so that it resides in the page cache.

    If there is sufficient unused physical memory, the operating system will keep
    the data in memory (hot) for subsequent reuse.
    """
    with open(path, "rb") as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            # Touching one byte per page is enough to fault the whole mapping in, and avoids
            # copying the file. (Avoid madvise because it is not portable.)
            mm[:: mmap.PAGESIZE]


def make_cold_copy(src):
    """Write a copy of the 2bit file whose contents are not in the page cache, and return its path.

    This is how a "cold" benchmark is set up without `sudo purge`, which would evict every other
    file on the machine as well. Caching is disabled on the destination descriptor, so the bytes
    written here go straight to disk and are never cached; truncating the file first discards
    whatever a previous iteration of the benchmark left behind. Reading the copy back therefore
    has to go to disk, exactly as reading a freshly downloaded 2bit file would.

    Returns None on targets where the page cache cannot be bypassed per file.
    """
    if sys.platform != "darwin":
        return None
    import fcntl

    with open(src, "rb") as fsrc, open(COLD_PATH, "wb") as fdst:
        fcntl.fcntl(fdst.fileno(), F_NOCACHE, 1)
        while chunk := fsrc.read(8 << 20):
            fdst.write(chunk)
        fdst.flush()
        os.fsync(fdst.fileno())
    return COLD_PATH


def download_hg38():
    """Return a local path to the hg38 2bit file, downloading it from UCSC if not present."""
    if not HG38_PATH.exists():
        # --fail and --location so an HTTP error or redirect is not written out as if it
        # were the genome; the Rust harness omits these but is exposed to the same hazard.
        subprocess.run(["curl", "--fail", "--location", "-o", str(HG38_PATH), HG38_URL], check=True)
    return HG38_PATH


def setup(cache):
    """Return the path to benchmark against, preparing the page cache as the scenario requires.

    Setup happens outside the timed region, so its cost is not part of any reported duration.
    """
    path = download_hg38()
    if cache is Cache.HOT:
        load_into_page_cache(path)
        return path
    return make_cold_copy(path)


def apply_strand(dna, strand):
    if strand == "-":
        return reverse_complement(dna)
    assert strand == "+", f"Invalid strand '{strand}'"
    return dna


def bench_hg38_open(cache):
    """Load hg38 and time processing the header."""
    path = setup(cache)
    tic = perf_counter()

    # Open file, either the hot version or the cold version.
    tbr = TwobitReader.open(str(path))

    duration = perf_counter() - tic
    del tbr
    return duration


def bench_hg38_exons_prefetch_only():
    """Time how long prefetch takes, without waiting for data to finish paging in.

    The goal is to estimate how much time in the prefetch cases is taken by the prefetch phase.
    """
    path = setup(Cache.COLD)
    exons = read_exons()
    tic = perf_counter()

    tbr = TwobitReader.open(str(path))
    tbr.prefetch(exons)

    return perf_counter() - tic


def bench_hg38_exons(cache, parallel):
    """Load hg38 exon intervals and time how long it takes to extract their DNA."""
    path = setup(cache)
    exons = read_exons()
    tic = perf_counter()

    # Open file and collect strings into a list, either across threads or sequentially.
    tbr = TwobitReader.open(str(path))
    if cache is Cache.PREFETCH:
        tbr.prefetch(exons)
    if parallel is Parallel.THREADS:
        dst = [seq for part in executor().map(tbr.get_batch, chunked(exons)) for seq in part]
    else:
        dst = tbr.get_batch(exons)

    duration = perf_counter() - tic
    del dst
    return duration


def concat_transcripts(tbr, transcripts):
    """Extract and strand-correct the DNA of each transcript in the given list."""
    return [apply_strand(tbr.concat(chrom, exons), strand) for (_id, chrom, strand, exons) in transcripts]


def bench_hg38_transcripts(cache, parallel):
    """Load hg38 transcript intervals (groups of exons) and time extracting their DNA."""
    path = setup(cache)
    transcripts = read_transcripts()
    tic = perf_counter()

    # Open file and collect transcript strings into a list, either across threads or sequentially.
    tbr = TwobitReader.open(str(path))
    if cache is Cache.PREFETCH:
        # Collect the exons of all transcripts.
        tbr.prefetch(
            (chrom, start, end) for (_id, chrom, _strand, exons) in transcripts for (start, end) in exons
        )
    if parallel is Parallel.THREADS:
        parts = executor().map(partial(concat_transcripts, tbr), chunked(transcripts))
        dst = [seq for part in parts for seq in part]
    else:
        dst = concat_transcripts(tbr, transcripts)

    duration = perf_counter() - tic
    del dst
    return duration


def main():
    # The _threads runs only scale without the GIL, so make its state part of the output.
    if GIL_ENABLED:
        print("     GIL enabled: skipping the _threads benchmarks")
    else:
        print(f"     GIL disabled: running the _threads benchmarks on {NUM_THREADS} threads")

    # Hot: the file is already in the page cache, so these measure decoding speed.
    bench("hg38_open_hot", lambda: bench_hg38_open(Cache.HOT), 3)
    bench("hg38_exons_hot_basic", lambda: bench_hg38_exons(Cache.HOT, Parallel.NONE), 10)
    bench("hg38_exons_hot_threads", lambda: bench_hg38_exons(Cache.HOT, Parallel.THREADS), 10)
    bench("hg38_transcripts_hot_basic", lambda: bench_hg38_transcripts(Cache.HOT, Parallel.NONE), 10)
    bench("hg38_transcripts_hot_threads", lambda: bench_hg38_transcripts(Cache.HOT, Parallel.THREADS), 10)

    # Cold: the file must be read from disk, so these measure how well the reads are overlapped.
    if make_cold_copy(download_hg38()) is not None:
        bench("hg38_open_cold", lambda: bench_hg38_open(Cache.COLD), 3)
        bench("hg38_exons_cold_prefetch_only", bench_hg38_exons_prefetch_only, 10)
        bench("hg38_exons_cold_prefetch", lambda: bench_hg38_exons(Cache.PREFETCH, Parallel.NONE), 10)
        bench("hg38_exons_cold_basic", lambda: bench_hg38_exons(Cache.COLD, Parallel.NONE), 10)
        bench("hg38_exons_cold_threads", lambda: bench_hg38_exons(Cache.COLD, Parallel.THREADS), 10)
        bench("hg38_transcripts_cold_prefetch", lambda: bench_hg38_transcripts(Cache.PREFETCH, Parallel.NONE), 10)
        bench("hg38_transcripts_cold_basic", lambda: bench_hg38_transcripts(Cache.COLD, Parallel.NONE), 10)
        bench("hg38_transcripts_cold_threads", lambda: bench_hg38_transcripts(Cache.COLD, Parallel.THREADS), 10)
        COLD_PATH.unlink()
    else:
        print("     skipping cold benchmarks: cannot write an uncached file on this target")

    if _EXECUTOR is not None:
        _EXECUTOR.shutdown()


if __name__ == "__main__":
    main()
