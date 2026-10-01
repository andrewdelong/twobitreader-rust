"""Type stubs for the twobitreader_rs extension module.

Signatures only. The full documentation, including Args/Returns/Raises for each method,
lives in the runtime docstrings generated from src/python.rs and is visible via help().
"""

from collections.abc import Iterable
from typing import Never, final

__all__ = ["TwobitReader", "reverse_complement"]

@final
class TwobitReader:
    """A reader for a single 2bit file.

    Not constructible directly; use TwobitReader.open or TwobitReader.open_masked.
    """

    # The extension type has no constructor, so calling it raises TypeError. An
    # unsatisfiable keyword makes a type checker report that at the call site.
    def __init__(self, *, _use_open_or_open_masked: Never) -> None: ...

    @staticmethod
    def open(path: str) -> TwobitReader:
        """Open a 2bit file for reading."""

    @staticmethod
    def open_masked(path: str) -> TwobitReader:
        """Open a 2bit file for reading, with the lowercase mask applied."""

    def names(self) -> list[str]:
        """Return the sequence record names, in the order they appear in the file."""

    def contains_name(self, name: str) -> bool:
        """Return whether `name` matches a sequence record in the file."""

    def seq_len(self, name: str) -> int:
        """Return the length in nucleotides of the named sequence record."""

    def get(self, name: str, start: int, end: int) -> str:
        """Extract the range `start..end` (0-based, exclusive end) from a record."""

    def get_inclusive(self, name: str, start: int, end: int) -> str:
        """Extract a 1-based inclusive range from a record."""

    def get_batch(self, batch: Iterable[tuple[str, int, int]]) -> list[str]:
        """Extract many `(name, start, end)` ranges in one call."""

    def get_batch_inclusive(self, batch: Iterable[tuple[str, int, int]]) -> list[str]:
        """Extract many 1-based inclusive ranges in one call."""

    def concat(self, name: str, ranges: Iterable[tuple[int, int]]) -> str:
        """Concatenate several `(start, end)` ranges of one record into one string."""

    def concat_inclusive(self, name: str, ranges: Iterable[tuple[int, int]]) -> str:
        """Concatenate several 1-based inclusive ranges of one record."""

    def prefetch(self, batch: Iterable[tuple[str, int, int]]) -> None:
        """Ask the operating system to begin paging in the data for a batch of ranges."""

    def prefetch_inclusive(self, batch: Iterable[tuple[str, int, int]]) -> None:
        """Ask the operating system to begin paging in a batch of 1-based ranges."""

def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence containing only ACGTNacgtn."""
