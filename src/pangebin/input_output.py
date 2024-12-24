"""Input-Output PlasBin-Flow module."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import gzip
from typing import IO


def is_gz_file(filepath: Path) -> bool:
    """Check if a file is gzipped.

    Parameters
    ----------
    filepath : Path
        Path of file to check.

    Returns
    -------
    bool
        True if file is gzipped.

    """
    if not filepath.exists():
        return filepath.suffix == ".gz"
    with filepath.open("rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


def open_file_read(filepath: Path) -> IO[str]:
    """Open a file for reading.

    Parameters
    ----------
    filepath : Path
        Path of file to read.

    Returns
    -------
    file object
        File to read.

    """
    if is_gz_file(filepath):
        return gzip.open(filepath, "rt")
    return filepath.open()


def open_file_write(filepath: Path) -> IO[str]:
    """Open a file for writing.

    Parameters
    ----------
    filepath : Path
        Path of file to write to.

    Returns
    -------
    file object
        File to write to.

    """
    if is_gz_file(filepath):
        return gzip.open(filepath, "wt")
    return filepath.open("w")
