"""Input-Output PlasBin-Flow module."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

import gzip
import io
from contextlib import contextmanager
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


@contextmanager
def open_file_read(filepath: Path) -> Generator[IO[str], None, None]:
    """Open a file for reading.

    Parameters
    ----------
    filepath : Path
        Path of file to read.

    Yields
    ------
    file object
        File to read.

    """
    if is_gz_file(filepath):
        with (
            filepath.open("rb") as f_in,
            gzip.GzipFile(fileobj=f_in, mode="rb") as g_in,
        ):
            text_in = io.TextIOWrapper(g_in, encoding="utf-8")
            yield text_in
            text_in.close()
    else:
        text_in = filepath.open()
        yield text_in
        text_in.close()


@contextmanager
def open_file_write(filepath: Path) -> Generator[IO[str], None, None]:
    """Open a file for writing.

    Parameters
    ----------
    filepath : Path
        Path of file to write to.

    Yields
    ------
    file object
        File to write to.

    """
    if is_gz_file(filepath):
        with (
            filepath.open("wb") as f_out,
            gzip.GzipFile(fileobj=f_out, mode="wb") as g_out,
        ):
            text_out = io.TextIOWrapper(g_out, encoding="utf-8")
            yield text_out
            text_out.close()
    else:
        text_out = filepath.open("w")
        yield text_out
        text_out.close()
