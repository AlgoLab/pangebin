"""Main input-output module."""

from __future__ import annotations

import datetime
import gzip
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

import io
import logging
import shutil
from contextlib import contextmanager
from typing import IO

import bgzip  # type: ignore[import-untyped]

_LOGGER = logging.getLogger(__name__)


def is_gz_file(filepath: Path) -> bool:
    """Check if a file is gzipped or bgzipped.

    The last extension must be ".gz"

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
        return filepath.suffix == ".gz" and test_f.read(2) == b"\x1f\x8b"


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
    text_in: IO[str]
    if is_gz_file(filepath):
        # XXX bug https://github.com/DataBiosphere/bgzip/issues/88
        # with (
        #     filepath.open("rb") as f_in,
        #     bgzip.BGZipReader(f_in) as g_in,
        # ):
        #     text_in = io.TextIOWrapper(g_in, encoding="utf-8")
        #     yield text_in
        #     text_in.close()
        with gzip.open(filepath, "rb") as f_in:
            text_in = io.TextIOWrapper(f_in, encoding="utf-8")
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
    text_out: IO[str]
    if is_gz_file(filepath):
        with (
            filepath.open("wb") as f_out,
            bgzip.BGZipWriter(fileobj=f_out) as g_out,
        ):
            text_out = io.TextIOWrapper(g_out, encoding="utf-8")
            yield text_out
            text_out.close()
    else:
        text_out = filepath.open("w")
        yield text_out
        text_out.close()


@contextmanager
def open_file_append(filepath: Path) -> Generator[IO[str], None, None]:
    """Open a file for appending.

    Parameters
    ----------
    filepath : Path
        Path of file to append to.

    Yields
    ------
    file object
        File to append to.

    """
    text_out: IO[str]
    if is_gz_file(filepath):
        with (
            filepath.open("ab") as f_out,
            bgzip.BGZipWriter(fileobj=f_out) as g_out,
        ):
            text_out = io.TextIOWrapper(g_out, encoding="utf-8")
            yield text_out
            text_out.close()
    else:
        text_out = filepath.open("a")
        yield text_out
        text_out.close()


@contextmanager
def possible_tmp_file(
    in_filepath: Path,
    out_filepath: Path | None = None,
) -> Generator[tuple[Path, Path], None, None]:
    """Possibly rename input file to replace it after.

    Parameters
    ----------
    in_filepath : Path
        Path of input file
    out_filepath : Path | None, optional
        Path of output file, by default None

    Yields
    ------
    Path
        Path of the input temporary file
    Path
        Path of the output file

    """
    replace_file = out_filepath is None

    if out_filepath is None:
        out_filepath = in_filepath
        in_filepath = in_filepath.rename(
            in_filepath.parent
            / (f"{datetime.datetime.now(tz=datetime.UTC).isoformat()}.gfa"),
        )
        _LOGGER.debug("Temporary input file: %s", in_filepath)
    elif in_filepath == out_filepath:
        _err_msg = f"Input and output files are the same: {in_filepath}"
        _LOGGER.error(_err_msg)
        raise ValueError(_err_msg)

    yield in_filepath, out_filepath

    if replace_file:
        _LOGGER.debug("Remove temporary input file: %s", in_filepath)
        in_filepath.unlink()


def bgzip_file(input_filepath: Path, compressed_filepath: Path | None = None) -> Path:
    """Compress a file using bgzip.

    Parameters
    ----------
    input_filepath : Path
        Path of file to compress.
    compressed_filepath : Path
        Path of compressed file.

    Returns
    -------
    Path
        Path of compressed file

    """
    if compressed_filepath is None:
        _LOGGER.info("Compressing file: %s", input_filepath)
        compressed_filepath = input_filepath.parent / f"{input_filepath.name}.gz"
        _LOGGER.debug("Output file: %s", compressed_filepath)
    else:
        if input_filepath == compressed_filepath:
            _err_msg = f"Input and output files are the same: {input_filepath}"
            _LOGGER.error(_err_msg)
            raise ValueError(_err_msg)
        _LOGGER.info("Compressing file: %s to %s", input_filepath, compressed_filepath)

    with (
        input_filepath.open("rb") as raw_in,
        compressed_filepath.open("wb") as raw_out,
        bgzip.BGZipWriter(raw_out) as f_out,
    ):
        shutil.copyfileobj(raw_in, f_out)

    return compressed_filepath


def bgunzip_file(
    compressed_filepath: Path,
    uncompressed_filepath: Path | None = None,
) -> Path:
    """Uncompress a file using bgzip.

    Parameters
    ----------
    compressed_filepath : Path
        Path of compressed file.
    uncompressed_filepath : Path, optional
        Path of uncompressed file, by default None.

    Returns
    -------
    Path
        Path of uncompressed file

    """
    if uncompressed_filepath is None:
        _LOGGER.info("Uncompressing file: %s", compressed_filepath)
        uncompressed_filepath = (
            compressed_filepath.parent / f"{compressed_filepath.name}.gz"
        )
        _LOGGER.debug("Output file: %s", uncompressed_filepath)
    else:
        if compressed_filepath == uncompressed_filepath:
            _err_msg = f"Input and output files are the same: {compressed_filepath}"
            _LOGGER.error(_err_msg)
            raise ValueError(_err_msg)
        _LOGGER.info(
            "Uncompressing file: %s to %s",
            compressed_filepath,
            uncompressed_filepath,
        )

    with (
        compressed_filepath.open("rb") as raw_in,
        uncompressed_filepath.open("wb") as raw_out,
        bgzip.BGZipReader(raw_in) as f_in,
    ):
        shutil.copyfileobj(f_in, raw_out)

    return uncompressed_filepath
