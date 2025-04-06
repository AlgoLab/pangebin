"""Bin input/output module."""

from __future__ import annotations

import csv
from contextlib import contextmanager
from enum import StrEnum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import _csv
    from collections.abc import Generator, Iterator
    from pathlib import Path


class Header(StrEnum):
    """Bin sequences normalized coverage TSV header."""

    SEQUENCE_ID = "Sequence_ID"
    NORMALIZED_COVERAGE = "Normalized_coverage"


class Writer:
    """Bin sequences normalized coverage TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        file: Path,
    ) -> Generator[Writer, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = Writer(csv.writer(f_out, delimiter="\t"))
            writer.__write_header()
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def __write_header(self) -> None:
        """Write header."""
        self.__csv_writer.writerow([Header.SEQUENCE_ID, Header.NORMALIZED_COVERAGE])

    def write_sequence_normcov(self, sequence_id: str, normcov: float) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow([sequence_id, normcov])


class Reader:
    """Bin sequences normalized coverage TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[Reader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = Reader(csv.reader(f_in, delimiter="\t"))
            next(reader.__csv_reader)  # skip header
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader

    def __iter__(self) -> Iterator[tuple[str, float]]:
        """Iterate sequence probability scores."""
        return ((row[0], float(row[1])) for row in self.__csv_reader)
