"""Bin input/output module."""

from __future__ import annotations

import csv
from contextlib import contextmanager
from enum import StrEnum
from typing import TYPE_CHECKING

from . import items

if TYPE_CHECKING:
    import _csv
    from collections.abc import Generator, Iterable, Iterator
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

    def write_sequence_normcov(self, seq_normcov: items.SequenceNormCoverage) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow(
            [seq_normcov.identifier(), seq_normcov.normalized_coverage()],
        )

    def write_bunch_sequences_normcov(
        self,
        bunch_sequences_normcov: Iterable[items.SequenceNormCoverage],
    ) -> None:
        """Write bunch sequences normalized coverage."""
        for seq_normcov in bunch_sequences_normcov:
            self.write_sequence_normcov(seq_normcov)


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

    def __iter__(self) -> Iterator[items.SequenceNormCoverage]:
        """Iterate sequence probability scores."""
        return (
            items.SequenceNormCoverage(row[0], float(row[1]))
            for row in self.__csv_reader
        )
