"""GC content input-output module."""

from __future__ import annotations

import csv
from contextlib import contextmanager
from enum import StrEnum
from typing import TYPE_CHECKING

from pangebin.gc_content import items

if TYPE_CHECKING:
    import _csv
    from collections.abc import Generator, Iterator
    from pathlib import Path


class Header(StrEnum):
    """GC probability scores TSV header."""

    SEQUENCE_ID = "sequence_id"

    @classmethod
    def fmt_interval(cls, interval: tuple[float, float]) -> str:
        """Format interval."""
        return f"{interval[0]}_{interval[1]}"

    @classmethod
    def str_to_interval(cls, interval_str: str) -> tuple[float, float]:
        """Parse interval."""
        l_interval_str = interval_str.split("_")
        return (float(l_interval_str[0]), float(l_interval_str[1]))


class Writer:
    """GC probability scores TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        intervals: items.Intervals,
        file: Path,
    ) -> Generator[Writer, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = Writer(csv.writer(f_out, delimiter="\t"))
            writer.__write_header(intervals)
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def __write_header(self, intervals: items.Intervals) -> None:
        """Write header."""
        self.__csv_writer.writerow(
            [
                Header.SEQUENCE_ID,
                *(Header.fmt_interval(interval) for interval in intervals),
            ],
        )

    def write_sequence_proba_scores(
        self,
        seq_proba_scores: items.SequenceProbabilityScores,
    ) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow(
            [
                seq_proba_scores.sequence_id(),
                *seq_proba_scores.probability_scores(),
            ],
        )


class Reader:
    """GC probability scores TSV reader."""

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
        self.__intervals: items.Intervals = self.__header_to_intervals()

    def intervals(self) -> items.Intervals:
        """Get intervals."""
        return self.__intervals

    def __iter__(self) -> Iterator[items.SequenceProbabilityScores]:
        """Iterate sequence probability scores."""
        for row in self.__csv_reader:
            yield items.SequenceProbabilityScores(
                row[0],
                (float(item) for item in row[1:]),
            )

    def __header_to_intervals(self) -> items.Intervals:
        """Parse header."""
        intervals = items.Intervals()
        intervals_details = [
            Header.str_to_interval(item) for item in next(self.__csv_reader)[1:]
        ]
        for interval in intervals_details:
            intervals.append(interval[0])
        intervals.append(intervals_details[-1][1])
        return intervals
