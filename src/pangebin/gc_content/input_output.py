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


class IntervalStepsHeader(StrEnum):
    """Interval steps TSV header."""

    INTERVAL_STEP = "Interval_step"


class IntervalStepsWriter:
    """Interval steps TSV writer."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[IntervalStepsWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = IntervalStepsWriter(csv.writer(f_out, delimiter="\t"))
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer
        self.__write_header()

    def __write_header(self) -> None:
        """Write header."""
        self.__csv_writer.writerow([IntervalStepsHeader.INTERVAL_STEP])

    def write_interval_step(self, interval_step: float) -> None:
        """Write interval step."""
        self.__csv_writer.writerow([interval_step])


class IntervalStepsReader:
    """Interval steps TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[IntervalStepsReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = IntervalStepsReader(csv.reader(f_in, delimiter="\t"))
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader
        next(self.__csv_reader)  # skip header

    def __iter__(self) -> Iterator[float]:
        """Iterate over the interval steps."""
        for interval_step in self.__csv_reader:
            yield float(interval_step[0])


class ScoresHeader(StrEnum):
    """GC scores TSV header."""

    SEQUENCE_ID = "Sequence_ID"


class ScoresWriter:
    """GC scores TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        intervals: items.Intervals,
        file: Path,
    ) -> Generator[ScoresWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = ScoresWriter(csv.writer(f_out, delimiter="\t"))
            writer.__write_header(intervals)
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def __write_header(self, intervals: items.Intervals) -> None:
        """Write header."""
        self.__csv_writer.writerow(
            [
                ScoresHeader.SEQUENCE_ID,
                *(items.IntervalFormatter.to_str(interval) for interval in intervals),
            ],
        )

    def write_sequence_gc_scores(self, seq_gc_scores: items.SequenceGCScores) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow(
            [
                seq_gc_scores.sequence_id(),
                *seq_gc_scores.probability_scores(),
            ],
        )


class ScoresReader:
    """GC scores TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[ScoresReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = ScoresReader(csv.reader(f_in, delimiter="\t"))
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader
        self.__intervals: items.Intervals = self.__header_to_intervals()

    def intervals(self) -> items.Intervals:
        """Get intervals."""
        return self.__intervals

    def __iter__(self) -> Iterator[items.SequenceGCScores]:
        """Iterate sequence GC scores."""
        return (
            items.SequenceGCScores(
                row[0],
                (float(item) for item in row[1:]),
            )
            for row in self.__csv_reader
        )

    def __header_to_intervals(self) -> items.Intervals:
        """Parse header."""
        intervals = items.Intervals()
        intervals_details = [
            items.IntervalFormatter.from_str(item)
            for item in next(self.__csv_reader)[1:]
        ]
        for interval in intervals_details:
            intervals.append(interval[0])
        intervals.append(intervals_details[-1][1])
        return intervals
