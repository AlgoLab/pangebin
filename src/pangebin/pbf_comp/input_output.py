"""PlasBin-flow compatibility input-output module."""

from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import _csv
    from collections.abc import Generator, Iterator
    from pathlib import Path
import csv
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import _csv
    from collections.abc import Generator, Iterator
    from pathlib import Path


class PBFPLMWriter:
    """PlasBin-flow plasmidness scores TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        file: Path,
    ) -> Generator[PBFPLMWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = PBFPLMWriter(csv.writer(f_out, delimiter="\t"))
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def write_sequence_plasmidness(self, sequence_id: str, plasmidness: float) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow([sequence_id, plasmidness])


class PBFPLMReader:
    """PlasBin-flow plasmidness scores TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[PBFPLMReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = PBFPLMReader(csv.reader(f_in, delimiter="\t"))
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader

    def __iter__(self) -> Iterator[tuple[str, float]]:
        """Iterate sequence probability scores."""
        return ((row[0], float(row[1])) for row in self.__csv_reader)


class PBFSeedWriter:
    """PlasBin-flow seed sequences TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        file: Path,
    ) -> Generator[PBFSeedWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = PBFSeedWriter(csv.writer(f_out, delimiter="\t"))
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def write_sequence(self, sequence_id: str) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow([sequence_id])


class PBFSeedReader:
    """PlasBin-flow seed sequences TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[PBFSeedReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = PBFSeedReader(csv.reader(f_in, delimiter="\t"))
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader

    def __iter__(self) -> Iterator[str]:
        """Iterate sequence probability scores."""
        return (row[0] for row in self.__csv_reader)
