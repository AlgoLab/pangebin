"""PlasBin-flow compatibility input-output module."""

from __future__ import annotations

import csv
from contextlib import contextmanager
from enum import StrEnum
from typing import TYPE_CHECKING

import pangebin.pbf_comp.items as pbf_items

if TYPE_CHECKING:
    import _csv
    from collections.abc import Generator, Iterator
    from pathlib import Path


class PlmWriter:
    """PlasBin-flow plasmidness scores TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        file: Path,
    ) -> Generator[PlmWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = PlmWriter(csv.writer(f_out, delimiter="\t"))
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def write_sequence_plasmidness(self, sequence_id: str, plasmidness: float) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow([sequence_id, plasmidness])


class PlmReader:
    """PlasBin-flow plasmidness scores TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[PlmReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = PlmReader(csv.reader(f_in, delimiter="\t"))
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader

    def __iter__(self) -> Iterator[tuple[str, float]]:
        """Iterate sequence probability scores."""
        return ((row[0], float(row[1])) for row in self.__csv_reader)


class SeedWriter:
    """PlasBin-flow seed sequences TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        file: Path,
    ) -> Generator[SeedWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = SeedWriter(csv.writer(f_out, delimiter="\t"))
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def write_sequence(self, sequence_id: str) -> None:
        """Write sequence probability scores."""
        self.__csv_writer.writerow([sequence_id])


class SeedReader:
    """PlasBin-flow seed sequences TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[SeedReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = SeedReader(csv.reader(f_in, delimiter="\t"))
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader

    def __iter__(self) -> Iterator[str]:
        """Iterate sequence probability scores."""
        return (row[0] for row in self.__csv_reader)


class BinsHeader(StrEnum):
    """PlasBin-flow bins header."""

    PLASMID_ID = "#Pls_ID"
    FLOW = "Flow"
    GC_INTERVAL = "GC_bin"
    CONTIGS = "Contigs"


class IntervalFormatter:
    """PlasBin-flow interval formatter."""

    SEP = "-"

    @classmethod
    def to_str(cls, interval: tuple[float, float]) -> str:
        """Format interval."""
        return f"{interval[0]}{cls.SEP}{interval[1]}"

    @classmethod
    def from_str(cls, interval_str: str) -> tuple[float, float]:
        """Parse interval."""
        l_split = interval_str.split(cls.SEP)
        return (float(l_split[0]), float(l_split[1]))


class ContigMultFormatter:
    """PlasBin-flow contig multiplier formatter."""

    _SEP = ":"

    @classmethod
    def to_str(cls, contig_mult: pbf_items.ContigMult) -> str:
        """Format contig multiplier."""
        return f"{contig_mult.identifier()}{cls._SEP}{contig_mult.multiplicity():.2f}"

    @classmethod
    def from_str(cls, contig_mult_str: str) -> pbf_items.ContigMult:
        """Parse contig multiplier."""
        l_split = contig_mult_str.split(cls._SEP)
        return pbf_items.ContigMult(l_split[0], float(l_split[1]))


class BinsWriter:
    """PlasBin-flow bins TSV writer."""

    @classmethod
    @contextmanager
    def open(
        cls,
        file: Path,
    ) -> Generator[BinsWriter, None, None]:
        """Open TSV file for writing."""
        with file.open("w") as f_out:
            writer = BinsWriter(csv.writer(f_out, delimiter="\t"))
            writer.__write_header()
            yield writer

    def __init__(self, csv_writer: _csv._writer) -> None:
        """Initialize object."""
        self.__csv_writer = csv_writer

    def __write_header(self) -> None:
        """Write header."""
        self.__csv_writer.writerow(
            [
                BinsHeader.PLASMID_ID,
                BinsHeader.FLOW,
                BinsHeader.GC_INTERVAL,
                BinsHeader.CONTIGS,
            ],
        )

    def write_bin_line(
        self,
        bin_info: pbf_items.PBFBinInfo,
    ) -> None:
        """Write bin line."""
        self.__csv_writer.writerow(
            [
                bin_info.identifier(),
                bin_info.flow(),
                IntervalFormatter.to_str(bin_info.gc_interval()),
                ",".join(
                    ContigMultFormatter.to_str(contig_mult)
                    for contig_mult in bin_info.contigs_mults()
                ),
            ],
        )


class BinsReader:
    """PlasBin-flow bin info TSV reader."""

    @classmethod
    @contextmanager
    def open(cls, file: Path) -> Generator[BinsReader, None, None]:
        """Open TSV file for reading."""
        with file.open() as f_in:
            reader = BinsReader(csv.reader(f_in, delimiter="\t"))
            next(reader.__csv_reader)  # skip header
            yield reader

    def __init__(self, csv_reader: _csv._reader) -> None:
        """Initialize object."""
        self.__csv_reader = csv_reader

    def __iter__(
        self,
    ) -> Iterator[pbf_items.PBFBinInfo]:
        """Iterate sequence probability scores."""
        return (
            pbf_items.PBFBinInfo(
                row[0],
                float(row[1]),
                IntervalFormatter.from_str(row[2]),
                (
                    ContigMultFormatter.from_str(ctg_id_mult)
                    for ctg_id_mult in row[3].split(",")
                ),
            )
            for row in self.__csv_reader
        )
