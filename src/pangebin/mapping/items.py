"""Blast mapping items."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


class BlastColumns(StrEnum):
    """Blast columns."""

    QSEQID = "qseqid"
    SSEQID = "sseqid"
    PIDENT = "pident"
    LENGTH = "length"
    MISMATCH = "mismatch"
    GAPOPEN = "gapopen"
    QSTART = "qstart"
    QEND = "qend"
    SSTART = "sstart"
    SEND = "send"
    EVALUE = "evalue"
    BITSCORE = "bitscore"


QSEQID_TYPE = str
SSEQID_TYPE = str
PIDENT_TYPE = float
LENGTH_TYPE = int
MISMATCH_TYPE = int
GAPOPEN_TYPE = int
QSTART_TYPE = int
QEND_TYPE = int
SSTART_TYPE = int
SEND_TYPE = int
EVALUE_TYPE = float
BITSCORE_TYPE = float

# REFACTOR use functions
BLAST6_COL_TYPES = {
    str(BlastColumns.QSEQID): QSEQID_TYPE,
    str(BlastColumns.SSEQID): SSEQID_TYPE,
    str(BlastColumns.PIDENT): PIDENT_TYPE,
    str(BlastColumns.LENGTH): LENGTH_TYPE,
    str(BlastColumns.MISMATCH): MISMATCH_TYPE,
    str(BlastColumns.GAPOPEN): GAPOPEN_TYPE,
    str(BlastColumns.QSTART): QSTART_TYPE,
    str(BlastColumns.QEND): QEND_TYPE,
    str(BlastColumns.SSTART): SSTART_TYPE,
    str(BlastColumns.SEND): SEND_TYPE,
    str(BlastColumns.EVALUE): EVALUE_TYPE,
    str(BlastColumns.BITSCORE): BITSCORE_TYPE,
}


class Mapping:
    """Mapping."""

    @classmethod
    def from_string(cls, line: str) -> Mapping:
        """Create mapping from string."""
        str_items: list[str] = line.split()
        return cls(
            QSEQID_TYPE(str_items[0]),
            SSEQID_TYPE(str_items[1]),
            PIDENT_TYPE(str_items[2]),
            LENGTH_TYPE(str_items[3]),
            MISMATCH_TYPE(str_items[4]),
            GAPOPEN_TYPE(str_items[5]),
            QSTART_TYPE(str_items[6]),
            QEND_TYPE(str_items[7]),
            SSTART_TYPE(str_items[8]),
            SEND_TYPE(str_items[9]),
            EVALUE_TYPE(str_items[10]),
            BITSCORE_TYPE(str_items[11]),
        )

    @classmethod
    def from_dataframe_row(cls, row: pd.Series) -> Mapping:
        """Create mapping from dataframe row."""
        return cls(
            QSEQID_TYPE(row[BlastColumns.QSEQID]),
            SSEQID_TYPE(row[BlastColumns.SSEQID]),
            PIDENT_TYPE(row[BlastColumns.PIDENT]),
            LENGTH_TYPE(row[BlastColumns.LENGTH]),
            MISMATCH_TYPE(row[BlastColumns.MISMATCH]),
            GAPOPEN_TYPE(row[BlastColumns.GAPOPEN]),
            QSTART_TYPE(row[BlastColumns.QSTART]),
            QEND_TYPE(row[BlastColumns.QEND]),
            SSTART_TYPE(row[BlastColumns.SSTART]),
            SEND_TYPE(row[BlastColumns.SEND]),
            EVALUE_TYPE(row[BlastColumns.EVALUE]),
            BITSCORE_TYPE(row[BlastColumns.BITSCORE]),
        )

    def __init__(  # noqa: PLR0913
        self,
        qseqid: str,
        sseqid: str,
        pident: float,
        length: int,
        mismatch: int,
        gapopen: int,
        qstart: int,
        qend: int,
        sstart: int,
        send: int,
        evalue: float,
        bitscore: float,
    ) -> None:
        """Initialize object."""
        self.__data = (
            qseqid,
            sseqid,
            pident,
            length,
            mismatch,
            gapopen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        )

    def qseqid(self) -> str:
        """Query ID."""
        return self.__data[0]

    def sseqid(self) -> str:
        """Subject ID."""
        return self.__data[1]

    def pident(self) -> float:
        """Percent identity."""
        return self.__data[2]

    def length(self) -> int:
        """Alignment length."""
        return self.__data[3]

    def mismatch(self) -> int:
        """Get number of mismatches."""
        return self.__data[4]

    def gapopen(self) -> int:
        """Get number of gap opens."""
        return self.__data[5]

    def qstart(self) -> int:
        """Query start."""
        return self.__data[6]

    def qend(self) -> int:
        """Query end."""
        return self.__data[7]

    def sstart(self) -> int:
        """Subject start."""
        return self.__data[8]

    def send(self) -> int:
        """Subject end."""
        return self.__data[9]

    def evalue(self) -> float:
        """Get e-value."""
        return self.__data[10]

    def bitscore(self) -> float:
        """Bit score."""
        return self.__data[11]

    def to_list(self) -> list:
        """Convert to list."""
        return list(self.__data)


class MappingInterval:
    """Mapping interval.

    Note
    ----
    The interval correspond to :math:`[start, end]`,
    with a length equals to :math:`end - start + 1`

    """

    @classmethod
    def from_string(cls, line: str) -> MappingInterval:
        """Create object from string."""
        str_items: list[str] = line.split(":")
        return cls(int(str_items[0]), int(str_items[1]))

    def __init__(self, start: int, end: int) -> None:
        """Initialize object."""
        self.__data = (start, end)

    def start(self) -> int:
        """Start."""
        return self.__data[0]

    def end(self) -> int:
        """End."""
        return self.__data[1]

    def set_start(self, start: int) -> None:
        """Set start."""
        self.__data = (start, self.__data[1])

    def set_end(self, end: int) -> None:
        """Set end."""
        self.__data = (self.__data[0], end)

    def __len__(self) -> int:
        """Get length."""
        return self.__data[1] - self.__data[0] + 1

    def __str__(self) -> str:
        """Get string representation."""
        return f"{self.__data[0]}:{self.__data[1]}"
