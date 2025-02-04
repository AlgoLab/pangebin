"""Blast mapping items."""

from enum import StrEnum


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


# REFACTOR use functions
BLAST6_COL_TYPES = {
    str(BlastColumns.QSEQID): str,
    str(BlastColumns.SSEQID): str,
    str(BlastColumns.PIDENT): float,
    str(BlastColumns.LENGTH): int,
    str(BlastColumns.MISMATCH): int,
    str(BlastColumns.GAPOPEN): int,
    str(BlastColumns.QSTART): int,
    str(BlastColumns.QEND): int,
    str(BlastColumns.SSTART): int,
    str(BlastColumns.SEND): int,
    str(BlastColumns.EVALUE): float,
    str(BlastColumns.BITSCORE): float,
}
