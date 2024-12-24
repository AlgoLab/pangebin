"""Assembler module."""

from enum import StrEnum


class Assembler(StrEnum):
    """Assemblers."""

    UNICYCLER = "unicycler"
    SKESA = "skesa"
    PANGENOME = "pangenome"


class ContigPrefix(StrEnum):
    """Contig prefixes."""

    UNICYCLER = "uni"
    SKESA = "ske"
    PANGENOME = "pan"
