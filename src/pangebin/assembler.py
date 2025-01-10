"""Assembler module."""

from enum import StrEnum


class Item(StrEnum):
    """Assemblers."""

    UNICYCLER = "unicycler"
    SKESA = "skesa"
    PANGENOME = "pangenome"


class Prefix(StrEnum):
    """Contig prefixes."""

    UNICYCLER = "uni"
    SKESA = "ske"
