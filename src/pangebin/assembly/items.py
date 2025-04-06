"""Assembly items."""

from __future__ import annotations

from enum import IntEnum, StrEnum


class Identifier(StrEnum):
    """Assemblers."""

    SKESA = "SKESA"
    UNICYCLER = "Unicycler"


class HaplotypeID(IntEnum):
    """Haplotype IDs."""

    SKESA = 1
    UNICYCLER = 2
