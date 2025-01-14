"""Assembler module."""

from __future__ import annotations

from enum import StrEnum


class Identifier(StrEnum):
    """Assemblers."""

    UNICYCLER = "Unicycler"
    SKESA = "Skesa"

    # FIXME check the removal of PANGENOME item
