"""GFA segment API wrapper for assembler."""

from __future__ import annotations

from enum import StrEnum

from pangebin import assembler


class NamePrefix(StrEnum):
    """Name prefixes for assemblers."""

    UNICYCLER = "uni"
    SKESA = "ske"

    @classmethod
    def from_assembler(cls, assembler_id: assembler.Identifier) -> NamePrefix:
        """Get prefix from assembler."""
        match assembler_id:
            case assembler.Identifier.UNICYCLER:
                return cls.UNICYCLER
            case assembler.Identifier.SKESA:
                return cls.SKESA
