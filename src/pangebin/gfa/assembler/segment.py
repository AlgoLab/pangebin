"""GFA segment API wrapper for assembler."""

from __future__ import annotations

from enum import StrEnum

import pangebin.assembly.items as asm_items


class NamePrefix(StrEnum):
    """Name prefixes for assemblers."""

    UNICYCLER = "uni"
    SKESA = "ske"

    @classmethod
    def from_assembler(cls, assembler_id: asm_items.Identifier) -> NamePrefix:
        """Get prefix from assembler."""
        match assembler_id:
            case asm_items.Identifier.UNICYCLER:
                return cls.UNICYCLER
            case asm_items.Identifier.SKESA:
                return cls.SKESA

    def to_assembler(self) -> asm_items.Identifier:
        """Get assembler from prefix."""
        match self:
            case NamePrefix.UNICYCLER:
                return asm_items.Identifier.UNICYCLER
            case NamePrefix.SKESA:
                return asm_items.Identifier.SKESA
