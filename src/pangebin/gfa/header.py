"""GFA header API wrapper."""

from __future__ import annotations

from enum import StrEnum

from pangebin.gfa.tag import FieldType


class Tag(StrEnum):
    """GFA header tag types."""

    SKESA_FIX = "FX"
    STANDARDIZED = "Sd"


class TagType(StrEnum):
    """GFA header tag types."""

    SKESA_FIX = FieldType.CHAR
    STANDARDIZED = FieldType.CHAR

    @classmethod
    def from_tag(cls, tag: Tag) -> TagType:
        """Get field type from tag."""
        return cls(tag.name)


class SkesaFixTagValue(StrEnum):
    """Skesa fix header tag values."""

    YES = "Y"
    NO = "N"


class StandardizedTagValue(StrEnum):
    """Standardized header tag values."""

    YES = "Y"
    NO = "N"
