"""GFA header API wrapper."""

from __future__ import annotations

from enum import StrEnum

import gfapy  # type: ignore[import-untyped]

from pangebin.gfa.tag import FieldType


class Tag(StrEnum):
    """GFA header tag types."""

    SKESA_FIX = "FX"
    STANDARDIZED = "Sd"
    PANASSEMBLY_TYPE = "PA"


class TagType(StrEnum):
    """GFA header tag types."""

    SKESA_FIX = FieldType.CHAR
    STANDARDIZED = FieldType.CHAR
    PANASSEMBLY_TYPE = FieldType.STRING

    @classmethod
    def from_tag(cls, tag: Tag) -> TagType:
        """Get field type from tag."""
        return cls(tag.name)


class SkesaFixTagValue(StrEnum):
    """Skesa fix header tag values."""

    YES = "Y"
    NO = "N"

    @classmethod
    def from_bool(
        cls,
        is_fixed: bool,  # noqa: FBT001
    ) -> SkesaFixTagValue:
        """Initialize from bool."""
        return cls.YES if is_fixed else cls.NO


class StandardizedTagValue(StrEnum):
    """Standardized header tag values."""

    YES = "Y"
    NO = "N"

    @classmethod
    def from_bool(
        cls,
        is_standardized: bool,  # noqa: FBT001
    ) -> StandardizedTagValue:
        """Initialize from bool."""
        return cls.YES if is_standardized else cls.NO


class PanAssemblyTypeTagValue(StrEnum):
    """PanAssemblyType header tag values."""

    BASE = "panassembly"
    WITH_CONTIGS = "panassembly_with_contigs"


def is_skesa_fixed(gfa: gfapy.Gfa) -> bool:
    """Check if a GFA graph is fixed."""
    return gfa.header.get(Tag.SKESA_FIX) == SkesaFixTagValue.YES


def set_skesa_fixed_header_tag(
    gfa: gfapy.Gfa,
    is_fixed: bool,  # noqa: FBT001
) -> None:
    """Set the Skesa fixed header tag in the GFA graph."""
    gfa.header.add(Tag.SKESA_FIX, SkesaFixTagValue.from_bool(is_fixed))


def is_standardized(gfa: gfapy.Gfa) -> bool:
    """Check if a GFA graph is standardized."""
    return gfa.header.get(Tag.STANDARDIZED) == StandardizedTagValue.YES


def set_standardized_header_tag(
    gfa: gfapy.Gfa,
    is_standardized: bool,  # noqa: FBT001
) -> None:
    """Set the standardized header tag in the GFA graph."""
    gfa.header.add(Tag.STANDARDIZED, StandardizedTagValue.from_bool(is_standardized))


def panassembly_type(gfa: gfapy.Gfa) -> PanAssemblyTypeTagValue:
    """Get PanAssemblyType header tag value."""
    return PanAssemblyTypeTagValue(gfa.header.get(Tag.PANASSEMBLY_TYPE))


def set_panassembly_type(gfa: gfapy.Gfa, value: PanAssemblyTypeTagValue) -> None:
    """Set PanAssemblyType header tag value."""
    gfa.header.add(Tag.PANASSEMBLY_TYPE, value)
