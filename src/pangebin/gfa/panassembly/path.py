"""Path GFA API wrapper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
        Path as GfaPath,
    )


class Tag(StrEnum):
    """Segment tags."""

    LENGTH = "LN"
    NORMALIZED_COVERAGE = "dp"


class TagType(StrEnum):
    """Segment tag types."""

    LENGTH = FieldType.SIGNED_INT
    NORMALIZED_COVERAGE = FieldType.FLOAT

    @classmethod
    def from_tag(cls, tag: Tag) -> TagType:
        """Get field type from tag."""
        return cls(tag.name)


def length(path: GfaPath) -> int:
    """Get path length."""
    return path.get(Tag.LENGTH)


def set_length(path: GfaPath, length: int) -> None:
    """Set path length."""
    path.set(Tag.LENGTH, length)


def normalized_coverage(path: GfaPath) -> float:
    """Get normalized coverage."""
    return path.get(Tag.NORMALIZED_COVERAGE)


def set_normalized_coverage(path: GfaPath, coverage: float) -> None:
    """Set normalized coverage."""
    path.set(Tag.NORMALIZED_COVERAGE, coverage)
