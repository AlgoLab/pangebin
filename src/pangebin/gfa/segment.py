"""GFA segment API wrapper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


class Orientation(StrEnum):
    """Orientation of a fragment."""

    FORWARD = "+"
    REVERSE = "-"

    @classmethod
    def from_reverse_string(cls, reverse_string: str) -> Orientation:
        """Get orientation from reverse string."""
        match reverse_string:
            case "+":
                return cls.REVERSE
            case "-":
                return cls.FORWARD
            case _:
                _err_msg = f"Invalid reverse string: {reverse_string}"
                raise ValueError(_err_msg)

    def reverse(self) -> Orientation:
        """Get reverse orientation."""
        match self:
            case Orientation.FORWARD:
                return Orientation.REVERSE
            case Orientation.REVERSE:
                return Orientation.FORWARD


class OrientedFragment:
    """Oriented fragment."""

    @classmethod
    def from_left_dovetail_line(
        cls,
        link_line: GfaLink,
        segment_name: str,
    ) -> OrientedFragment:
        """Get fragment from left dovetail line."""
        if link_line.to_segment.name == segment_name:
            return cls(link_line.from_segment.name, Orientation(link_line.from_orient))
        return cls(
            link_line.to_segment.name,
            Orientation.from_reverse_string(link_line.to_orient),
        )

    @classmethod
    def from_right_dovetail_line(
        cls,
        link_line: GfaLink,
        segment_name: str,
    ) -> OrientedFragment:
        """Get fragment from right dovetail line."""
        if link_line.from_segment.name == segment_name:
            return cls(link_line.to_segment.name, Orientation(link_line.to_orient))
        return cls(
            link_line.from_segment.name,
            Orientation.from_reverse_string(link_line.from_orient),
        )

    def __init__(self, identifier: str, orientation: Orientation) -> None:
        """Initialize object."""
        self.__identifier = identifier
        self.__orientation = orientation

    def identifier(self) -> str:
        """Get identifier."""
        return self.__identifier

    def orientation(self) -> Orientation:
        """Get orientation."""
        return self.__orientation

    def to_reverse(self) -> OrientedFragment:
        """Get reverse fragment."""
        return OrientedFragment(self.__identifier, self.__orientation.reverse())

    def is_forward(self) -> bool:
        """Check if fragment is forward."""
        return self.__orientation == Orientation.FORWARD

    def is_reverse(self) -> bool:
        """Check if fragment is reverse."""
        return self.__orientation == Orientation.REVERSE

    def __str__(self) -> str:
        """Get string representation."""
        return f"{self.__identifier}\t{self.__orientation}"


class Tag(StrEnum):
    """Segment tags."""

    LENGTH = "LN"
    KMER_COVERAGE = "KC"
    NORMALIZED_COVERAGE = "dp"


class TagType(StrEnum):
    """Segment tag types."""

    LENGTH = FieldType.SIGNED_INT
    KMER_COVERAGE = FieldType.SIGNED_INT
    NORMALIZED_COVERAGE = FieldType.FLOAT

    @classmethod
    def from_tag(cls, tag: Tag) -> TagType:
        """Get field type from tag."""
        return cls(tag.name)


def length(segment: GfaSegment) -> int:
    """Get segment length."""
    return segment.get(Tag.LENGTH)


def set_length(segment: GfaSegment, length: int | None = None) -> None:
    """Set segment length."""
    if length is None:
        length = len(segment.sequence)
    segment.set(Tag.LENGTH, length)


def kmer_coverage(segment: GfaSegment) -> int:
    """Get kmer coverage."""
    return segment.get(Tag.KMER_COVERAGE)


def set_kmer_coverage(segment: GfaSegment, coverage: int) -> None:
    """Set kmer coverage."""
    segment.set(Tag.KMER_COVERAGE, coverage)


def normalized_coverage(segment: GfaSegment) -> float:
    """Get normalized coverage."""
    return segment.get(Tag.NORMALIZED_COVERAGE)


def set_normalized_coverage(segment: GfaSegment, coverage: float) -> None:
    """Set normalized coverage."""
    segment.set(Tag.NORMALIZED_COVERAGE, coverage)


__SEP = "_"


def format_name(segment_prefix: str, integer: int) -> str:
    """Format segment name."""
    return f"{segment_prefix}{__SEP}{integer}"


def name_to_prefix(segment_name: str) -> str:
    """Get assembler prefix from segment name."""
    return segment_name.split(__SEP)[0]


def name_to_integer(segment_name: str) -> int:
    """Get integer from segment name."""
    return int(segment_name.split(__SEP)[1])
