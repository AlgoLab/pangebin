"""GFA items module."""

from __future__ import annotations

from enum import StrEnum

import gfapy  # type: ignore[import-untyped]


class GFALineType(StrEnum):
    """GFA line types."""

    COMMENT = "#"
    HEADER = "H"
    SEGMENT = "S"
    LINK = "L"
    JUMP = "J"
    CONTAINMENT = "C"
    PATH = "P"
    WALK = "W"


class GFAFieldType(StrEnum):
    """GFA field types."""

    CHAR = "A"
    SIGNED_INT = "i"
    FLOAT = "f"
    STRING = "Z"
    JSON = "J"
    BYTE_ARRAY = "H"
    INT_OR_FLOAT_ARRAY = "B"


SKESA_FIX_HEADER_TAG = "FX"
SKESA_FIX_HEADER_TAG_TYPE = GFAFieldType.CHAR


class SkesaFixHeaderTagValue(StrEnum):
    """Skesa fix header tag values."""

    YES = "Y"
    NO = "N"


class SequenceTag(StrEnum):
    """Sequence tags."""

    LENGTH = "LN"
    NORMALIZED_COVERAGE = "dp"


class SequenceTagType(StrEnum):
    """Sequence tag types."""

    LENGTH = GFAFieldType.SIGNED_INT
    NORMALIZED_COVERAGE = GFAFieldType.FLOAT


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
        link_line: gfapy.Line,
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
        link_line: gfapy.Line,
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

    def __str__(self) -> str:
        """Get string representation."""
        return f"{self.__identifier}\t{self.__orientation}"


class Link:
    """Link."""

    def __init__(
        self,
        predecessor: OrientedFragment,
        successor: OrientedFragment,
    ) -> None:
        """Initialize object."""
        self.__predecessor = predecessor
        self.__successor = successor

    def predecessor(self) -> OrientedFragment:
        """Get predecessor."""
        return self.__predecessor

    def successor(self) -> OrientedFragment:
        """Get successor."""
        return self.__successor

    def to_reverse(self) -> Link:
        """Get reverse link."""
        return Link(self.__successor.to_reverse(), self.__predecessor.to_reverse())

    def to_gfa_link_line(self) -> gfapy.Line:
        """Get GFA link line."""
        return gfapy.Line(
            f"L\t{self.__predecessor}\t{self.__successor}\t0M",
        )
