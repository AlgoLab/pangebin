"""GFA segment API wrapper."""

from __future__ import annotations

import logging
from enum import StrEnum
from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]

import pangebin.gfa.line as gfa_line
from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gfapy  # type: ignore[import-untyped]
    from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


_LOGGER = logging.getLogger(__name__)


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

    @classmethod
    def from_segment_line(
        cls,
        segment_line: GfaSegment,
        orientation: Orientation,
    ) -> OrientedFragment:
        """Get oriented fragment from segment line, by default the forward."""
        return cls(segment_line.name, orientation)

    @classmethod
    def from_segment_line_to_forward(cls, segment_line: GfaSegment) -> OrientedFragment:
        """Get forward fragment from segment line."""
        return cls.from_segment_line(segment_line, Orientation.FORWARD)

    @classmethod
    def from_segment_line_to_reverse(cls, segment_line: GfaSegment) -> OrientedFragment:
        """Get reverse fragment from segment line."""
        return cls.from_segment_line(segment_line, Orientation.REVERSE)

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
        """Get string representation, e.g. `frag+` or `frag-`."""
        return f"{self.__identifier}{self.__orientation}"


def get_segment_line_by_name(gfa: gfapy.Gfa, name: str) -> GfaSegment:
    """Get segment by name.

    Parameters
    ----------
    gfa : gfapy.Gfa
        GFA graph
    name : str
        Segment name

    Returns
    -------
    GfaSegment
        GFA segment line

    Raises
    ------
    ValueError
        Invalid segment name

    """
    # REFACTOR simplify using gfa.segment() method
    line: gfapy.Line | None = gfa.line(str(name))
    if line is None or line.record_type != gfa_line.Type.SEGMENT:
        _err_msg = f"Invalid segment name: {name}, line is {line}"
        _LOGGER.error(_err_msg)
        raise ValueError(_err_msg)
    return line


def gfa_oriented_fragments(
    gfa: gfapy.Gfa,
) -> Iterator[tuple[OrientedFragment, OrientedFragment]]:
    """Iterate over GFA forward and reverse fragments.

    Yields
    ------
    OrientedFragment
        Forward fragment
    OrientedFragment
        Reverse fragment
    """
    for segment_line in gfa.segments:
        yield (
            OrientedFragment.from_segment_line(segment_line, Orientation.FORWARD),
            OrientedFragment.from_segment_line(segment_line, Orientation.REVERSE),
        )


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


DEFAULT_ATTRIBUTE_STR_SEP = " "


def to_sequence_record(
    segment: GfaSegment,
    sep: str = DEFAULT_ATTRIBUTE_STR_SEP,
) -> SeqRecord:
    """Convert a GFA segment line to a sequence record.

    Parameters
    ----------
    segment: GfaSegment
        GFA segment line
    sep: str, optional
        string for separating GFA attributes, default is space

    Return
    ------
    SeqRecord
        Sequence record

    """
    return SeqRecord(
        Seq(segment.sequence),
        id=segment.name,
        name=segment.name,
        description=f"{__format_attributes_string(segment, sep=sep)}",
    )


def __format_attributes_string(
    segment: GfaSegment,
    sep: str = DEFAULT_ATTRIBUTE_STR_SEP,
) -> str:
    """Format GFA segment attributes as a string.

    Parameters
    ----------
    segment: GfaSegment
        GFA segment
    sep: str, optional
        string for separating GFA attributes, default is space

    Returns
    -------
    str
        GFA segment attributes in format sep.join(key:value)

    """
    return sep.join(
        [f"{tag_name}:{segment.get(tag_name)}" for tag_name in segment.tagnames],
    )
