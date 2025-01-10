"""GFA segment API wrapper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]

from pangebin import assembler
from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


class Tag(StrEnum):
    """Segment tags."""

    LENGTH = "LN"
    KMER_COVERAGE = "KC"
    NORMALIZED_COVERAGE = "dp"
    OCCURENCE_IN_PATHS = "OC"
    FROM_CONTIGS = "cl"
    CONTIG_PERCENTAGES = "ll"
    FROM_ASSEMBLER = "aa"


class TagType(StrEnum):
    """Segment tag types."""

    LENGTH = FieldType.SIGNED_INT
    KMER_COVERAGE = FieldType.SIGNED_INT
    NORMALIZED_COVERAGE = FieldType.FLOAT
    OCCURENCE_IN_PATH = FieldType.SIGNED_INT
    FROM_CONTIGS = FieldType.STRING
    CONTIG_PERCENTAGES = FieldType.INT_OR_FLOAT_ARRAY
    FROM_ASSEMBLER = FieldType.CHAR


class FromAssemblerTagValue(StrEnum):
    """From assembler tag values.

    The string representation is only one letter.
    """

    SUB_UNICYCLER_AND_SKESA_CONTIGS = "b"
    SUB_UNICYCLER_CONTIG = "u"
    SUB_SKESA_CONTIG = "s"
    WHOLE_UNICYCLER_CONTIG = "U"
    WHOLE_SKESA_CONTIG = "S"

    @classmethod
    def from_panassembly(
        cls,
        panassembly_assembler: assembler.Item,
    ) -> FromAssemblerTagValue:
        """Get from panassembly assembler tag value."""
        # FIXME not really this, you can also have UU or SS
        match panassembly_assembler:
            case assembler.Item.PANGENOME:
                return cls.SUB_UNICYCLER_AND_SKESA_CONTIGS
            case assembler.Item.SKESA:
                return cls.SUB_SKESA_CONTIG
            case assembler.Item.UNICYCLER:
                return cls.SUB_UNICYCLER_CONTIG


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


def occurence_in_paths(segment: GfaSegment) -> int:
    """Get occurence in paths."""
    return segment.get(Tag.OCCURENCE_IN_PATHS)


def set_occurence_in_paths(segment: GfaSegment, occurence: int) -> None:
    """Set occurence in paths."""
    segment.set(Tag.OCCURENCE_IN_PATHS, occurence)


def from_contigs(segment: GfaSegment) -> list[str]:
    """Get from contigs."""
    return segment.get(Tag.FROM_CONTIGS).split(",")


def set_from_contigs(segment: GfaSegment, contig_names: list[str]) -> None:
    """Set from contigs."""
    segment.set(Tag.FROM_CONTIGS, ",".join(contig_names))


def append_from_contigs(segment: GfaSegment, contig_name: str) -> None:
    """Append from contigs."""
    segment.set(
        Tag.FROM_CONTIGS,
        segment.get(Tag.FROM_CONTIGS) + "," + contig_name,
    )


def contig_percentages(segment: GfaSegment) -> gfapy.NumericArray:
    """Get contig percentages."""
    return segment.get(Tag.CONTIG_PERCENTAGES)


def set_contig_percentages(
    segment: GfaSegment,
    contig_percentages: list[float],
) -> gfapy.NumericArray:
    """Set contig percentages."""
    numeric_array = gfapy.NumericArray(contig_percentages)
    segment.set(Tag.CONTIG_PERCENTAGES, numeric_array)
    return numeric_array


def append_contig_percentage(segment: GfaSegment, contig_percentage: float) -> None:
    """Append contig percentage."""
    contig_percentages(segment).append(contig_percentage)


def from_assembler(segment: GfaSegment) -> FromAssemblerTagValue:
    """Get from assembler tag value."""
    return FromAssemblerTagValue(segment.get(Tag.FROM_ASSEMBLER))


def set_from_assembler(segment: GfaSegment, value: FromAssemblerTagValue) -> None:
    """Set from assembler tag value."""
    segment.set(Tag.FROM_ASSEMBLER, value)


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
