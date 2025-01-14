"""Segment API wrapper for GFA pangenome."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.assembler.segment as gfa_asm_segment
import pangebin.gfa.segment as gfa_segment
from pangebin import assembler
from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]

_COMMON_PREFIX = "pan"


class NamePrefix(StrEnum):
    """Segment prefix for GFA pangenome."""

    SUB_SKESA = f"{_COMMON_PREFIX}{gfa_asm_segment.NamePrefix.SKESA}"
    SUB_UNICYCLER = f"{_COMMON_PREFIX}{gfa_asm_segment.NamePrefix.UNICYCLER}"
    SUB_SKESA_AND_UNICYCLER = f"{_COMMON_PREFIX}both"

    @classmethod
    def from_segment_nature(cls, segment_nature: NatureTagValue) -> NamePrefix:
        """Get prefix from assembler."""
        match segment_nature:
            case NatureTagValue.SUB_SKESA_CONTIG:
                return cls.SUB_SKESA
            case NatureTagValue.SUB_UNICYCLER_CONTIG:
                return cls.SUB_UNICYCLER
            case NatureTagValue.SUB_SKESA_AND_UNICYCLER_CONTIGS:
                return cls.SUB_SKESA_AND_UNICYCLER


class Tag(StrEnum):
    """Segment tags.

    The tags from standardized assembly graphs are also defined for panassembly graphs.
    See `pangebin.gfa.assembler.segment.Tag`
    """

    OCCURENCE_IN_PATHS = "OC"
    CONTIG_LIST = "cl"
    CONTIG_PERCENTAGES = "cp"
    SEGMENT_NATURE = "sn"
    PANGENOME_PENALTY = "pp"


class TagType(StrEnum):
    """Segment tag types."""

    OCCURENCE_IN_PATH = FieldType.SIGNED_INT
    FROM_CONTIGS = FieldType.STRING
    CONTIG_PERCENTAGES = FieldType.INT_OR_FLOAT_ARRAY
    # FIXME this tag is renamed
    SEGMENT_NATURE = FieldType.CHAR
    PANGENOME_PENALTY = FieldType.FLOAT


# FIXME this tag is renamed
class NatureTagValue(StrEnum):
    """Segment nature tag values.

    The string representation is only one letter.
    """

    SUB_SKESA_CONTIG = "s"
    SUB_UNICYCLER_CONTIG = "u"
    SUB_SKESA_AND_UNICYCLER_CONTIGS = "b"

    @classmethod
    def from_name(cls, segment_name: str) -> NatureTagValue:
        """Get segment nature tag value from segment name."""
        match NamePrefix(gfa_segment.name_to_prefix(segment_name)):
            case NamePrefix.SUB_SKESA:
                return NatureTagValue.SUB_SKESA_CONTIG
            case NamePrefix.SUB_UNICYCLER:
                return NatureTagValue.SUB_UNICYCLER_CONTIG
            case NamePrefix.SUB_SKESA_AND_UNICYCLER:
                return NatureTagValue.SUB_SKESA_AND_UNICYCLER_CONTIGS

    @classmethod
    def from_contig_name(cls, contig_name: str) -> NatureTagValue:
        """Get segment nature tag value from contig name.

        Only return one of `SUB_SKESA_CONTIG` and `SUB_UNICYCLER_CONTIG`.
        """
        match gfa_asm_segment.NamePrefix(gfa_segment.name_to_prefix(contig_name)):
            case gfa_asm_segment.NamePrefix.SKESA:
                return NatureTagValue.SUB_SKESA_CONTIG
            case gfa_asm_segment.NamePrefix.UNICYCLER:
                return NatureTagValue.SUB_UNICYCLER_CONTIG

    @classmethod
    def from_assembler(cls, assembler_id: assembler.Identifier) -> NatureTagValue:
        """Get segment nature tag value from assembler.

        It returns only the segment nature for the two assemblers, not the pangenome.
        """
        match assembler_id:
            case assembler.Identifier.SKESA:
                return NatureTagValue.SUB_SKESA_CONTIG
            case assembler.Identifier.UNICYCLER:
                return NatureTagValue.SUB_UNICYCLER_CONTIG


def occurence_in_paths(segment: GfaSegment) -> int:
    """Get occurence in paths."""
    return segment.get(Tag.OCCURENCE_IN_PATHS)


def set_occurence_in_paths(segment: GfaSegment, occurence: int) -> None:
    """Set occurence in paths."""
    segment.set(Tag.OCCURENCE_IN_PATHS, occurence)


def contig_list(segment: GfaSegment) -> list[str]:
    """Get from contigs."""
    if segment.get(Tag.CONTIG_LIST) == ",":
        return []
    return segment.get(Tag.CONTIG_LIST).split(",")


def set_contig_list(segment: GfaSegment, contig_names: list[str]) -> None:
    """Set from contigs."""
    if not contig_names:
        segment.set(Tag.CONTIG_LIST, ",")
        return
    segment.set(Tag.CONTIG_LIST, ",".join(contig_names))


def append_to_contig_list(segment: GfaSegment, contig_name: str) -> None:
    """Append from contigs."""
    if segment.get(Tag.CONTIG_LIST) == ",":
        segment.set(Tag.CONTIG_LIST, contig_name)
        return

    segment.set(
        Tag.CONTIG_LIST,
        segment.get(Tag.CONTIG_LIST) + "," + contig_name,
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


def segment_nature(segment: GfaSegment) -> NatureTagValue:
    """Get segment nature tag value."""
    return NatureTagValue(segment.get(Tag.SEGMENT_NATURE))


def set_segment_nature(segment: GfaSegment, value: NatureTagValue) -> None:
    """Set segment nature tag value."""
    segment.set(Tag.SEGMENT_NATURE, value)


def pangenome_penalty(segment: GfaSegment) -> float:
    """Get pangenome penalty."""
    return segment.get(Tag.PANGENOME_PENALTY)


def set_pangenome_penalty(segment: GfaSegment, value: float) -> None:
    """Set pangenome penalty."""
    segment.set(Tag.PANGENOME_PENALTY, value)
