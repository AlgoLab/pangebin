"""Link GFA API wrapper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

import pangebin.assembly.items as asm_items
import pangebin.gfa.panassembly.segment as gfa_pan_segment
from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]


class Tag(StrEnum):
    """Link tags."""

    LINK_ORIGIN = "lo"
    LINK_TYPE = "lt"


class TagType(StrEnum):
    """Link tag types."""

    LINK_ORIGIN = FieldType.CHAR
    LINK_TYPE = FieldType.STRING


class LinkOriginTagValue(StrEnum):
    """Link origin tag values."""

    SKESA_LINK = "s"
    UNICYCLER_LINK = "u"
    PANGENOME_LINK = "p"

    @classmethod
    def from_assembler(cls, assembler_id: asm_items.Identifier) -> LinkOriginTagValue:
        """Get link origin tag value from assembler.

        It only returns the link origin for the two assemblers, not the pangenome.
        """
        match assembler_id:
            case asm_items.Identifier.SKESA:
                return cls.SKESA_LINK
            case asm_items.Identifier.UNICYCLER:
                return cls.UNICYCLER_LINK


def link_origin(link: GfaLink) -> LinkOriginTagValue:
    """Get link origin tag value."""
    return LinkOriginTagValue(link.get(Tag.LINK_ORIGIN))


def set_link_origin(link: GfaLink, value: LinkOriginTagValue) -> None:
    """Set link origin tag value."""
    link.set(Tag.LINK_ORIGIN, value)


def link_type(
    link: GfaLink,
) -> tuple[
    gfa_pan_segment.NatureTagValue,
    gfa_pan_segment.NatureTagValue,
]:
    """Get link type tag value."""
    return (
        gfa_pan_segment.NatureTagValue(link.get(Tag.LINK_TYPE)[0]),
        gfa_pan_segment.NatureTagValue(link.get(Tag.LINK_TYPE)[1]),
    )


def set_link_type(
    link: GfaLink,
    pred_segment_nature: gfa_pan_segment.NatureTagValue,
    succ_segment_nature: gfa_pan_segment.NatureTagValue,
) -> str:
    """Set from assembler tag value."""
    value = f"{pred_segment_nature.value}{succ_segment_nature.value}"
    link.set(Tag.LINK_TYPE, value)
    return value
