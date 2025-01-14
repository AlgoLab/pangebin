"""Link GFA API wrapper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

import pangebin.gfa.panassembly.segment as gfa_pan_segment
from pangebin import assembler
from pangebin.gfa.segment import OrientedFragment, get_segment_line_by_name
from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gfapy  # type: ignore[import-untyped]
    from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]

    from pangebin.gfa.link import Link


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
    def from_assembler(cls, assembler_id: assembler.Identifier) -> LinkOriginTagValue:
        """Get link origin tag value from assembler.

        It only returns the link origin for the two assemblers, not the pangenome.
        """
        match assembler_id:
            case assembler.Identifier.SKESA:
                return cls.SKESA_LINK
            case assembler.Identifier.UNICYCLER:
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


def get_link_line_by_link_definition(gfa: gfapy.Gfa, link: Link) -> GfaLink | None:
    """Get link line by link definition.

    Parameters
    ----------
    gfa : gfapy.Gfa
        GFA graph
    link : Link
        Link

    Returns
    -------
    GfaLink | None
        GFA link line, None if not found

    """
    pred = link.predecessor()
    if pred.is_forward():
        succ = link.successor()
        r_link_line_iter: Iterator[GfaLink] = iter(
            get_segment_line_by_name(gfa, pred.identifier()).dovetails_R,
        )
        link_line = next(r_link_line_iter, None)
        while link_line is not None:
            if (
                OrientedFragment.from_right_dovetail_line(link_line, pred.identifier())
                == succ
            ):
                return link_line
            link_line = next(r_link_line_iter, None)
    else:
        succ_rev = link.successor().to_reverse()
        l_link_line_iter: Iterator[GfaLink] = iter(
            get_segment_line_by_name(gfa, pred.identifier()).dovetails_L,
        )
        link_line = next(l_link_line_iter, None)
        while link_line is not None:
            if (
                OrientedFragment.from_left_dovetail_line(link_line, pred.identifier())
                == succ_rev
            ):
                return link_line
            link_line = next(l_link_line_iter, None)
    return None
