"""Link GFA API wrapper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]

import pangebin.gfa.segment as gfa_segment
from pangebin import assembler
from pangebin.gfa import line as gfa_line
from pangebin.gfa.segment import Orientation, OrientedFragment
from pangebin.gfa.tag import FieldType

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gfapy  # type: ignore[import-untyped]


class Link:
    """Link."""

    @classmethod
    def from_link_line(cls, link_line: GfaLink) -> Link:
        """Get link from link line."""
        return cls(
            OrientedFragment(link_line.from_segment.name, link_line.from_orient),
            OrientedFragment(link_line.to_segment.name, link_line.to_orient),
            link_line.overlap,
        )

    def __init__(
        self,
        predecessor: OrientedFragment,
        successor: OrientedFragment,
        overlap_match: int = 0,
    ) -> None:
        """Initialize object."""
        self.__predecessor = predecessor
        self.__successor = successor
        self.__overlap_match = overlap_match

    def predecessor(self) -> OrientedFragment:
        """Get predecessor."""
        return self.__predecessor

    def successor(self) -> OrientedFragment:
        """Get successor."""
        return self.__successor

    def overlap_match(self) -> int:
        """Get match."""
        return self.__overlap_match

    def to_reverse(self) -> Link:
        """Get reverse link."""
        return Link(
            self.__successor.to_reverse(),
            self.__predecessor.to_reverse(),
            self.__overlap_match,
        )

    def to_gfa_link_line(self) -> GfaLink:
        """Get GFA link line."""
        return GfaLink(
            f"{gfa_line.Type.LINK}\t{self.__predecessor}\t{self.__successor}\t{self.__overlap_match}M",
        )

    def simplify(self) -> None:
        """Humanly simplify the link.

        Just in the case the two orientations are reversed.
        """
        if (
            self.__predecessor.orientation() == Orientation.REVERSE
            and self.__successor.orientation() == Orientation.REVERSE
        ):
            self.__predecessor, self.__successor = (
                self.__successor.to_reverse(),
                self.__predecessor.to_reverse(),
            )


class Tag(StrEnum):
    """Link tags."""

    FROM_ASSEMBLER = "aa"
    FROM_ASSEMBLER_TO_ASSEMBLER = "lt"


class TagType(StrEnum):
    """Link tag types."""

    FROM_ASSEMBLER = FieldType.CHAR
    FROM_ASSEMBLER_TO_ASSEMBLER = FieldType.STRING


class FromAssemblerTagValue(StrEnum):
    """From assembler tag values."""

    UNICYCLER_LINK = "u"
    SKESA_LINK = "s"
    PANGENOME_LINK = "p"

    @classmethod
    def from_assembler(cls, assembler_value: assembler.Item) -> FromAssemblerTagValue:
        """Get from assembler tag value."""
        match assembler_value:
            case assembler.Item.UNICYCLER:
                return FromAssemblerTagValue.UNICYCLER_LINK
            case assembler.Item.SKESA:
                return FromAssemblerTagValue.SKESA_LINK
            case assembler.Item.PANGENOME:
                return FromAssemblerTagValue.PANGENOME_LINK


def from_assembler(link: GfaLink) -> FromAssemblerTagValue:
    """Get from assembler tag value."""
    return FromAssemblerTagValue(link.get(Tag.FROM_ASSEMBLER))


def set_from_assembler(link: GfaLink, value: FromAssemblerTagValue) -> None:
    """Set from assembler tag value."""
    link.set(Tag.FROM_ASSEMBLER, value)


def from_assembler_to_assembler(
    link: GfaLink,
) -> tuple[gfa_segment.FromAssemblerTagValue, gfa_segment.FromAssemblerTagValue]:
    """Get from assembler tag value."""
    return (
        gfa_segment.FromAssemblerTagValue(link.get(Tag.FROM_ASSEMBLER_TO_ASSEMBLER)[0]),
        gfa_segment.FromAssemblerTagValue(link.get(Tag.FROM_ASSEMBLER_TO_ASSEMBLER)[1]),
    )


def set_from_assembler_to_assembler(
    link: GfaLink,
    from_assembler: gfa_segment.FromAssemblerTagValue,
    to_assembler: gfa_segment.FromAssemblerTagValue,
) -> str:
    """Set from assembler tag value."""
    value = f"{from_assembler.value}{to_assembler.value}"
    link.set(Tag.FROM_ASSEMBLER_TO_ASSEMBLER, value)
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
            gfa.segments[pred.identifier()].dovetails_R,
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
            gfa.segments[pred.identifier()].dovetails_L,
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


def link_or_its_reversed_exists(
    gfa: gfapy.Gfa,
    link: Link,
) -> tuple[GfaLink, Orientation] | None:
    """Check if link or its reversed exists.

    Parameters
    ----------
    gfa : gfapy.Gfa
        GFA graph
    link : Link
        Link

    Returns
    -------
    GfaLink, Orientation
        The corresponding gfa link line
        and if it corresponds to the fowward or the reverse of the given link
    None
        If the link or its reversed does not exist in the graph

    """
    pred = link.predecessor()
    if pred.is_forward():
        succ = link.successor()
        r_link_line_iter: Iterator[GfaLink] = iter(
            gfa.segments[pred.identifier()].dovetails_R,
        )
        link_line = next(r_link_line_iter, None)
        while link_line is not None:
            if (
                OrientedFragment.from_right_dovetail_line(link_line, pred.identifier())
                == succ
            ):
                return link_line, Orientation.FORWARD
            link_line = next(r_link_line_iter, None)
    else:
        succ_rev = link.successor().to_reverse()
        l_link_line_iter: Iterator[GfaLink] = iter(
            gfa.segments[pred.identifier()].dovetails_L,
        )
        link_line = next(l_link_line_iter, None)
        while link_line is not None:
            if (
                OrientedFragment.from_left_dovetail_line(link_line, pred.identifier())
                == succ_rev
            ):
                return link_line, Orientation.REVERSE
            link_line = next(l_link_line_iter, None)
    return None
