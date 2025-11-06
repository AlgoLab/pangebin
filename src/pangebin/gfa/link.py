"""Link GFA API wrapper."""

from __future__ import annotations

from typing import TYPE_CHECKING

from gfapy.line.edge import Link as GfaLink  # type: ignore[import-untyped]

from pangebin.gfa import line as gfa_line
from pangebin.gfa.segment import Orientation, OrientedFragment, get_segment_line_by_name

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gfapy  # type: ignore[import-untyped]


class Link:
    """Link."""

    @classmethod
    def from_link_line(cls, link_line: GfaLink) -> Link:
        """Get link from link line."""
        return cls(
            OrientedFragment(
                link_line.field_to_s("from_segment"),
                Orientation(link_line.from_orient),
            ),
            OrientedFragment(
                link_line.field_to_s("to_segment"),
                Orientation(link_line.to_orient),
            ),
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

    def to_gfa_line(self) -> GfaLink:
        """Get GFA link line."""
        return GfaLink(
            (
                f"{gfa_line.Type.LINK}"
                f"\t{self.__predecessor.identifier()}\t{self.__predecessor.orientation()}"
                f"\t{self.__successor.identifier()}\t{self.__successor.orientation()}"
                f"\t{self.__overlap_match}M"
            ),
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

    def __str__(self) -> str:
        """Get string representation, e.g. `u+_v-`."""
        return f"{self.__predecessor}_{self.__successor}"


def link_line_or_its_reversed_from_link(
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
            get_segment_line_by_name(gfa, pred.identifier()).dovetails_R,
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
            get_segment_line_by_name(gfa, pred.identifier()).dovetails_L,
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


def gfa_links(gfa: gfapy.Gfa) -> Iterator[Link]:
    """Iterate over only one version of each link.

    Yields
    ------
    Link
        One version of the link
    """
    for link_line in gfa.dovetails:
        yield Link.from_link_line(link_line)


def predecessors(
    gfa: gfapy.Gfa,
    oriented_fragment: OrientedFragment,
) -> Iterator[OrientedFragment]:
    """Iterate over all predecessors.

    Yields
    ------
    OrientedFragment
        Oriented predecessor
    """
    if oriented_fragment.is_forward():
        return (
            OrientedFragment.from_left_dovetail_line(
                link_line,
                oriented_fragment.identifier(),
            )
            for link_line in get_segment_line_by_name(
                gfa,
                oriented_fragment.identifier(),
            ).dovetails_L
        )
    return (
        OrientedFragment.from_right_dovetail_line(
            link_line,
            oriented_fragment.identifier(),
        ).to_reverse()
        for link_line in get_segment_line_by_name(
            gfa,
            oriented_fragment.identifier(),
        ).dovetails_R
    )


def incoming_links(
    gfa: gfapy.Gfa,
    oriented_fragment: OrientedFragment,
) -> Iterator[Link]:
    """Get incoming links.

    Iterate over only one version of the link.

    Yields
    ------
    Link
        Link where the successor is the given oriented fragment
    """
    return (
        Link(predecessor, oriented_fragment)
        for predecessor in predecessors(gfa, oriented_fragment)
    )


def successors(
    gfa: gfapy.Gfa,
    oriented_fragment: OrientedFragment,
) -> Iterator[OrientedFragment]:
    """Iterate over all successors.

    Yields
    ------
    OrientedFragment
        Oriented successor
    """
    if oriented_fragment.is_forward():
        return (
            OrientedFragment.from_right_dovetail_line(
                link_line,
                oriented_fragment.identifier(),
            )
            for link_line in get_segment_line_by_name(
                gfa,
                oriented_fragment.identifier(),
            ).dovetails_R
        )
    return (
        OrientedFragment.from_left_dovetail_line(
            link_line,
            oriented_fragment.identifier(),
        ).to_reverse()
        for link_line in get_segment_line_by_name(
            gfa,
            oriented_fragment.identifier(),
        ).dovetails_L
    )


def outgoing_links(
    gfa: gfapy.Gfa,
    oriented_fragment: OrientedFragment,
) -> Iterator[Link]:
    """Get outgoing links.

    Iterate over only one version of the link.

    Yields
    ------
    Link
        Link where the predecessor is the given oriented fragment
    """
    return (
        Link(oriented_fragment, successor)
        for successor in successors(gfa, oriented_fragment)
    )
