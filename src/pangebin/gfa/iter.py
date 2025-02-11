"""GFA iterator module."""

from __future__ import annotations

from typing import TYPE_CHECKING

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pangebin.gfa.segment import Orientation, OrientedFragment

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gfapy  # type: ignore[import-untyped]
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


def __gfa_segment_attribute_to_string(segment: GfaSegment, sep: str = " ") -> str:
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


def sequence_records(
    graph: gfapy.Gfa,
    sep: str = " ",
) -> Iterator[SeqRecord]:
    """Iterate over sequence records corresponding to GFA segment lines.

    Parameters
    ----------
    graph: gfapy.Gfa
        GFA graph
    sep: str, optional
        string for separating GFA attributes, default is space

    Yields
    ------
    SeqRecord
        Sequence record

    """
    seg: GfaSegment
    for seg in graph.segments:
        yield SeqRecord(
            Seq(seg.sequence),
            id=seg.name,
            name=seg.name,
            description=f"{__gfa_segment_attribute_to_string(seg, sep=sep)}",
        )


def left_dovetails(
    graph: gfapy.Gfa,
    fragment: OrientedFragment,
) -> Iterator[OrientedFragment]:
    """Iterate over left dovetails.

    Parameters
    ----------
    graph: gfapy.Gfa
        GFA graph
    fragment: OrientedFragment
        fragment

    Yields
    ------
    OrientedFragment
        left dovetail

    """
    link_line: str
    for link_line in graph.segments[fragment.identifier()].dovetails_L:
        tab_link_line = link_line.split("\t")
        if tab_link_line[3] == fragment.identifier():
            yield OrientedFragment(tab_link_line[1], Orientation(tab_link_line[2]))
        yield OrientedFragment(
            tab_link_line[3],
            Orientation.from_reverse_string(tab_link_line[4]),
        )


def right_dovetails(
    graph: gfapy.Gfa,
    fragment: OrientedFragment,
) -> Iterator[OrientedFragment]:
    """Iterate over right dovetails.

    Parameters
    ----------
    graph: gfapy.Gfa
        GFA graph
    fragment: OrientedFragment
        fragment

    Yields
    ------
    OrientedFragment
        right dovetail

    """
    link_line: str
    for link_line in graph.segments[fragment.identifier()].dovetails_R:
        tab_link_line = link_line.split("\t")
        if tab_link_line[1] == fragment.identifier():
            yield OrientedFragment(tab_link_line[3], Orientation(tab_link_line[4]))
        yield OrientedFragment(
            tab_link_line[1],
            Orientation.from_reverse_string(tab_link_line[2]),
        )
