"""GFA iterator module."""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING

import pangebin.gfa.line as gfa_line
import pangebin.input_output as io

if TYPE_CHECKING:
    from pathlib import Path

    from Bio.SeqRecord import SeqRecord


import gfapy  # type: ignore[import-untyped]
from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
    Path as GfaPath,
)
from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]

import pangebin.gfa.segment as gfa_segment

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    import gfapy  # type: ignore[import-untyped]


def segment_lines(gfa_file: Path) -> Iterator[GfaSegment]:
    """Get a segment line iterator."""
    with io.open_file_read(gfa_file) as f_in:
        yield from (
            GfaSegment(line) for line in f_in if line[0] == gfa_line.Type.SEGMENT
        )


def sequence_records(
    gfa_file: Path,
    sep: str = " ",
) -> Iterator[SeqRecord]:
    """Iterate over sequence records corresponding to GFA segment lines.

    Parameters
    ----------
    gfa_file: Path
        GFA file
    sep: str, optional
        string for separating GFA attributes, default is space

    Yields
    ------
    SeqRecord
        Sequence record

    """
    for segment_line in segment_lines(gfa_file):
        yield gfa_segment.to_sequence_record(segment_line, sep=sep)


def contig_path(gfa_file: Path) -> Iterator[GfaPath]:
    """Get contig paths from a GFA graph."""
    with io.open_file_read(gfa_file) as f_in:
        yield from (GfaPath(line) for line in f_in if line[0] == gfa_line.Type.PATH)


def left_dovetails(
    graph: gfapy.Gfa,
    fragment: gfa_segment.OrientedFragment,
) -> Iterator[gfa_segment.OrientedFragment]:
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
        tab_link_line = link_line.split()
        if tab_link_line[3] == fragment.identifier():
            yield gfa_segment.OrientedFragment(
                tab_link_line[1],
                gfa_segment.Orientation(tab_link_line[2]),
            )
        yield gfa_segment.OrientedFragment(
            tab_link_line[3],
            gfa_segment.Orientation.from_reverse_string(tab_link_line[4]),
        )


def right_dovetails(
    graph: gfapy.Gfa,
    fragment: gfa_segment.OrientedFragment,
) -> Iterator[gfa_segment.OrientedFragment]:
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
        tab_link_line = link_line.split()
        if tab_link_line[1] == fragment.identifier():
            yield gfa_segment.OrientedFragment(
                tab_link_line[3],
                gfa_segment.Orientation(tab_link_line[4]),
            )
        yield gfa_segment.OrientedFragment(
            tab_link_line[1],
            gfa_segment.Orientation.from_reverse_string(tab_link_line[2]),
        )
