"""Create preprocessed files."""

import contextlib
from itertools import product

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.segment as gfa_segment
from pangebin.gfa.link import Link
from pangebin.gfa.segment import OrientedFragment


def transform_small_contigs_into_links(
    gfa: gfapy.Gfa,
    min_contig_length: int,
) -> None:
    """Remove small contigs and transform them into links.

    Parameters
    ----------
    gfa : gfapy.Gfa
        GFA graph
    min_contig_length : int
        Minimum contig length

    Warnings
    --------
    This function mutates the GFA graph

    """
    for segment in gfa.segments:
        if gfa_segment.length(segment) <= min_contig_length:
            left_edge_lines = list(segment.dovetails_L)
            right_edge_lines = list(segment.dovetails_R)
            predecessors: list[OrientedFragment] = []
            successors: list[OrientedFragment] = []

            for left_link_line in left_edge_lines:
                pred = OrientedFragment.from_left_dovetail_line(
                    left_link_line,
                    segment.name,
                )
                if pred.identifier() != segment.name:
                    predecessors.append(pred)
                    gfa.rm(left_link_line)

            for right_link_line in right_edge_lines:
                succ = OrientedFragment.from_right_dovetail_line(
                    right_link_line,
                    segment.name,
                )
                if succ.identifier() != segment.name:
                    successors.append(succ)
                    gfa.rm(right_link_line)

            gfa.rm(segment)
            gfa.validate()
            for link in (
                Link(pred, succ) for pred, succ in product(predecessors, successors)
            ):
                with contextlib.suppress(gfapy.NotUniqueError):
                    gfa.add_line(link.to_gfa_link_line())
            gfa.validate()
