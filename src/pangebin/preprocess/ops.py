"""Create preprocessed files."""

import contextlib
from itertools import product

import gfapy  # type: ignore[import-untyped]

from pangebin.gfa.items import Link, OrientedFragment


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
    for seg in gfa.segments:
        if (seg.LN is not None and min_contig_length >= seg.LN) or (
            len(seg.sequence) <= min_contig_length
        ):
            left_edge_lines = list(seg.dovetails_L)
            right_edge_lines = list(seg.dovetails_R)
            predecessors: list[OrientedFragment] = []
            successors: list[OrientedFragment] = []

            for left_link_line in left_edge_lines:
                pred = OrientedFragment.from_left_dovetail_line(
                    left_link_line,
                    seg.name,
                )
                if pred.identifier() != seg.name:
                    predecessors.append(pred)
                    gfa.rm(left_link_line)

            for right_link_line in right_edge_lines:
                succ = OrientedFragment.from_right_dovetail_line(
                    right_link_line,
                    seg.name,
                )
                if succ.identifier() != seg.name:
                    successors.append(succ)
                    gfa.rm(right_link_line)

            gfa.rm(seg)
            gfa.validate()
            for link in (
                Link(pred, succ) for pred, succ in product(predecessors, successors)
            ):
                with contextlib.suppress(gfapy.NotUniqueError):
                    gfa.add_line(link.to_gfa_link_line())
            gfa.validate()
