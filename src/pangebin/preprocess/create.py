"""Create preprocessed files."""

import contextlib
from itertools import product

import gfapy  # type: ignore[import-untyped]

from pangebin.graph_utils import etos, extract_node


def remove_nodes(
    gfa: gfapy.Gfa,
    min_contig_length: int,
) -> None:
    """Remove nodes that are smaller than min_contig_length.

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
            pass
        else:
            continue
        right_edges = list(seg.dovetails_R)
        left_edges = list(seg.dovetails_L)
        right_nodes = set()  # edge(node, orient, 'r')
        left_nodes = set()

        for e in right_edges:
            n = extract_node(e, "r", seg.name)
            if n is not None:
                right_nodes.add(n)
                gfa.rm(e)
        for e in left_edges:
            n = extract_node(e, "l", seg.name)
            if n is not None:
                left_nodes.add(n)
                gfa.rm(e)

        gfa.rm(seg)
        gfa.validate()
        for edge in product(left_nodes, right_nodes):
            new_edge = etos(edge)
            with contextlib.suppress(gfapy.NotUniqueError):
                gfa.add_line(new_edge)
        gfa.validate()
