"""Create standardized files."""

import contextlib
import logging
from itertools import product

import gfapy  # type: ignore[import-untyped]

import pangebin.assembly.items as asm_items
import pangebin.gfa.ops as gfa_ops
import pangebin.gfa.segment as gfa_segment
from pangebin.gfa.assembler import segment as gfa_asm_segment
from pangebin.gfa.link import Link
from pangebin.gfa.segment import OrientedFragment
from pangebin.std_asm_graph.config import Config

_LOGGER = logging.getLogger(__name__)


def standardize_assembly_graph(
    graph: gfapy.Gfa,
    assembler_id: asm_items.Identifier,
    config: Config,
) -> None:
    """Standardize assembly graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph
    assembler_id : asm_items.Identifier
        Assembler
    config : Config
        Configuration

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Standardizing assembly graph: %s", assembler_id)

    gfa_ops.set_segment_length_tags(graph)

    gfa_ops.rename_segments(
        graph,
        gfa_asm_segment.NamePrefix.from_assembler(assembler_id),
    )

    if assembler_id == asm_items.Identifier.SKESA:
        gfa_ops.convert_kmer_coverage_to_normalized_coverage(graph)

    transform_small_contigs_into_links(graph, config.min_contig_length())

    gfa_ops.set_standardized_header_tag(graph)


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
                    # XXX the exception should not happen
                    gfa.add_line(link.to_gfa_link_line())
            gfa.validate()
