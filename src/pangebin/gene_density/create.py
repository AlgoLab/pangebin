"""Gene density."""

from pathlib import Path

import gfapy  # type: ignore[import-untyped]
from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
    Path as GfaPath,
)

import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.segment as gfa_segment
import pangebin.mapping.items as map_items
import pangebin.mapping.ops as map_ops
import pangebin.std_asm_graph.fasta as std_asm_graph_fasta


def fragment_gene_density(
    panassembly_gfa_file: Path,
    gene_mapping_on_contigs_sam: Path,
) -> tuple[
    dict[str, float],
    dict[str, list[map_items.MappingInterval]],
]:
    """Return the ammping intervals for each fragment in a GFA graph.

    Parameters
    ----------
    panassembly_gfa_file : Path
        Path to the pan-assembly GFA file.
    gene_mapping_on_contigs_sam : Path
        Path to the Blast6 SAM mapping file.

    Returns
    -------
    dict[str, float]
        Fragment ID (str): Gene density
    dict[str, list[MappingInterval]]
        Fragment ID (str): List of MappingIntervals

    """
    panasm_graph = gfa_io.from_file(panassembly_gfa_file)
    all_contig_disjoint_intervals = {
        std_asm_graph_fasta.pansn_name_to_contig_name(
            ctg_name,
        ): map_ops.interval_union(intervals)
        for ctg_name, intervals in map_ops.subject_intervals_from_sam(
            gene_mapping_on_contigs_sam,
        ).items()
    }
    fragment_lengths: dict[str, int] = {
        segment_line.name: gfa_segment.length(segment_line)
        for segment_line in panasm_graph.segments
    }
    fragment_disjoint_mapping_intervals = {
        frag_name: map_ops.interval_union(intervals)
        for frag_name, intervals in fragment_mapping_intervals(
            panasm_graph,
            all_contig_disjoint_intervals,
            fragment_lengths,
        ).items()
    }
    return (
        {
            frag_name: gene_density(
                intervals,
                fragment_lengths[frag_name],
            )
            for frag_name, intervals in fragment_disjoint_mapping_intervals.items()
        },
        fragment_disjoint_mapping_intervals,
    )


def gene_density(
    disjoint_mapping_intervals: list[map_items.MappingInterval],
    sequence_length: int,
) -> float:
    """Return the gene density."""
    return (
        sum(len(interval) for interval in disjoint_mapping_intervals) / sequence_length
    )


def fragment_mapping_intervals(
    graph: gfapy.Gfa,
    all_contig_disjoint_intervals: dict[str, list[map_items.MappingInterval]],
    fragment_lengths: dict[str, int],
) -> dict[str, list[map_items.MappingInterval]]:
    """Return the mapping intervals for each fragment in a GFA graph.

    Parameters
    ----------
    graph : Gfa
        GFA graph
    all_contig_disjoint_intervals : list[MappingInterval]
        Sorted disjoint contig intervals
    fragment_lengths : dict[str, int]
        The length of each fragments

    Returns
    -------
    dict[str, list[MappingInterval]]
        Fragment ID (str): List of MappingIntervals

    """
    fragment_intervals: dict[str, list[map_items.MappingInterval]] = {
        segment_line.name: [] for segment_line in graph.segments
    }

    path_line: GfaPath
    for path_line in graph.paths:
        if path_line.name in all_contig_disjoint_intervals:
            dispatch_contig_intervals_to_fragment_intervals(
                path_line,
                all_contig_disjoint_intervals[path_line.name],
                fragment_lengths,
                fragment_intervals,
            )

    return fragment_intervals


def dispatch_contig_intervals_to_fragment_intervals(
    contig_path_line: GfaPath,
    contig_disjoint_intervals: list[map_items.MappingInterval],
    fragment_lengths: dict[str, int],
    fragment_intervals: dict[str, list[map_items.MappingInterval]],
) -> None:
    """Dispatch the contig sorted disjoint mapping interval to fragment intervals.

    Parameters
    ----------
    contig_path_line : GfaPath
        GFA path line corresponding to the contig
    contig_disjoint_intervals : list[MappingInterval]
        Sorted disjoint contig intervals
    fragment_lengths : dict[str, int]
        The length of each fragments
    fragment_intervals : dict[str, list[MappingInterval]]
        Mpping intervals for each fragment (unsorted, and not necessary disjoint)
    """
    if not contig_disjoint_intervals:
        return

    contig_fragment_names: list[str] = [
        or_id.name for or_id in contig_path_line.segment_names
    ]
    ctg_interval_idx = 0
    frag_in_ctg_idx = 0
    frag_start_pos_on_ctg = 1
    while ctg_interval_idx < len(contig_disjoint_intervals) and frag_in_ctg_idx < len(
        contig_fragment_names,
    ):
        frag_name = contig_fragment_names[frag_in_ctg_idx]
        frag_length = fragment_lengths[frag_name]
        ctg_interval = contig_disjoint_intervals[ctg_interval_idx]
        frag_end_pos_on_ctg = frag_start_pos_on_ctg + frag_length - 1
        if frag_end_pos_on_ctg < ctg_interval.start():
            frag_in_ctg_idx += 1
        elif frag_start_pos_on_ctg > ctg_interval.end():
            frag_start_pos_on_ctg += frag_length
            ctg_interval_idx += 1
        else:
            fragment_intervals[frag_name].append(
                map_items.MappingInterval(
                    max(ctg_interval.start() - frag_start_pos_on_ctg + 1, 1),
                    min(
                        frag_length - (frag_end_pos_on_ctg - ctg_interval.end()),
                        frag_length,
                    ),
                ),
            )
            if frag_end_pos_on_ctg == ctg_interval.end():
                frag_start_pos_on_ctg += frag_length
                frag_in_ctg_idx += 1
                ctg_interval_idx += 1
            elif frag_end_pos_on_ctg < ctg_interval.end():
                frag_start_pos_on_ctg += frag_length
                frag_in_ctg_idx += 1
            else:
                ctg_interval_idx += 1
