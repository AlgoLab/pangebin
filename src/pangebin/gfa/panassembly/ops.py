"""GFA pan-assembly operations."""

from collections.abc import Iterator

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.panassembly.path as gfa_pan_path
import pangebin.gfa.panassembly.segment as gfa_pan_segment
import pangebin.gfa.path as gfa_path


def fragment_max_contig_coverages(
    panasm_graph: gfapy.Gfa,
) -> Iterator[tuple[str, float]]:
    """Iterate over the maximum of the contig coverages for each fragment."""
    return (
        (
            segment_line.name,
            max(
                gfa_pan_path.normalized_coverage(
                    gfa_path.get_path_line_by_name(panasm_graph, ctg_id),
                )
                for ctg_id in gfa_pan_segment.contig_list(segment_line)
            ),
        )
        for segment_line in panasm_graph.segments
    )
