"""Standardize assembly GFA graph operations."""

from collections.abc import Iterator

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.segment as gfa_segment


def contig_coverages(asm_graph: gfapy.Gfa) -> Iterator[tuple[str, float]]:
    """Iterate over the contigs and their coverages."""
    return (
        (segment_line.name, gfa_segment.normalized_coverage(segment_line))
        for segment_line in asm_graph.segments
    )
