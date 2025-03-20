"""Seed sequence creation module."""

from collections.abc import Iterator
from pathlib import Path

import gfapy  # type: ignore[import-untyped]

import pangebin.gene_density.input_output as gd_io
import pangebin.gfa.panassembly.segment as gfa_pan_segment


def from_gene_density(gene_density_file: Path) -> Iterator[str]:
    """Extract seed sequence from gene density file."""
    for sequence_id, gene_density in gd_io.from_file(gene_density_file):
        if gene_density > 0:
            yield sequence_id


def from_contigs_to_fragment_seeds(
    unicycler_seeds: Iterator[str],
    skesa_seeds: Iterator[str],
    panassembly_gfa: gfapy.Gfa,
) -> Iterator[str]:
    """Iterate over fragments contained in seed contigs."""
    seed_contigs = set(unicycler_seeds) | set(skesa_seeds)
    for fragment_line in panassembly_gfa.segments:
        fragment_contig_list = gfa_pan_segment.contig_list(fragment_line)
        if any((contig in seed_contigs) for contig in fragment_contig_list):
            yield fragment_line.name
