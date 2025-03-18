"""Seed sequence creation module."""

from collections.abc import Iterator
from pathlib import Path

import pangebin.gene_density.input_output as gd_io


def from_gene_density(gene_density_file: Path) -> Iterator[str]:
    """Extract seed sequence from gene density file."""
    for sequence_id, gene_density in gd_io.from_file(gene_density_file):
        if gene_density > 0:
            yield sequence_id
