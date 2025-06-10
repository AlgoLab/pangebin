"""Create standardized files."""

import logging

import gfapy  # type: ignore[import-untyped]

import pangebin.assembly.items as asm_items
import pangebin.gfa.ops as gfa_ops
from pangebin.gfa.assembler import segment as gfa_asm_segment
from pangebin.std_asm_graph.config import Config

_LOGGER = logging.getLogger(__name__)


def standardize_assembly_graph(
    graph: gfapy.Gfa,
    assembler_id: asm_items.Identifier,
    config: Config,
) -> None:
    """Standardize assembly graph.

    * Add `LN` tag to segments (sequence lengths)
    * Rename segments according to the assembler
    * Set (if not exist) `dp` segment tag (normalized coverage)
    * Set header tag `Sd` to `Y` (gfa is standardized)
    * Transform small contigs into links

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph
    assembler_id : asm_items.Identifier
        Assembler identifier
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

    gfa_ops.transform_small_contigs_into_links(graph, config.min_contig_length())

    gfa_ops.set_standardized_header_tag(graph)

    # TODO get map between original and new segment names
