"""Plasmidness creation module."""

from collections.abc import Iterable

import gfapy  # type: ignore[import-untyped]

import pangebin.assembly.items as asm_items
import pangebin.gfa.panassembly.path as gfa_pan_path
import pangebin.gfa.panassembly.segment as gfa_pan_segment
import pangebin.gfa.segment as gfa_segment


def from_contigs_to_fragments_plasmidness(
    skesa_plasmidness: Iterable[tuple[str, float]],
    unicycler_plasmidness: Iterable[tuple[str, float]],
    panasm_graph: gfapy.Gfa,
) -> Iterable[tuple[str, float]]:
    """Create plasmidness from contigs to fragments."""
    d_ske_plms = dict(skesa_plasmidness)
    d_uni_plms = dict(unicycler_plasmidness)
    for frag_id in panasm_graph.segment_names:
        frag_len = gfa_segment.length(frag_id)
        len_normalized_sum_plasmidness = 0.0
        frag_ctg_list = gfa_pan_segment.contig_list(frag_id)
        for contig in frag_ctg_list:
            contig_len = gfa_pan_path.length(contig)
            match gfa_pan_path.name_to_assembler_id(contig):
                case asm_items.Identifier.SKESA:
                    len_normalized_sum_plasmidness += (
                        d_ske_plms[contig] * frag_len / contig_len
                    )
                case asm_items.Identifier.UNICYCLER:
                    len_normalized_sum_plasmidness += (
                        d_uni_plms[contig] * frag_len / contig_len
                    )
        yield (frag_id, len_normalized_sum_plasmidness / len(frag_ctg_list))
