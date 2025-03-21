"""Plasbin network graph module."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pangebin.gc_content.items as gc_items
import pangebin.gfa.assembler.ops as gfa_asm_ops
import pangebin.gfa.link as gfa_link
import pangebin.gfa.panassembly.ops as gfa_pan_ops
import pangebin.gfa.segment as gfa_segment

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    import gfapy  # type: ignore[import-untyped]

_LOGGER = logging.getLogger(__name__)


class Network:
    """Plasbin network graph."""

    SOURCE_VERTEX = "S"
    SINK_VERTEX = "T"

    @classmethod
    def from_asm_graph(
        cls,
        asm_graph: gfapy.Gfa,
        seed_contigs: Iterable[str],
        contig_gc_scores: Iterable[gc_items.SequenceGCScores],
        contig_plasmidness: Iterable[tuple[str, float]],
    ) -> Network:
        """Create plasbin network graph from a standardized assembly graph."""
        return cls(
            asm_graph,
            seed_contigs,
            gfa_asm_ops.contig_coverages(asm_graph),
            contig_gc_scores,
            contig_plasmidness,
        )

    @classmethod
    def from_panasm_graph(
        cls,
        panasm_graph: gfapy.Gfa,
        seed_fragments: Iterable[str],
        fragment_gc_scores: Iterable[gc_items.SequenceGCScores],
        fragment_plasmidness: Iterable[tuple[str, float]],
    ) -> Network:
        """Create plasbin network graph from pan-assembly graph."""
        return cls(
            panasm_graph,
            seed_fragments,
            gfa_pan_ops.fragment_max_contig_coverages(panasm_graph),
            fragment_gc_scores,
            fragment_plasmidness,
        )

    def __init__(
        self,
        gfa_graph: gfapy.Gfa,
        seeds: Iterable[str],
        coverages: Iterable[tuple[str, float]],
        gc_scores: Iterable[gc_items.SequenceGCScores],
        plasmidness: Iterable[tuple[str, float]],
    ) -> None:
        """Initialize a plasbin network graph.

        GFA graph enriched with:

        * Source vertex
        * Sink vertex
        * Links from the source to seeds and from fragments to the sink
        """
        self.__gfa_graph: gfapy.Gfa = gfa_graph
        self.__seeds = set(seeds)
        self.__coverages = dict(coverages)
        self.__gc_scores = {
            frag_gc_score.sequence_id(): frag_gc_score.probability_scores()
            for frag_gc_score in gc_scores
        }
        self.__plasmidness = dict(plasmidness)

    def number_of_fragments(self) -> int:
        """Get number of fragments."""
        return len(self.__gfa_graph.segment_names)

    def number_of_vertices(self) -> int:
        """Get number of vertices.

        :math:`|V| = 2 x |Fragments| + 2`
        """
        return 2 * self.number_of_fragments() + 2

    def fragment_ids(self) -> Iterator[str]:
        """Get fragment ids."""
        return self.__gfa_graph.segment_names

    def oriented_fragments(self) -> Iterator[gfa_segment.OrientedFragment]:
        """Create oriented fragments."""
        return (
            fragment
            for oriented_fragments in gfa_segment.gfa_oriented_fragments(
                self.__gfa_graph,
            )
            for fragment in oriented_fragments
        )

    def source_arcs(
        self,
    ) -> Iterator[tuple[str, gfa_segment.OrientedFragment]]:
        """Create source to seed arcs."""
        for seed in self.__seeds:
            yield (
                self.SOURCE_VERTEX,
                gfa_segment.OrientedFragment(seed, gfa_segment.Orientation.FORWARD),
            )
            yield (
                self.SOURCE_VERTEX,
                gfa_segment.OrientedFragment(seed, gfa_segment.Orientation.REVERSE),
            )

    def link_arcs(self) -> Iterator[gfa_link.Link]:
        """Create arcs between fragments."""
        for link in gfa_link.gfa_links(self.__gfa_graph):
            yield link
            rev_link = link.to_reverse()
            if rev_link.predecessor() != link.predecessor():
                yield rev_link

    def sink_arcs(
        self,
    ) -> Iterator[tuple[gfa_segment.OrientedFragment, str]]:
        """Create fragment to sink arcs."""
        return ((fragment, self.SINK_VERTEX) for fragment in self.oriented_fragments())

    #
    # Attributes
    #
    def gfa_graph(self) -> gfapy.Gfa:
        """Get the GFA graph."""
        return self.__gfa_graph

    def seeds(self) -> set[str]:
        """Get seeds."""
        return self.__seeds

    def coverage(self, fragment_id: str) -> float:
        """Get coverage."""
        return self.__coverages[fragment_id]

    def gc_score(self, fragment_id: str) -> list[float]:
        """Get GC score."""
        return self.__gc_scores[fragment_id]

    def plasmidness(self, fragment_id: str) -> float:
        """Get plasmidness."""
        return self.__plasmidness[fragment_id]

    def cap_s(self, source_arc: tuple[str, gfa_segment.OrientedFragment]) -> float:
        """Get capacity for source arc."""
        return self.coverage(source_arc[1].identifier())

    def cap(self, arc: gfa_link.Link) -> float:
        """Get capacity for regular arc."""
        return min(
            self.coverage(arc.predecessor().identifier()),
            self.coverage(arc.successor().identifier()),
        )

    def cap_t(self, sink_arc: tuple[gfa_segment.OrientedFragment, str]) -> float:
        """Get capacity for sink arc."""
        return self.coverage(sink_arc[0].identifier())

    def reduce_coverage(self, fragment_id: str, coverage_to_substrack: float) -> None:
        """Reduce coverage."""
        if coverage_to_substrack > self.coverage(fragment_id):
            _crit_msg = (
                f"Coverage to substrack {coverage_to_substrack}"
                f" > fragment coverage {self.coverage(fragment_id)}"
            )
            _LOGGER.critical(_crit_msg)
            raise ValueError(_crit_msg)

        if coverage_to_substrack == self.coverage(fragment_id):
            # XXX contig depedending of the frag are also removed, problem?
            self.__gfa_graph.segment(fragment_id).disconnect()
            if fragment_id in self.__seeds:
                self.__seeds.remove(fragment_id)
            del self.__coverages[fragment_id]
            del self.__gc_scores[fragment_id]
            del self.__plasmidness[fragment_id]
        else:
            self.__coverages[fragment_id] -= coverage_to_substrack
