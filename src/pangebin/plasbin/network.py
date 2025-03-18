"""Plasbin network graph module."""

from __future__ import annotations

import logging
from itertools import chain
from typing import TYPE_CHECKING

import pangebin.gc_content.items as gc_items
import pangebin.gfa.link as gfa_link
import pangebin.gfa.panassembly.path as gfa_pan_path
import pangebin.gfa.panassembly.segment as gfa_pan_segment
import pangebin.gfa.path as gfa_path
import pangebin.gfa.segment as gfa_segment

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    import gfapy  # type: ignore[import-untyped]

_LOGGER = logging.getLogger(__name__)


class Network:
    """Plasbin network graph."""

    SOURCE_VERTEX = "S"
    SINK_VERTEX = "T"

    def __init__(
        self,
        panasm_graph: gfapy.Gfa,
        seeds: Iterable[str],
        gc_scores: Iterable[gc_items.SequenceProbabilityScores],
        plasmidness: Iterable[tuple[str, float]],
    ) -> None:
        """Initialize a plasbin network graph.

        Pan-assembly graph enriched with:

        * Source vertex
        * Sink vertex
        * Links from the source to seeds and from fragments to the sink
        """
        self.__panasm_graph: gfapy.Gfa = panasm_graph
        self.__seeds = set(seeds)
        self.__coverages: dict[str, float] = dict(self.__init_coverages())
        self.__gc_scores = {
            frag_gc_score.sequence_id(): frag_gc_score.probability_scores()
            for frag_gc_score in gc_scores
        }
        self.__plasmidness = dict(plasmidness)

    def number_of_vertices(self) -> int:
        """Get number of vertices.

        :math:`|V| = 2 x |Fragments| + 2`
        """
        return 2 * len(self.__panasm_graph.segment_names) + 2

    def fragment_ids(self) -> Iterator[str]:
        """Get fragment ids."""
        return self.__panasm_graph.segment_names

    def oriented_fragments(self) -> Iterator[gfa_segment.OrientedFragment]:
        """Create oriented fragments."""
        return (
            fragment
            for oriented_fragments in gfa_segment.gfa_oriented_fragments(
                self.__panasm_graph,
            )
            for fragment in chain(oriented_fragments)
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
        for link in gfa_link.gfa_links(self.__panasm_graph):
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
    def __init_coverages(self) -> Iterator[tuple[str, float]]:
        for segment in self.__panasm_graph.segments:
            yield (
                segment.name,
                max(
                    gfa_pan_path.normalized_coverage(
                        gfa_path.get_path_line_by_name(self.__panasm_graph, ctg_id),
                    )
                    for ctg_id in gfa_pan_segment.contig_list(segment)
                ),
            )

    def panasm_graph(self) -> gfapy.Gfa:
        """Get pan-assembly graph."""
        return self.__panasm_graph

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
            self.__panasm_graph.segment(fragment_id).disconnect()
            if fragment_id in self.__seeds:
                self.__seeds.remove(fragment_id)
            del self.__coverages[fragment_id]
            del self.__gc_scores[fragment_id]
            del self.__plasmidness[fragment_id]
        else:
            self.__coverages[fragment_id] -= coverage_to_substrack
