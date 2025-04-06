"""Plasbin network graph module."""

from __future__ import annotations

import contextlib
from enum import StrEnum
from itertools import chain
from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]

import pangebin.gc_content.items as gc_items
import pangebin.gfa.assembler.ops as gfa_asm_ops
import pangebin.gfa.link as gfa_link
import pangebin.gfa.panassembly.ops as gfa_pan_ops
import pangebin.gfa.path as gfa_path
import pangebin.gfa.segment as gfa_segment

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator


class SinkArcsDomain(StrEnum):
    """Sink-arcs definition."""

    ALL = "all"
    SEEDS = "seeds"


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
        sink_arcs_definition: SinkArcsDomain,
    ) -> Network:
        """Create plasbin network graph from a standardized assembly graph."""
        return cls(
            asm_graph,
            seed_contigs,
            gfa_asm_ops.contig_coverages(asm_graph),
            contig_gc_scores,
            contig_plasmidness,
            sink_arcs_definition,
        )

    @classmethod
    def from_panasm_graph(
        cls,
        panasm_graph: gfapy.Gfa,
        seed_fragments: Iterable[str],
        fragment_gc_scores: Iterable[gc_items.SequenceGCScores],
        fragment_plasmidness: Iterable[tuple[str, float]],
        sink_arcs_definition: SinkArcsDomain,
    ) -> Network:
        """Create plasbin network graph from pan-assembly graph."""
        return cls(
            panasm_graph,
            seed_fragments,
            gfa_pan_ops.fragment_max_contig_coverages(panasm_graph),
            fragment_gc_scores,
            fragment_plasmidness,
            sink_arcs_definition,
        )

    def __init__(  # noqa: PLR0913
        self,
        gfa_graph: gfapy.Gfa,
        seeds: Iterable[str],
        coverages: Iterable[tuple[str, float]],
        gc_scores: Iterable[gc_items.SequenceGCScores],
        plasmidness: Iterable[tuple[str, float]],
        sink_arcs_definition: SinkArcsDomain,
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
        self.__sink_arc_gen_fn = (
            self.fragment_ids
            if sink_arcs_definition == SinkArcsDomain.ALL
            else self.seeds
        )
        self.__is_sink_connected_fn = (
            (lambda _: True)
            if sink_arcs_definition == SinkArcsDomain.ALL
            else (lambda frag_id: frag_id in self.__seeds)
        )

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

    def to_oriented(
        self,
        fragment_id: str,
    ) -> tuple[gfa_segment.OrientedFragment, gfa_segment.OrientedFragment]:
        """Give the two oriented fragments corresponding to the fragment id.

        Returns
        -------
        OrientedFragment
            Forward fragment
        OrientedFragment
            Reverse fragment
        """
        return (
            gfa_segment.OrientedFragment(fragment_id, gfa_segment.Orientation.FORWARD),
            gfa_segment.OrientedFragment(fragment_id, gfa_segment.Orientation.REVERSE),
        )

    def is_source_connected(self, fragment_id: str) -> bool:
        """Check if the source goes into the corresponding oriented fragment."""
        return fragment_id in self.__seeds

    def is_sink_connected(self, fragment_id: str) -> bool:
        """Check if the corresponding oriented fragment goes into the sink."""
        return self.__is_sink_connected_fn(fragment_id)

    def source_connected_fragment_ids(self) -> Iterator[str]:
        """Get source connected fragment ids."""
        return iter(self.__seeds)

    def sink_connected_fragment_ids(self) -> Iterator[str] | Iterable[str]:
        """Get sink connected fragment ids."""
        return self.__sink_arc_gen_fn()

    def source_arcs(
        self,
    ) -> Iterator[tuple[str, gfa_segment.OrientedFragment]]:
        """Create source to seed arcs."""
        return (
            (self.SOURCE_VERTEX, frag)
            for seed in self.__seeds
            for frag in self.to_oriented(seed)
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
        """Create fragment to sink-arcs."""
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
        if coverage_to_substrack >= self.coverage(fragment_id):
            # XXX All the references of a fragment is not deleted
            # see https://github.com/ggonnella/gfapy/issues/33
            #
            # Remove path
            #
            path_line: gfa_path.GfaPath
            # Use list constructor to avoid concurrent modification
            for path_line in list(self.__gfa_graph.paths):
                if fragment_id in (
                    seg_line.name for seg_line in path_line.segment_names
                ):
                    path_line.disconnect()
            #
            # Remove links
            #
            link_line: gfa_link.GfaLink
            for link_line in list(
                gfa_segment.get_segment_line_by_name(
                    self.__gfa_graph,
                    fragment_id,
                ).dovetails,
            ):
                with contextlib.suppress(gfapy.RuntimeError):
                    link_line.disconnect()  # FIXME potential bug https://github.com/ggonnella/gfapy/issues/36
            self.__gfa_graph.segment(fragment_id).disconnect()
            if fragment_id in self.__seeds:
                self.__seeds.remove(fragment_id)
            del self.__coverages[fragment_id]
            del self.__gc_scores[fragment_id]
            del self.__plasmidness[fragment_id]
        else:
            self.__coverages[fragment_id] -= coverage_to_substrack


class StrFormatter:
    """Network string formatter."""

    @classmethod
    def arc(cls, arc: gfa_link.Link) -> str:
        """Get arc string."""
        return str(arc)

    @classmethod
    def s_arc(cls, link: tuple[str, gfa_segment.OrientedFragment]) -> str:
        """Get source arc string."""
        return f"{link[0]}_{link[1]}"

    @classmethod
    def t_arc(cls, link: tuple[gfa_segment.OrientedFragment, str]) -> str:
        """Get sink arc string."""
        return f"{link[0]}_{link[1]}"

    @classmethod
    def arc_ids(cls, network: Network) -> Iterator[str]:
        """Get iterator over source, arc and sink arc ids."""
        return chain(
            (cls.s_arc(s_link) for s_link in network.source_arcs()),
            (cls.arc(arc) for arc in network.link_arcs()),
            (cls.t_arc(t_link) for t_link in network.sink_arcs()),
        )
