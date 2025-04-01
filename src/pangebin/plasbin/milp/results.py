"""MILP results module.

After solving a model, and if the optimal solution is found,
gurobipy adds the X attribute.
However, gurobi seems to free this attribute while it should not,
for example after yielding the solution.

This module unfortunately copy these values but ensure their validity.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

_LOGGER = logging.getLogger(__name__)


def active_fragments(
    network: net.Network,
    frag_vars: pb_lp_var.SubFragments,
) -> Iterator[str]:
    """Get active fragments."""
    return (f_id for f_id in network.fragment_ids() if frag_vars.x(f_id).X > 0)


def active_gc_content_interval(
    intervals: gc_items.Intervals,
    gc_vars: pb_lp_var.GCIntervals,
) -> tuple[float, float]:
    """Get active GC content interval."""
    for interval in intervals:
        if gc_vars.gc(interval).X > 0:
            return interval
    _crt_msg = "Could not find active GC content interval"
    _LOGGER.critical(_crt_msg)
    raise ValueError(_crt_msg)


class Pangebin:
    """Pangebin results."""

    @classmethod
    def from_optimal_variables(
        cls,
        network: net.Network,
        intervals: gc_items.Intervals,
        flow_vars: pb_lp_var.Flow,
        frag_vars: pb_lp_var.SubFragments,
        gc_vars: pb_lp_var.GCIntervals,
    ) -> Pangebin:
        """Get result from variable values."""
        return cls(
            (
                (
                    frag_id,
                    flow_vars.incoming_forward_reverse(network, frag_id).getValue(),
                )
                for frag_id in active_fragments(network, frag_vars)
            ),
            flow_vars.total().X,
            sum(
                gfa_segment.length(network.gfa_graph().segment(frag_id))
                for frag_id in active_fragments(network, frag_vars)
            ),
            active_gc_content_interval(intervals, gc_vars),
        )

    def __init__(
        self,
        fragments_incoming_flow: Iterable[tuple[str, float]],
        total_flow: float,
        cumulative_length: int,
        gc_interval: tuple[float, float],
    ) -> None:
        self.__frags_incoming_flow = dict(fragments_incoming_flow)
        self.__total_flow = total_flow
        self.__cumulative_length = cumulative_length
        self.__gc_interval = gc_interval

    def fragments_incoming_flow(self) -> Iterator[tuple[str, float]]:
        """Get fragments incoming flow."""
        yield from self.__frags_incoming_flow.items()

    def fragment_incoming_flow(self, frag_id: str) -> float:
        """Get fragment incoming flow."""
        return self.__frags_incoming_flow[frag_id]

    def number_of_active_fragments(self) -> int:
        """Get number of active fragments."""
        return len(self.__frags_incoming_flow)

    def active_fragments(self) -> Iterator[str]:
        """Get active fragments."""
        yield from self.__frags_incoming_flow

    def total_flow(self) -> float:
        """Get total flow value."""
        return self.__total_flow

    def cumulative_length(self) -> int:
        """Get cumulative length."""
        return self.__cumulative_length

    def gc_interval(self) -> tuple[float, float]:
        """Get GC interval."""
        return self.__gc_interval
