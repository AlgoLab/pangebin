"""MILP results module.

After solving a model, and if the optimal solution is found,
gurobipy adds the X attribute.
However, gurobi seems to free this attribute while it should not,
for example after yielding the solution.

This module unfortunately copy these values but ensure their validity.
"""

from __future__ import annotations

import logging
from contextlib import suppress
from typing import TYPE_CHECKING

import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.network as pb_network

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

_LOGGER = logging.getLogger(__name__)


# TODO move all result method here from milp_vars module


def active_fragments(
    network: pb_network.Network,
    variables: milp_vars.MaxCovFlow,
) -> Iterator[str]:
    """Get active fragments."""
    # BUG tmp print
    _number_of_active_fragments = 0
    for frag_id in network.fragment_ids():
        is_active = False
        frag_f = gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.FORWARD)
        frag_r = gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.REVERSE)
        # OPTIMIZE can I use a conditionnal rather than catching the exception?
        with suppress(AttributeError):
            if variables.x(frag_f).X > 0:
                is_active = True
        with suppress(AttributeError):
            if variables.x(frag_r).X > 0:
                is_active = True
        if is_active:
            _number_of_active_fragments += 1
            yield frag_id
    _LOGGER.debug("Number of active fragments: %d", _number_of_active_fragments)


class MaxCovFlow:
    """Maximum Coverage Flow results."""

    @classmethod
    def from_optimal_variables(
        cls,
        network: pb_network.Network,
        var: milp_vars.MaxCovFlow,
    ) -> MaxCovFlow:
        """Get result from variable values."""
        return cls(
            (
                (
                    frag_id,
                    milp_vars.incoming_flow_forward_reverse(
                        frag_id,
                        network,
                        var,
                    ).getValue(),
                )
                for frag_id in active_fragments(network, var)
            ),
            var.total_flow().X,
        )

    def __init__(
        self,
        fragments_incoming_flow: Iterable[tuple[str, float]],
        total_flow: float,
    ) -> None:
        self.__frags_incoming_flow = dict(fragments_incoming_flow)
        self.__total_flow = total_flow

    def fragments_incoming_flow(self) -> Iterator[tuple[str, float]]:
        """Get fragments incoming flow."""
        yield from self.__frags_incoming_flow.items()

    def fragment_incoming_flow(self, frag_id: str) -> float:
        """Get fragment incoming flow."""
        return self.__frags_incoming_flow[frag_id]

    def active_fragments(self) -> Iterator[str]:
        """Get active fragments."""
        yield from self.__frags_incoming_flow

    def total_flow(self) -> float:
        """Get total flow value."""
        return self.__total_flow
