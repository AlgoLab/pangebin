"""PangeBin-flow MILP objectives."""

from collections.abc import Callable, Iterable, Iterator
from enum import StrEnum

import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.network as net


class ObjectiveFunctionDomain(StrEnum):
    """Objective function domain."""

    ALL = "all"
    SEEDS = "seeds"

    def to_fn(self) -> Callable[[net.Network], Iterator[str] | Iterable[str]]:
        """Get objective function domain function."""
        match self:
            case ObjectiveFunctionDomain.ALL:
                return net.Network.fragment_ids
            case ObjectiveFunctionDomain.SEEDS:
                return net.Network.seeds
            case _:
                raise ValueError


def max_frag_length(
    network: net.Network,
    frag_set_fn: Callable[[net.Network], Iterable[str]],
) -> int:
    """Get maximum fragment length."""
    return max(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in frag_set_fn(network)
    )


def zeta_i(network: net.Network, frag_id: str, max_frag_length: int) -> float:
    """Get zeta_i coefficient."""
    return gfa_segment.length(network.gfa_graph().segment(frag_id)) / max_frag_length
