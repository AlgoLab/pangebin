"""PangeBin-flow MILP objectives."""

from collections.abc import Callable, Iterable, Iterator
from enum import StrEnum

import pangebin.plasbin.network as net

# REFACTOR to different class
# REFACTOR move to cli the choice of a final class


class ObjectiveFunctionDomain(StrEnum):
    """Objective function domain."""

    ALL = "all"
    SEEDS = "seeds"

    def to_fn(self) -> Callable[[net.Network], Iterator[str] | Iterable[str]]:
        """Get objective function domain function."""
        match self:
            case ObjectiveFunctionDomain.ALL:
                # XXX PAER TESTS
                # return all_but_list
                return net.Network.fragment_ids
            case ObjectiveFunctionDomain.SEEDS:
                return net.Network.seeds
            case _:
                raise ValueError


def all_but_list(network: net.Network, refused: Iterable[str]) -> Iterator[str]:
    """Get all but list."""
    return (frag_id for frag_id in network.fragment_ids() if frag_id not in refused)


def zeta_i(network: net.Network, frag_id: str, max_frag_length: int) -> float:
    """Get zeta_i coefficient."""
    return net.length(network, frag_id) / max_frag_length
