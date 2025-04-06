"""PangeBin-flow MILP objectives."""

from collections.abc import Callable, Iterable, Iterator
from enum import StrEnum

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
