"""Multi-flow managers."""

from __future__ import annotations

from abc import abstractmethod
from itertools import product
from typing import TYPE_CHECKING, Self

import gurobipy as gp

import pangebin.plasbin.network as net

from . import bins
from . import config as cfg
from .milp import config as lp_cfg
from .milp import constraints as lp_csts
from .milp import models as lp_models

if TYPE_CHECKING:
    from collections.abc import Iterator


class Stats:
    """Stats container."""

    @classmethod
    def new(
        cls,
        max_number_of_bins: int,
    ) -> Stats:
        """Create new stats container."""
        return cls(0, max_number_of_bins)

    def __init__(
        self,
        number_of_active_bins: int,
        max_number_of_bins: int,
    ) -> None:
        self.__number_of_active_bins = number_of_active_bins
        self.__max_number_of_bins = max_number_of_bins

    def number_of_active_bins(self) -> int:
        """Get number of active bins."""
        return self.__number_of_active_bins

    def max_number_of_bins(self) -> int:
        """Get maximum number of bins."""
        return self.__max_number_of_bins

    def add_bin(self) -> None:
        """Add seed bin."""
        self.__number_of_active_bins += 1


class BinClass:
    """Base manager."""

    @classmethod
    @abstractmethod
    def topology(cls) -> bins.Topology:
        """Get topology."""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def seed_constraint(cls) -> bins.SeedConstraint:
        """Get seed constraint."""
        raise NotImplementedError

    @classmethod
    def new(
        cls,
        network: net.Network,
        model_config: lp_cfg.Config,
        max_number_of_bins: int,
        topology: bins.Topology,
        seed_constraint: bins.SeedConstraint,
    ) -> Self:
        """Create new circular manager."""
        model = lp_models.Model(
            network,
            model_config,
            max_number_of_bins,
            topology,
            seed_constraint,
        )
        return cls(model, Stats.new(max_number_of_bins))

    def __init__(
        self,
        model: lp_models.Model,
        stats: Stats,
    ) -> None:
        self.__model = model
        self.__stats = stats
        self.__obj_lb_cst: gp.Constr | None = None

    def model(self) -> lp_models.Model:
        """Get HMF model."""
        return self.__model

    def stats(self) -> Stats:
        """Get stats."""
        return self.__stats

    def new_bin(self) -> int:
        """Add new bin."""
        new_mf_idx = self.__stats.number_of_active_bins()
        self.__model.activate_flow(new_mf_idx)
        return new_mf_idx

    def set_objective_lb(self, lower_bound: float, network: net.Network) -> None:
        """Set objective lower bound."""
        if self.__obj_lb_cst is None:
            # FIXME specialize for PC
            self.__obj_lb_cst = lp_csts.objective_lower_bound(
                self.__model.gurobi_model(),
                lower_bound,
            )
        else:
            # FIXME specialize for PC
            lp_csts.update_binning_objective_lower_bound(
                self.__model.gurobi_model(),
                self.__obj_lb_cst,
                lower_bound,
            )


class _Circular(BinClass):
    """Circular multi-flows manager."""

    @classmethod
    def topology(cls) -> bins.Topology:
        """Get topology."""
        return bins.Topology.CIRCULAR

    def _define_topology(self) -> None:
        """Define topology."""
        self.__model.state_constraints()[
            self.__stats.number_of_active_bins()
        ].define_circular()


def iter_bin_class_manager(
    bin_properties: cfg.BinProperties,
    network: net.Network,
    model_config: lp_cfg.Config,
) -> Iterator[BinClass]:
    """Iterate over bin class managers."""
    for topology, seed_constraint in product(
        (bins.Topology.CIRCULAR, bins.Topology.PARTIALLY_CIRCULAR),
        (bins.SeedConstraint.REQUIRED, bins.SeedConstraint.NOT_REQUIRED),
    ):
        match seed_constraint:
            case bins.SeedConstraint.REQUIRED:
                max_number_of_bins = (
                    len(network.seeds())
                    if bin_properties.seed_constraint_for_topology(
                        topology,
                    ).with_seeds()
                    else 0
                )
            case bins.SeedConstraint.NOT_REQUIRED:
                max_number_of_bins = min(
                    bin_properties.seed_constraint_for_topology(
                        topology,
                    ).max_number_free_of_seeds(),
                    network.number_of_fragments(),
                )

        if max_number_of_bins > 0:
            yield BinClass.new(
                network,
                model_config,
                max_number_of_bins,
                topology,
                seed_constraint,
            )
