"""MILP model."""

from __future__ import annotations

import gurobipy as gp

import pangebin.plasbin.hmf.bins as hmf_bins
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net

from . import config as lp_cfg
from . import constraints as lp_csts
from . import objectives as lp_objs
from . import variables as lp_vars


class Model:
    """HMF base model."""

    def __init__(
        self,
        network: net.Network,
        config: lp_cfg.Config,
        max_number_of_bins: int,
        topology: hmf_bins.Topology,
        seed_constraint: hmf_bins.SeedConstraint,
    ) -> None:
        self.__config = config
        self.__topology = topology
        self.__seed_constraint = seed_constraint
        self.__model = gp.Model("Hierarchical Multi Flow")
        self.__bins_variables, self.__flow_union_frag_vars = lp_vars.init_binning(
            network,
            self.__model,
            max_number_of_bins,
        )
        self.__state_constraints = lp_csts.set_constraints(
            self.__model,
            self.__bins_variables,
            self.__flow_union_frag_vars,
            network,
            config,
        )
        self.__model.setObjective(self._init_obj_linexpr(network), gp.GRB.MAXIMIZE)

    def config(self) -> lp_cfg.Config:
        """Config."""
        return self.__config

    def gurobi_model(self) -> gp.Model:
        """Gurobi model."""
        return self.__model

    def bins_variables(self) -> list[lp_vars.BinVariables]:
        """Binning variables."""
        return self.__bins_variables

    def flow_union_frag_vars(self) -> cmn_lp_vars.SubFragments:
        """Flow union fragment variables."""
        return self.__flow_union_frag_vars

    def state_constraints(self) -> list[lp_csts.BinStateConstraints]:
        """State constraints."""
        return self.__state_constraints

    def topology(self) -> hmf_bins.Topology:
        """Topology."""
        return self.__topology

    def seed_constraint(self) -> hmf_bins.SeedConstraint:
        """Seed constraint."""
        return self.__seed_constraint

    def activate_flow(self, index: int) -> None:
        """Activate the flow at the index."""
        self.__state_constraints[index].activate_bin()
        self.__state_constraints[index].define_topology(self.__topology)
        self.__state_constraints[index].define_seed_constraint(
            self.__seed_constraint,
        )

    def _init_obj_linexpr(self, network: net.Network) -> gp.LinExpr:
        """Objective."""
        match self.__topology:
            case hmf_bins.Topology.CIRCULAR:
                return lp_objs.circular_objective(
                    self.__bins_variables,
                    self.__flow_union_frag_vars,
                    network,
                    self.__config.obj_fun_domain(),
                )
            case hmf_bins.Topology.PARTIALLY_CIRCULAR:
                return lp_objs.partially_circular_objective(
                    self.__bins_variables,
                    self.__flow_union_frag_vars,
                    network,
                    self.__config.obj_fun_domain(),
                )
