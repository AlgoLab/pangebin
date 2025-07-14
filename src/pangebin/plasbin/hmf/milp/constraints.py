"""PangeBin-flow classbin MILP constraints."""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, cast, final, runtime_checkable

import gurobipy as gp

import pangebin.plasbin.hmf.bins as hmf_bins
import pangebin.plasbin.milp.connected_component.constraints as ccomp_cst
import pangebin.plasbin.milp.constraints as cmn_lp_cst
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net

from . import config as lp_cfg
from . import objectives as lp_obj
from . import variables as lp_vars

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

# REFACTOR use config dataclass for each model


# OPTIMIZE use common constraints/constants
# REFACTOR use wrapper for obj lower bound cst
def set_constraints(
    m: gp.Model,
    bins_vars: list[lp_vars.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
    config: lp_cfg.Config,
) -> list[BinStateConstraints]:
    """Set binning constraints."""
    m.addConstrs(
        (
            gp.quicksum(
                var.flows().incoming_forward_reverse(network, frag_id)
                for var in bins_vars
            )
            <= network.coverage(frag_id)
            for frag_id in network.gfa_graph().segment_names
        ),
        name="fragment_coverages_limit_cumulative_flows",
    )

    # DOCU define use of frag in at least one bin
    m.addConstrs(
        (
            bin_vars.sub_frag().frag(frag_id) <= flow_union_frag_vars.frag(frag_id)
            for bin_vars in bins_vars
            for frag_id in network.fragment_ids()
        ),
        name="active_union_flow_frag_if_frag_in_bin",
    )
    m.addConstrs(
        (
            flow_union_frag_vars.frag(frag_id)
            <= gp.quicksum(bin_vars.sub_frag().frag(frag_id) for bin_vars in bins_vars)
            for frag_id in network.fragment_ids()
        ),
        name="frag_active_in_at_least_one_bin_if_active_in_any_bin",
    )

    multi_flow_state_csts = []
    for b, var in enumerate(bins_vars):
        constraints, bin_state_csts = _set_constraints_for_one_bin(
            m,
            var,
            network,
            config,
        )
        for const in constraints:
            const.ConstrName = f"bin_{b}_{const.ConstrName}"

        m.update()
        bin_state_csts.deactivate_bin()

        multi_flow_state_csts.append(bin_state_csts)

    return multi_flow_state_csts


def _set_constraints_for_one_bin(
    m: gp.Model,
    var: lp_vars.BinVariables,
    network: net.Network,
    config: lp_cfg.Config,
) -> tuple[list[gp.Constr], BinStateConstraints]:
    """Set MGCLB constraints."""
    constraints: list[gp.Constr] = []
    structural_constraints: list[StructuralCst] = []

    constraints += cmn_lp_cst.active_fragments_active_one_of_their_orientations(
        m,
        var.sub_vertices(),
        var.sub_frag(),
        network,
    )

    structural_constraints.append(SeedsLB(m, var, network))

    #
    # Source-sink constraints
    #
    structural_constraints.append(SourceArcsUB(m, var, network))
    structural_constraints.append(SinkArcsUB(m, var, network))

    structural_constraints.append(CycleBeforeOut(m, var, network))
    structural_constraints.append(CycleBeforeIn(m, var, network))

    # DOCU trick on parcirc cst for s and t connected repeated frag
    # In fact, if we begin with a repeat looping at the beginning,
    # we cannot loop over it because the beginning is connected to the source
    # However, we can out at the end of the loop with T
    # * PC: no other incoming link if source connected
    structural_constraints.append(SourceVerticesInArcUB(m, var, network))
    # * PC: no other outgoing link if sink connected
    structural_constraints.append(SinkVerticesOutArcUB(m, var, network))

    # Two inequalities to change
    structural_constraints.append(STConnexionLB(m, var, network))
    structural_constraints.append(STConnexionUB(m, var, network))

    structural_constraints.append(NonSeedSourceArcsUB(m, var, network))
    structural_constraints.append(NonSeedSinkArcsUB(m, var, network))

    #
    # Flow constraints
    #
    # DOCU LATEX CONTINUE HERE
    constraints += cmn_lp_cst.flow_conservation(m, var.flows(), network)
    constraints += cmn_lp_cst.total_flow_value(m, var.flows(), network)
    constraints += cmn_lp_cst.active_arcs_implies_active_fragments(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    constraints += cmn_lp_cst.active_fragments_imply_at_least_one_active_arc(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )

    # DOCU let flow be 0 on source and sink arc when circular bin
    # * C: the flow value on source and sink arcs equal is 0 AND y (must) can still be 1

    structural_constraints.append(
        ActiveSourceArcFlowLB(m, var, network, config.min_flow()),
    )
    structural_constraints.append(SourceArcFlowUB(m, var, network))
    constraints += cmn_lp_cst.active_link_arc_has_strict_positive_flow(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
        config.min_flow(),
    )
    structural_constraints.append(
        ActiveSinkArcFlowLB(m, var, network, config.min_flow()),
    )
    structural_constraints.append(SinkArcFlowUB(m, var, network))

    constraints += cmn_lp_cst.arc_capacities_limit_arc_flows(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
    )
    constraints += cmn_lp_cst.fragment_coverages_limit_cumulative_flows(
        m,
        var.flows(),
        network,
    )

    #
    # Connectivity
    #
    constraints += ccomp_cst.arcs_in_tree_are_active(
        m,
        var.sub_arcs(),
        var.tree_edges().dtree(),
        network,
    )

    # DOCU C: (opti) force each rev beta to be 0
    structural_constraints.append(BetaRevUB(m, var, network))

    constraints += ccomp_cst.rev_link_arcs_in_tree_are_active(
        m,
        var.sub_arcs(),
        var.tree_edges(),
        network,
    )

    # DOCU change RHS weither bin is active or not
    structural_constraints.append(SourceSubtreeSize(m, var, network))
    # DOCU change RHS weither bin is active or not
    structural_constraints.append(SinkSubtreeSizeRevVersion(m, var, network))

    constraints += ccomp_cst.subtree_size_fragment_rev_version(
        m,
        var.sub_vertices(),
        var.tree_edges(),
        network,
    )

    # DOCU C: [opti] only the seeds can be connected to the source and the sink
    structural_constraints.append(NonSeedRootUB(m, var, network))

    # DOCU change RHS weither bin is active or not
    structural_constraints.append(NumberOfRootVertices(m, var, network))

    constraints += ccomp_cst.at_most_one_source_arc_in_tree(
        m,
        var.tree_edges(),
        var.root(),
        network,
    )
    constraints += ccomp_cst.root_priority(m, var.root(), var.sub_arcs(), network)
    #
    # Plasmid property
    #
    # DOCU Plasmidness lower bound still correct if bin deactivate
    # BUG TMP REMOVE PLASMIDNESS LOWER BOUND CST
    constraints.append(
        _plasmidness_lower_bound(
            m,
            var.flows(),
            network,
            config.obj_fun_domain(),
            config.plasmidness_coefficient(),
        ),
    )
    # DOCU min cumulative len constraint weither bin active or not
    structural_constraints.append(
        CumulativeLengthLB(m, var, network, config.min_cumulative_len()),
    )

    for structural_cst in structural_constraints:
        constraints.extend(structural_cst.constraints())

    m.update()  # necessary to change the rhs
    # REFACTOR use compositions (and perhaps sub compositions): Protocol
    return constraints, BinStateConstraints(structural_constraints)


def _plasmidness_lower_bound(
    m: gp.Model,
    flow_vars: cmn_lp_vars.Flow,
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
    plasmidness_coefficient: float,
) -> gp.Constr:
    """Set plasmidness lower bound."""
    # DOCU plasmidness lower bound with plm coeff
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return m.addConstr(
        lp_obj.plasmidness_score(network, flow_vars, obj_fun_domain)
        >= gp.quicksum(
            lp_obj.frag_coeff(network, frag_id)
            * plasmidness_coefficient
            * flow_vars.incoming_forward_reverse(network, frag_id)
            for frag_id in frag_set_fn(network)
        ),
        name="plasmidness_lower_bound",
    )


def objective_lower_bound(
    m: gp.Model,
    obj_lower_bound: float,
) -> gp.Constr:
    """Set objective lower bound."""
    cst = m.addConstr(
        m.getObjective() >= obj_lower_bound,
        name="objective_function_lower_bound",
    )
    m.update()
    return cast("gp.Constr", cst)


def update_binning_objective_lower_bound(
    m: gp.Model,
    obj_lower_bound_cst: gp.Constr,
    new_obj_lower_bound: float,
) -> None:
    """Update binning objective lower bound."""
    # DOCU update objective lower bound
    m.update()
    obj_lower_bound_cst.RHS = new_obj_lower_bound
    m.update()


def plasmidness_score_lb_for_multi_flow_subset(
    m: gp.Model,
    bins_vars: list[lp_vars.BinVariables],
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
    from_until: tuple[int, int],
    last_circular_plasmidness_score: float,
) -> gp.Constr:
    """Set binning objective lower bound."""
    return m.addConstr(
        lp_obj.multi_flow_plasmidness_score(
            bins_vars,
            network,
            obj_fun_domain,
            from_until,
        )
        >= last_circular_plasmidness_score,
        name=f"plasmidness_score_lb_for_multi_flows_from_{from_until[0]}_until_{from_until[1]}",
    )


def plasmidness_score_order(
    m: gp.Model,
    bins_vars: list[lp_vars.BinVariables],
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
    bin_interval: tuple[int, int],
) -> list[gp.Constr]:
    """Set plasmidness score order for circular flows."""
    bin_interval_str = f"{bin_interval[0]}_{bin_interval[1]}"
    return list(
        m.addConstrs(
            (
                lp_obj.plasmidness_score(network, bins_vars[k].flows(), obj_fun_domain)
                >= lp_obj.plasmidness_score(
                    network,
                    bins_vars[k + 1].flows(),
                    obj_fun_domain,
                )
                for k in range(bin_interval[0], bin_interval[1] - 1)
            ),
            name=f"circular_flow_plasmidness_score_order_{bin_interval_str}",
        ).values(),
    )


def each_seed_must_be_in_at_least_one_bin(
    m: gp.Model,
    bins_vars: list[lp_vars.BinVariables],
    network: net.Network,
) -> list[gp.Constr]:
    """Set each seed must be in at least one bin."""
    return list(
        m.addConstrs(
            gp.quicksum(bin_vars.sub_frag().frag(seed_id) for bin_vars in bins_vars)
            >= 1
            for seed_id in network.seeds()
        ).values(),
    )


@runtime_checkable
class ContainsCsts(Protocol):
    """Class that can return constraints."""

    def constraints(self) -> Iterator[gp.Constr]:
        """Iterate over the constraints."""
        ...


@runtime_checkable
class Circular(Protocol, ContainsCsts):
    """Constraints that define a circular bin."""

    def define_circular(self) -> None:
        """Change the constraints to define a circular bin."""


@runtime_checkable
class PartiallyCircular(Protocol, ContainsCsts):
    """Constraints that define a partially circular bin."""

    def define_partially_circular(self) -> None:
        """Change the constraints to define a partially circular bin."""


@runtime_checkable
class Seed(Protocol, ContainsCsts):
    """Constraints that define a seed constraint."""

    def define_must_have_a_seed(self) -> None:
        """Change the constraints to define a must-have-a-seed bin."""


@runtime_checkable
class FreeOfSeed(Protocol, ContainsCsts):
    """Constraints that define a free-of-seed constraint."""

    def define_can_be_free_of_seed(self) -> None:
        """Change the constraints to define a free-of-seed bin."""


@runtime_checkable
class Activate(Protocol, ContainsCsts):
    """Constraints that define a seed constraint."""

    def activate(self) -> None:
        """Change the constraints to define a must-have-a-seed bin."""


@runtime_checkable
class Deactivate(Protocol, ContainsCsts):
    """Constraints that define a free-of-seed constraint."""

    def deactivate(self) -> None:
        """Change the constraints to define a free-of-seed bin."""


@runtime_checkable
class CircularWithSeed(Protocol, ContainsCsts):
    """Constraints that define a circular bin with a seed."""

    def define_circular_with_seed(self) -> None:
        """Change the constraints to define a circular bin."""
        ...

    def define_other(self) -> None:
        """Change the constraints to define a circular bin."""
        ...


StructuralCst = (
    Circular
    | PartiallyCircular
    | Seed
    | FreeOfSeed
    | Activate
    | Deactivate
    | CircularWithSeed
)


@final
class SeedsLB(Seed, FreeOfSeed, Deactivate):
    """Active seeds lower bound constraints (s + fs + deact)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__cst = cmn_lp_cst.active_seeds_lower_bound(m, var.sub_frag(), network, 0)

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield self.__cst

    def define_must_have_a_seed(self) -> None:
        """Define the bin as must-have-a-seed."""
        self.__cst.RHS = 1

    def define_can_be_free_of_seed(self) -> None:
        """Define the bin as free-of-seed."""
        self.__cst.RHS = 0

    def deactivate(self) -> None:
        """Deactivate the constraint."""
        self.__cst.RHS = 0


@final
class SourceArcsUB(Circular, PartiallyCircular, Deactivate):
    """Source arcs upper bound constraints (c + pc + deact)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__cst = cmn_lp_cst.source_arcs_upper_bound(m, var.sub_arcs(), network, 0)
        self.__number_of_source_arcs = sum(1 for _ in network.source_arcs())

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield self.__cst

    def define_circular(self) -> None:
        """Define the bin as circular."""
        self.__cst.RHS = 1

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        self.__cst.RHS = self.__number_of_source_arcs

    def deactivate(self) -> None:
        """Deactivate the constraint."""
        self.__cst.RHS = 0


@final
class SinkArcsUB(Circular, PartiallyCircular, Deactivate):
    """Sink arcs upper bound constraints (c + pc + deact)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__cst = cmn_lp_cst.sink_arcs_upper_bound(m, var.sub_arcs(), network, 0)
        self.__number_of_sink_arcs = sum(1 for _ in network.sink_arcs())

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield self.__cst

    def define_circular(self) -> None:
        """Define the bin as circular."""
        self.__cst.RHS = 1

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        self.__cst.RHS = self.__number_of_sink_arcs

    def deactivate(self) -> None:
        """Deactivate the constraint."""
        self.__cst.RHS = 0


@final
class CycleBeforeOut(Circular, PartiallyCircular):
    """Cycle before out constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.cycle_before_out(
            m,
            var.sub_vertices(),
            var.sub_arcs(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = 1


@final
class CycleBeforeIn(Circular, PartiallyCircular):
    """Cycle before in constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.cycle_before_in(
            m,
            var.sub_vertices(),
            var.sub_arcs(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = 1


@final
class SourceVerticesInArcUB(Circular, PartiallyCircular):
    """Source vertices in arc upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts, self.__number_of_preds = (
            cmn_lp_cst.s_connected_orfrag_incoming_arcs_ub(m, var.sub_arcs(), network)
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for k, cst in enumerate(self.__csts):
            cst.RHS = 1 + 2 * self.__number_of_preds[k]

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for k, cst in enumerate(self.__csts):
            cst.RHS = 1 + self.__number_of_preds[k]


@final
class SinkVerticesOutArcUB(Circular, PartiallyCircular):
    """Sink vertices out arc upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts, self.__number_of_succs = (
            cmn_lp_cst.t_connected_orfrag_outgoing_arcs_ub(m, var.sub_arcs(), network)
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for k, cst in enumerate(self.__csts):
            cst.RHS = 1 + 2 * self.__number_of_succs[k]

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for k, cst in enumerate(self.__csts):
            cst.RHS = 1 + self.__number_of_succs[k]


@final
class STConnexionLB(Circular, PartiallyCircular):
    """ST connexion lower bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.same_orfrag_connects_s_and_t_lb(
            m,
            var.sub_arcs(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = -1


@final
class STConnexionUB(Circular, PartiallyCircular):
    """ST connexion upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.same_orfrag_connects_s_and_t_ub(
            m,
            var.sub_arcs(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = 1


@final
class NonSeedSourceArcsUB(CircularWithSeed):
    """Non seed source arcs upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.non_seed_source_arc_ub(m, var.sub_arcs(), network, 1)

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular_with_seed(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_other(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 1


@final
class NonSeedSinkArcsUB(CircularWithSeed):
    """Non seed sink arcs upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.non_seed_sink_arc_ub(m, var.sub_arcs(), network, 1)

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular_with_seed(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_other(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 1


@final
class ActiveSourceArcFlowLB(Circular, PartiallyCircular):
    """Source arc flow lower bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
        min_active_flow: float,
    ) -> None:
        self.__csts = cmn_lp_cst.active_source_arc_has_strict_positive_flow(
            m,
            var.flows(),
            var.sub_arcs(),
            network,
            min_active_flow,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 1

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = 0


@final
class SourceArcFlowUB(Circular, PartiallyCircular):
    """Source arc flow upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.source_arc_flow_upper_bound(
            m,
            var.flows(),
            network,
            0,
        )
        self.__source_arc_flow_ub = sum(network.cap_s(a) for a in network.source_arcs())

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = self.__source_arc_flow_ub


@final
class ActiveSinkArcFlowLB(Circular, PartiallyCircular):
    """Sink arc flow lower bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
        min_active_flow: float,
    ) -> None:
        self.__csts = cmn_lp_cst.active_sink_arc_has_strict_positive_flow(
            m,
            var.flows(),
            var.sub_arcs(),
            network,
            min_active_flow,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 1

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = 0


@final
class SinkArcFlowUB(Circular, PartiallyCircular):
    """Sink arc flow upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = cmn_lp_cst.sink_arc_flow_upper_bound(
            m,
            var.flows(),
            network,
            0,
        )
        self.__sink_arc_flow_ub = sum(network.cap_t(a) for a in network.sink_arcs())

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = self.__sink_arc_flow_ub


@final
class BetaRevUB(Circular, PartiallyCircular):
    """Beta_rev upper bound constraints (c + pc)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = ccomp_cst.beta_rev_upper_bound(m, var.tree_edges(), network)
        self.__number_of_vertices = network.number_of_vertices()

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__csts:
            cst.RHS = self.__number_of_vertices


@final
class SourceSubtreeSize(Activate, Deactivate):
    """Source subtree size constraints."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = ccomp_cst.subtree_size_from_source(
            m,
            var.sub_vertices(),
            var.tree_edges().dtree(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def activate(self) -> None:
        """Define the bin as active."""
        for cst in self.__csts:
            cst.RHS = -1

    def deactivate(self) -> None:
        """Define the bin as inactive."""
        for cst in self.__csts:
            cst.RHS = 0


@final
class SinkSubtreeSizeRevVersion(Activate, Deactivate):
    """Sink subtree size constraints."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = ccomp_cst.subtree_size_for_sink_rev_version(
            m,
            var.tree_edges(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def activate(self) -> None:
        """Define the bin as active."""
        for cst in self.__csts:
            cst.RHS = 1

    def deactivate(self) -> None:
        """Define the bin as inactive."""
        for cst in self.__csts:
            cst.RHS = 0


@final
class NonSeedRootUB(CircularWithSeed):
    """Non seed root upper bound constraints (c & s)."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__csts = ccomp_cst.non_seed_root_upper_bound(
            m,
            var.root(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield from self.__csts

    def define_circular_with_seed(self) -> None:
        """Define the bin as circular with seed."""
        for cst in self.__csts:
            cst.RHS = 0

    def define_other(self) -> None:
        """Define the bin as circular free of seed or partially circular."""
        for cst in self.__csts:
            cst.RHS = 1


@final
class NumberOfRootVertices(Activate, Deactivate):
    """Number of root vertices constraints."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
    ) -> None:
        self.__cst = ccomp_cst.only_one_oriented_fragment_connects_the_source(
            m,
            var.root(),
            network,
        )

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield self.__cst

    def activate(self) -> None:
        """Define the bin as active."""
        self.__cst.RHS = 1

    def deactivate(self) -> None:
        """Define the bin as inactive."""
        self.__cst.RHS = 0


@final
class CumulativeLengthLB(Activate, Deactivate):
    """Cumulative length lower bound constraints."""

    def __init__(
        self,
        m: gp.Model,
        var: lp_vars.BinVariables,
        network: net.Network,
        min_cumulative_len: int,
    ) -> None:
        self.__cst = cmn_lp_cst.minimum_cumulative_length(
            m,
            var.sub_frag(),
            network,
            min_cumulative_len,
        )
        self.__min_cumulative_len = min_cumulative_len

    def constraints(self) -> Iterator[gp.Constr]:
        """Get the constraint."""
        yield self.__cst

    def activate(self) -> None:
        """Define the bin as active."""
        self.__cst.RHS = self.__min_cumulative_len

    def deactivate(self) -> None:
        """Define the bin as inactive."""
        self.__cst.RHS = 0


# REFACTOR use State: Inactive, otherwise Active with properties
class BinStateConstraints:
    """Bin state constraints.

    Constraints to change in order to activate the bin
    and to set it circular or partially circular
    """

    def __init__(self, constraints: Iterable[StructuralCst]) -> None:
        self.__constraints = list(constraints)

        self.__c_csts = [cst for cst in self.__constraints if isinstance(cst, Circular)]
        self.__pc_csts = [
            cst for cst in self.__constraints if isinstance(cst, PartiallyCircular)
        ]
        self.__s_csts = [cst for cst in self.__constraints if isinstance(cst, Seed)]
        self.__fs_csts = [
            cst for cst in self.__constraints if isinstance(cst, FreeOfSeed)
        ]
        self.__activate_csts = [
            cst for cst in self.__constraints if isinstance(cst, Activate)
        ]
        self.__deactivate_csts = [
            cst for cst in self.__constraints if isinstance(cst, Deactivate)
        ]
        self.__c_s_csts = [
            cst for cst in self.__constraints if isinstance(cst, CircularWithSeed)
        ]

        self.__topology = hmf_bins.Topology.CIRCULAR
        self.define_topology(self.__topology)

        self.__seed_constraint = hmf_bins.SeedConstraint.REQUIRED
        self.define_seed_constraint(self.__seed_constraint)

        self.__activate = False
        self.deactivate_bin()

    def topology(self) -> hmf_bins.Topology:
        """Get the topology of the bin."""
        return self.__topology

    def seed_constraint(self) -> hmf_bins.SeedConstraint:
        """Get the seed constraint."""
        return self.__seed_constraint

    def activate_bin(
        self,
        topology: hmf_bins.Topology,
        seed_constraint: hmf_bins.SeedConstraint,
    ) -> None:
        """Activate the bin."""
        for cst in self.__activate_csts:
            cst.activate()

        self.define_topology(topology)
        self.define_seed_constraint(seed_constraint)

        if (topology, seed_constraint) == (
            hmf_bins.Topology.CIRCULAR,
            hmf_bins.SeedConstraint.REQUIRED,
        ):
            for c_s_cst in self.__c_s_csts:
                c_s_cst.define_circular_with_seed()
        else:
            for c_s_cst in self.__c_s_csts:
                c_s_cst.define_other()

        self.__activate = True

    def deactivate_bin(self) -> None:
        """Deactivate the bin."""
        for cst in self.__deactivate_csts:
            cst.deactivate()
        self.__activate = False

    def define_circular(self) -> None:
        """Define the bin as circular."""
        for cst in self.__c_csts:
            cst.define_circular()
        self.__topology = hmf_bins.Topology.CIRCULAR

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        for cst in self.__pc_csts:
            cst.define_partially_circular()
        self.__topology = hmf_bins.Topology.PARTIALLY_CIRCULAR

    def define_topology(self, topology: hmf_bins.Topology) -> None:
        """Define the topology of the bin."""
        match topology:
            case hmf_bins.Topology.CIRCULAR:
                self.define_circular()
            case hmf_bins.Topology.PARTIALLY_CIRCULAR:
                self.define_partially_circular()

    def define_must_have_a_seed(self) -> None:
        """Constrain the bin to have a seed."""
        for cst in self.__s_csts:
            cst.define_must_have_a_seed()
        self.__seed_constraint = hmf_bins.SeedConstraint.REQUIRED

    def define_can_be_free_of_seed(self) -> None:
        """Relax the must-have-a-seed constraint."""
        for cst in self.__fs_csts:
            cst.define_can_be_free_of_seed()
        self.__seed_constraint = hmf_bins.SeedConstraint.NOT_REQUIRED

    def define_seed_constraint(self, seed_constraint: hmf_bins.SeedConstraint) -> None:
        """Set the seed constraint."""
        match seed_constraint:
            case hmf_bins.SeedConstraint.REQUIRED:
                self.define_must_have_a_seed()
            case hmf_bins.SeedConstraint.NOT_REQUIRED:
                self.define_can_be_free_of_seed()

    def is_activate(self) -> bool:
        """Check if the bin is activated."""
        return self.__activate
