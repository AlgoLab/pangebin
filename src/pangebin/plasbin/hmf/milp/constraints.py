"""PangeBin-flow classbin MILP constraints."""

from __future__ import annotations

from collections.abc import Iterable
from typing import cast

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
    constraints += cmn_lp_cst.active_fragments_active_one_of_their_orientations(
        m,
        var.sub_vertices(),
        var.sub_frag(),
        network,
    )
    # DOCU change RHS weither bin is active or not
    lb_seeds_cst = cmn_lp_cst.active_seeds_lower_bound(m, var.sub_frag(), network, 0)
    constraints.append(lb_seeds_cst)

    #
    # Source-sink constraints
    #
    # DOCU when C, only one source arc
    source_arcs_ub_cst = cmn_lp_cst.source_arcs_upper_bound(
        m,
        var.sub_arcs(),
        network,
        1,
    )
    constraints.append(source_arcs_ub_cst)

    # DOCU when C, only one sink arc
    sink_arcs_ub_cst = cmn_lp_cst.sink_arcs_upper_bound(m, var.sub_arcs(), network, 1)
    constraints.append(sink_arcs_ub_cst)

    # DOCU C: RHS = 0, PC: RHS = 1
    cycle_before_out_csts = cmn_lp_cst.cycle_before_out(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    cycle_before_in_csts = cmn_lp_cst.cycle_before_in(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    constraints += cycle_before_out_csts
    constraints += cycle_before_in_csts

    # DOCU when PC, no other in/out link if source/sink connected
    # * PC: no other incoming link if source connected
    incoming_arcs_ub_for_s_connected_orfrag_csts, _nb_preds = (
        cmn_lp_cst.s_connected_orfrag_incoming_arcs_ub(m, var.sub_arcs(), network)
    )
    # * PC: no other outgoing link if sink connected
    outgoing_arcs_ub_for_t_connected_orfrag_csts, _nb_succs = (
        cmn_lp_cst.t_connected_orfrag_outgoing_arcs_ub(m, var.sub_arcs(), network)
    )
    constraints += incoming_arcs_ub_for_s_connected_orfrag_csts
    constraints += outgoing_arcs_ub_for_t_connected_orfrag_csts

    # DOCU C: the seed connected to s is also connected to t
    # Two inequalities to change
    same_seed_connects_s_and_t_lb_csts = cmn_lp_cst.same_seed_connects_s_and_t_lb(
        m,
        var.sub_arcs(),
        network,
    )
    same_seed_connects_s_and_t_ub_csts = cmn_lp_cst.same_seed_connects_s_and_t_ub(
        m,
        var.sub_arcs(),
        network,
    )
    constraints += same_seed_connects_s_and_t_lb_csts
    constraints += same_seed_connects_s_and_t_ub_csts

    # DOCU let flow be 0 on source and sink arc when circular bin
    # * C: the flow value on source and sink arcs equal is 0 AND y (must) can still be 1
    source_arc_flow_ub_csts = cmn_lp_cst.source_arc_flow_upper_bound(
        m,
        var.flows(),
        network,
        0,
    )
    sink_arc_flow_ub_csts = cmn_lp_cst.sink_arc_flow_upper_bound(
        m,
        var.flows(),
        network,
        0,
    )
    source_arc_flow_lb_csts = cmn_lp_cst.active_source_arc_has_strict_positive_flow(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
        config.min_flow(),
    )
    constraints += cmn_lp_cst.active_link_arc_has_strict_positive_flow(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
        config.min_flow(),
    )
    sink_arc_flow_lb_csts = cmn_lp_cst.active_sink_arc_has_strict_positive_flow(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
        config.min_flow(),
    )
    constraints += source_arc_flow_lb_csts
    constraints += sink_arc_flow_lb_csts
    constraints += source_arc_flow_ub_csts
    constraints += sink_arc_flow_ub_csts
    source_flow_ub = sum(network.cap_s(a) for a in network.source_arcs())
    sink_flow_ub = sum(network.cap_t(a) for a in network.sink_arcs())

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
    # constraints += lp_cst.active_arcs_imply_vertices_in_same_component(
    #     m,
    #     var.sub_vertices(),
    #     var.sub_arcs(),
    #     network,
    # )
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
    beta_rev_ub_csts = ccomp_cst.beta_rev_upper_bound(m, var.tree_edges(), network)
    constraints += beta_rev_ub_csts

    constraints += ccomp_cst.rev_link_arcs_in_tree_are_active(
        m,
        var.sub_arcs(),
        var.tree_edges(),
        network,
    )

    # DOCU change RHS weither bin is active or not
    subtree_size_from_source_csts = ccomp_cst.subtree_size_from_source(
        m,
        var.sub_vertices(),
        var.tree_edges().dtree(),
        network,
    )
    # DOCU change RHS weither bin is active or not
    subtree_size_for_sink_rev_version_csts = (
        ccomp_cst.subtree_size_for_sink_rev_version(
            m,
            var.tree_edges(),
            network,
        )
    )
    constraints += subtree_size_from_source_csts
    constraints += subtree_size_for_sink_rev_version_csts

    constraints += ccomp_cst.subtree_size_fragment_rev_version(
        m,
        var.sub_vertices(),
        var.tree_edges(),
        network,
    )

    # DOCU C: [opti] only the seeds can be connected to the source and the sink
    non_seed_root_ub_csts = ccomp_cst.non_seed_root_upper_bound(
        m,
        var.root(),
        network,
    )
    constraints += non_seed_root_ub_csts

    # DOCU change RHS weither bin is active or not
    nb_of_source_connected_orientfrag_cst = (
        ccomp_cst.only_one_oriented_fragment_connects_the_source(
            m,
            var.root(),
            network,
        )
    )
    constraints.append(nb_of_source_connected_orientfrag_cst)
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
    min_cumulative_len_cst = cmn_lp_cst.minimum_cumulative_length(
        m,
        var.sub_frag(),
        network,
        config.min_cumulative_len(),
    )
    constraints.append(min_cumulative_len_cst)

    m.update()  # necessary to change the rhs
    # REFACTOR use compositions (and perhaps sub compositions)
    return constraints, BinStateConstraints(
        lb_seeds_cst,
        source_arcs_ub_cst,
        sink_arcs_ub_cst,
        network.number_of_vertices(),
        cycle_before_out_csts,
        cycle_before_in_csts,
        incoming_arcs_ub_for_s_connected_orfrag_csts,
        outgoing_arcs_ub_for_t_connected_orfrag_csts,
        _nb_preds,
        _nb_succs,
        network.ub_number_of_edges(),
        same_seed_connects_s_and_t_lb_csts,
        same_seed_connects_s_and_t_ub_csts,
        subtree_size_from_source_csts,
        subtree_size_for_sink_rev_version_csts,
        beta_rev_ub_csts,
        non_seed_root_ub_csts,
        source_arc_flow_lb_csts,
        sink_arc_flow_lb_csts,
        source_arc_flow_ub_csts,
        sink_arc_flow_ub_csts,
        source_flow_ub,
        sink_flow_ub,
        nb_of_source_connected_orientfrag_cst,
        min_cumulative_len_cst,
        config.min_cumulative_len(),
    )


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
                for k in range(bin_interval[0], bin_interval[1])
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


class BinStateConstraints:
    """Bin state constraints.

    Constraints to change in order to activate the bin
    and to set it circular or partially circular
    """

    def __init__(
        self,
        lb_seeds_cst: gp.Constr,
        source_arcs_ub_cst: gp.Constr,
        sink_arcs_ub_cst: gp.Constr,
        number_of_vertices: int,
        cycle_before_out_csts: Iterable[gp.Constr],
        cycle_before_in_csts: Iterable[gp.Constr],
        incoming_arcs_ub_for_s_connected_orfrag_csts: Iterable[gp.Constr],
        outgoing_arcs_ub_for_t_connected_orfrag_csts: Iterable[gp.Constr],
        nb_pred: Iterable[int],
        nb_succ: Iterable[int],
        ub_number_of_edges: int,
        same_seed_connects_s_and_t_lb_csts: Iterable[gp.Constr],
        same_seed_connects_s_and_t_ub_csts: Iterable[gp.Constr],
        subtree_size_from_source_csts: Iterable[gp.Constr],
        subtree_size_for_sink_rev_version_csts: Iterable[gp.Constr],
        beta_rev_ub_csts: Iterable[gp.Constr],
        non_seed_root_ub_csts: Iterable[gp.Constr],
        source_arc_flow_lb_csts: Iterable[gp.Constr],
        sink_arc_flow_lb_csts: Iterable[gp.Constr],
        source_arc_flow_ub_csts: Iterable[gp.Constr],
        sink_arc_flow_ub_csts: Iterable[gp.Constr],
        source_flow_ub: float,
        sink_flow_ub: float,
        nb_of_source_connected_orientfrag_csts: gp.Constr,
        min_cumulative_len_cst: gp.Constr,
        min_cumulative_len: int,
    ) -> None:
        self.__lb_seeds_cst = lb_seeds_cst

        self.__source_arcs_ub_cst = source_arcs_ub_cst
        self.__sink_arcs_ub_cst = sink_arcs_ub_cst
        self.__number_of_vertices = number_of_vertices

        self.__cycle_before_out_csts = list(cycle_before_out_csts)
        self.__cycle_before_in_csts = list(cycle_before_in_csts)

        self.__incoming_arcs_ub_for_s_connected_orfrag_csts = list(
            incoming_arcs_ub_for_s_connected_orfrag_csts,
        )
        self.__outgoing_arcs_ub_for_t_connected_orfrag_csts = list(
            outgoing_arcs_ub_for_t_connected_orfrag_csts,
        )
        self.__nb_pred = list(nb_pred)
        self.__nb_succ = list(nb_succ)
        self.__ub_number_of_edges = ub_number_of_edges

        self.__same_seed_connects_s_and_t_lb_csts = list(
            same_seed_connects_s_and_t_lb_csts,
        )
        self.__same_seed_connects_s_and_t_ub_csts = list(
            same_seed_connects_s_and_t_ub_csts,
        )

        self.__subtree_size_from_source_csts = list(subtree_size_from_source_csts)
        self.__subtree_size_for_sink_rev_version_csts = list(
            subtree_size_for_sink_rev_version_csts,
        )
        self.__beta_rev_ub_csts = list(beta_rev_ub_csts)
        self.__non_seed_root_ub_csts = list(non_seed_root_ub_csts)

        self.__source_arc_flow_lb_csts = list(source_arc_flow_lb_csts)
        self.__sink_arc_flow_lb_csts = list(sink_arc_flow_lb_csts)
        self.__source_arc_flow_ub_csts = list(source_arc_flow_ub_csts)
        self.__sink_arc_flow_ub_csts = list(sink_arc_flow_ub_csts)
        self.__source_flow_ub = source_flow_ub
        self.__sink_flow_ub = sink_flow_ub

        self.__nb_of_source_connected_orientfrag_cst = (
            nb_of_source_connected_orientfrag_csts
        )

        self.__min_cumulative_len_cst = min_cumulative_len_cst
        self.__min_cumulative_len = min_cumulative_len

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

    def activate_bin(self) -> None:
        """Activate the bin."""
        self.__lb_seeds_cst.RHS = 1
        self.__source_arcs_ub_cst.RHS = self.__number_of_vertices
        self.__sink_arcs_ub_cst.RHS = self.__number_of_vertices
        for constraint in self.__subtree_size_from_source_csts:
            constraint.RHS = -1
        for constraint in self.__subtree_size_for_sink_rev_version_csts:
            constraint.RHS = 1
        self.__nb_of_source_connected_orientfrag_cst.RHS = 1
        self.__min_cumulative_len_cst.RHS = self.__min_cumulative_len
        self.__activate = True

    def deactivate_bin(self) -> None:
        """Deactivate the bin."""
        self.__lb_seeds_cst.RHS = 0
        self.__source_arcs_ub_cst.RHS = 0
        self.__sink_arcs_ub_cst.RHS = 0
        for constraint in self.__subtree_size_from_source_csts:
            constraint.RHS = 0
        for constraint in self.__subtree_size_for_sink_rev_version_csts:
            constraint.RHS = 0
        self.__nb_of_source_connected_orientfrag_cst.RHS = 0
        self.__min_cumulative_len_cst.RHS = 0
        self.__activate = False

    def define_circular(self) -> None:
        """Define the bin as circular."""
        self.__source_arcs_ub_cst.RHS = 1
        self.__sink_arcs_ub_cst.RHS = 1
        for constraint in self.__cycle_before_out_csts:
            constraint.RHS = 0
        for constraint in self.__cycle_before_in_csts:
            constraint.RHS = 0
        for k, constraint in enumerate(
            self.__incoming_arcs_ub_for_s_connected_orfrag_csts,
        ):
            constraint.RHS = 1 + 2 * self.__nb_pred[k]
        for k, constraint in enumerate(
            self.__outgoing_arcs_ub_for_t_connected_orfrag_csts,
        ):
            constraint.RHS = 1 + 2 * self.__nb_succ[k]
        for constraint in self.__same_seed_connects_s_and_t_lb_csts:
            constraint.RHS = 0
        for constraint in self.__same_seed_connects_s_and_t_ub_csts:
            constraint.RHS = 0
        for constraint in self.__beta_rev_ub_csts:
            constraint.RHS = 0
        for constraint in self.__non_seed_root_ub_csts:
            constraint.RHS = 0
        for constraint in self.__source_arc_flow_lb_csts:
            constraint.RHS = 1
        for constraint in self.__sink_arc_flow_lb_csts:
            constraint.RHS = 1
        for constraint in self.__source_arc_flow_ub_csts:
            constraint.RHS = 0
        for constraint in self.__sink_arc_flow_ub_csts:
            constraint.RHS = 0
        self.__topology = hmf_bins.Topology.CIRCULAR

    def define_partially_circular(self) -> None:
        """Define the bin as partially circular."""
        self.__source_arcs_ub_cst.RHS = self.__ub_number_of_edges
        self.__sink_arcs_ub_cst.RHS = self.__ub_number_of_edges
        for constraint in self.__cycle_before_out_csts:
            constraint.RHS = 1
        for constraint in self.__cycle_before_in_csts:
            constraint.RHS = 1
        # DOCU trick on parcirc cst for s and t connected repeated frag
        # In fact, if we begin with a repeat looping at the beginning,
        # we cannot loop over it because the beginning is connected to the source
        # However, we can out at the end of the loop with T
        for k, constraint in enumerate(
            self.__incoming_arcs_ub_for_s_connected_orfrag_csts,
        ):
            constraint.RHS = 1 + self.__nb_pred[k]
        for k, constraint in enumerate(
            self.__outgoing_arcs_ub_for_t_connected_orfrag_csts,
        ):
            constraint.RHS = 1 + self.__nb_succ[k]
        for constraint in self.__same_seed_connects_s_and_t_lb_csts:
            constraint.RHS = -1
        for constraint in self.__same_seed_connects_s_and_t_ub_csts:
            constraint.RHS = 1
        for constraint in self.__beta_rev_ub_csts:
            constraint.RHS = self.__number_of_vertices
        for constraint in self.__non_seed_root_ub_csts:
            constraint.RHS = 1
        for constraint in self.__source_arc_flow_lb_csts:
            constraint.RHS = 0
        for constraint in self.__sink_arc_flow_lb_csts:
            constraint.RHS = 0
        for constraint in self.__source_arc_flow_ub_csts:
            constraint.RHS = self.__source_flow_ub
        for constraint in self.__sink_arc_flow_ub_csts:
            constraint.RHS = self.__sink_flow_ub
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
        self.__lb_seeds_cst.RHS = 1
        self.__seed_constraint = hmf_bins.SeedConstraint.REQUIRED

    def define_can_be_free_of_seed(self) -> None:
        """Relax the must-have-a-seed constraint."""
        self.__lb_seeds_cst.RHS = 0
        if self.__topology == hmf_bins.Topology.CIRCULAR:
            # DOCU C-FS: A non-seed can connect the source in the tree
            for constraint in self.__non_seed_root_ub_csts:
                constraint.RHS = 1
            # OPTIMIZE C-FS: help the solver to logic tree root and only one source-arc
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
