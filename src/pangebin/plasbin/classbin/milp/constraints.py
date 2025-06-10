"""PangeBin-flow classbin MILP constraints."""

from __future__ import annotations

import gurobipy as gp

import pangebin.plasbin.classbin.milp.variables as lp_vars
import pangebin.plasbin.milp.constraints as lp_cst
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net


# REFACTOR use config dataclass for each model
# REFACTOR move to classify MILP module
def set_classify_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.Classify,
    network: net.Network,
    min_flow: float,
    plasmidness_coefficient: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    min_cumulative_len: int,
) -> list[gp.Constr]:
    """Set class constraints."""
    constraints: list[gp.Constr] = []
    constraints += lp_cst.active_fragments_active_one_of_their_orientations(
        m,
        var.sub_vertices(),
        var.frag(),
        network,
    )

    constraints += lp_cst.arc_capacities_limit_arc_flows(
        m,
        var.flow(),
        var.sub_arcs(),
        network,
    )
    constraints += lp_cst.fragment_coverages_limit_cumulative_flows(
        m,
        var.flow(),
        network,
    )
    constraints += lp_cst.active_arcs_imply_strict_positive_flow(
        m,
        var.flow(),
        var.sub_arcs(),
        network,
        min_flow,
    )
    constraints += lp_cst.flow_conservation(m, var.flow(), network)
    constraints += lp_cst.total_flow_value(m, var.flow(), network)
    constraints += lp_cst.active_arcs_implies_active_fragments(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    constraints += lp_cst.active_fragments_imply_at_least_one_active_arc(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    constraints += lp_cst.active_arcs_imply_vertices_in_same_component(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    #
    # Plasmid property
    #
    # DOCU Plasmidness lower bound
    # BUG tmp
    # constraints += _plasmidness_lower_bound(
    #     m,
    #     var.flow(),
    #     network,
    #     obj_fun_domain,
    #     plasmidness_coefficient,
    # )
    # DOCU MCF: + minimum flow value constraint
    constraints += lp_cst.total_flow_is_strictly_positive(m, var.flow(), min_flow)
    # DOCU MCF: + min cumulatie len constraint
    constraints.append(
        lp_cst.minimum_cumulative_length(
            m,
            var.frag(),
            network,
            min_cumulative_len,
        ),
    )
    # DOCU new constraint classify: use all seeds
    constraints.append(
        m.addConstr(
            gp.quicksum(var.frag().frag(seed) for seed in network.seeds())
            == len(network.seeds()),
        ),
    )
    #
    # Connectivity
    #
    constraints += lp_cst.subtree_depth_min_distance(m, var.ccomp(), network)
    constraints += lp_cst.arcs_in_tree_are_active(
        m,
        var.sub_arcs(),
        var.ccomp(),
        network,
    )
    constraints += lp_cst.alpha_is_the_number_of_vertices_in_the_solution(
        m,
        var.sub_vertices(),
        var.ccomp(),
        network,
    )
    #
    # Seed tree constraints
    #
    constraints += lp_cst.arc_in_seed_tree_is_active(
        m,
        var.seed_tree_arcs(),
        var.sub_arcs(),
        network,
    )
    constraints += lp_cst.size_of_seed_tree_under_the_source(
        m,
        var.sub_vertices(),
        var.seed_tree_arcs(),
        network,
    )
    constraints += lp_cst.size_of_seed_subtree(
        m,
        var.sub_vertices(),
        var.seed_tree_arcs(),
        network,
    )
    constraints += lp_cst.active_source_arc_implies_subgraph_contains_seed(
        m,
        var.sub_arcs(),
        var.seed_tree_arcs(),
        network,
    )

    return constraints
