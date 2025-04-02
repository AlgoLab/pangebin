"""PangeBin-flow binlab MILP constraints."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.binlab.milp.objectives as lp_obj
import pangebin.plasbin.binlab.milp.variables as lp_vars
import pangebin.plasbin.milp.constraints as lp_cst
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net


# ------------------------------------------------------------------------------------ #
#                                          MBS                                         #
# ------------------------------------------------------------------------------------ #
def set_mbs_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.MaxBinScore,
    network: net.Network,
    min_flow: float,
    min_cumulative_len: int,
    circular: bool,  # noqa: FBT001
) -> None:
    """Set MBS constraints."""
    lp_cst.active_fragments_active_one_of_their_orientations(
        m,
        var.sub_v(),
        var.frag(),
        network,
    )
    lp_cst.exactly_one_active_source_arc(m, var.sub_arc(), network)
    lp_cst.arc_capacities_limit_arc_flows(m, var.flow(), var.sub_arc(), network)
    lp_cst.fragment_coverages_limit_cumulative_flows(m, var.flow(), network)
    lp_cst.flow_conservation(m, var.flow(), network)
    lp_cst.total_flow_value(m, var.flow(), network)
    lp_cst.active_arcs_implies_active_fragments(
        m,
        var.sub_v(),
        var.sub_arc(),
        network,
    )
    lp_cst.active_fragments_imply_at_least_one_active_arc(
        m,
        var.sub_v(),
        var.sub_arc(),
        network,
    )
    lp_cst.active_arcs_imply_vertices_in_same_component(
        m,
        var.sub_v(),
        var.sub_arc(),
        network,
    )
    #
    # Connectivity
    #
    lp_cst.subtree_depth_min_distance(m, var.ccomp(), network)
    lp_cst.arcs_in_tree_are_active(m, var.sub_arc(), var.ccomp(), network)
    lp_cst.alpha_is_the_number_of_vertices_in_the_solution(
        m,
        var.sub_v(),
        var.ccomp(),
        network,
    )
    #
    # Arc flow lower bound
    #
    lp_cst.active_arcs_have_flow_at_least_total_flow(
        m,
        var.sub_arc(),
        var.flow(),
        var.pos_flow(),
        network,
    )
    #
    # Plasmid property
    #
    # DOCU MCF: + minimum flow value constraint
    lp_cst.total_flow_is_strictly_positive(m, var.flow(), min_flow)
    # DOCU MCF: + min cumulatie len constraint
    lp_cst.minimum_cumulative_length(m, var.frag(), network, min_cumulative_len)
    if circular:
        # DOCU circularity remove a little flow to the first seed but no problems(?)
        lp_cst.circularity(m, var.sub_arc(), network)


# ------------------------------------------------------------------------------------ #
#                                          MLS                                         #
# ------------------------------------------------------------------------------------ #
def add_mls_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.MaxLabScore,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    network: net.Network,
    intervals: gc_items.Intervals,
    gamma_mbs: float,
    previous_binning_score: float,
) -> None:
    """Add MLS constraints."""
    _binning_score_lower_bound(
        m,
        var,
        obj_fun_domain,
        gamma_mbs,
        previous_binning_score,
        network,
    )
    lp_cst.exactly_one_interval_is_active(m, var.gc(), intervals)
    lp_cst.define_frag_gc(m, var.frag(), var.gc(), var.frag_gc(), network, intervals)


def _binning_score_lower_bound(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.MaxLabScore,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    coefficient: float,
    previous_binning_score: float,
    network: net.Network,
) -> None:
    """Binning score lower bound."""
    # DOCU MGC: be carefull when previous_binning_score is < 0
    m.addConstr(
        previous_binning_score - (1 - coefficient) * abs(previous_binning_score)
        <= lp_obj.binning_score(network, var.flow(), var.frag(), obj_fun_domain),
        name="binning_score_lower_bound",
    )


# ------------------------------------------------------------------------------------ #
#                                         MRBS                                         #
# ------------------------------------------------------------------------------------ #
def add_mrbs_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.MaxRefBinScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    mls_obj_value: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Add MRBS constraints."""
    _fix_gc_score(m, var, network, intervals, mls_obj_value, obj_fun_domain)


def _fix_gc_score(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.MaxRefBinScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    mls_obj_value: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    # FIXME numerical problem can occur
    # Warning: max constraint violation (9.5055e-06) exceeds tolerance
    m.addConstr(
        lp_obj.gc_score(network, intervals, var.frag_gc(), obj_fun_domain)
        >= mls_obj_value,
        name="fix_gc_score",
    )
