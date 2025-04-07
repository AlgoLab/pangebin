"""Hiearchical binning constraints."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.decomp.milp.objectives as lp_obj
import pangebin.plasbin.decomp.milp.variables as lp_var
import pangebin.plasbin.milp.constraints as lp_cst
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net


# ------------------------------------------------------------------------------------ #
#                                          MCF                                         #
# ------------------------------------------------------------------------------------ #
def set_mcf_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_var.MaxCovFlow,
    network: net.Network,
    min_flow: float,
    min_cumulative_len: int,
    circular: bool,  # noqa: FBT001
) -> None:
    """Set MCF constraints."""
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
#                                          MGC                                         #
# ------------------------------------------------------------------------------------ #
def add_mgc_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_var.MaxGC,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    network: net.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
    previous_coverage_score: float,
) -> None:
    """Add MGC constraints."""
    _coverage_score_lower_bound(
        m,
        var,
        obj_fun_domain,
        coefficient,
        previous_coverage_score,
        network,
    )
    lp_cst.exactly_one_interval_is_active(m, var.gc(), intervals)
    lp_cst.define_frag_gc(m, var.frag(), var.gc(), var.frag_gc(), network, intervals)


def _coverage_score_lower_bound(  # noqa: PLR0913
    m: gp.Model,
    var: lp_var.MaxGC,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    coefficient: float,
    previous_coverage_score: float,
    network: net.Network,
) -> None:
    """Coverage score lower bound."""
    # DOCU MGC: be carefull when previous_coverage_score is < 0
    m.addConstr(
        previous_coverage_score - (1 - coefficient) * abs(previous_coverage_score)
        <= lp_obj.coverage_score(network, var.flow(), var.frag(), obj_fun_domain),
        name="coverage_score_lower_bound",
    )


# ------------------------------------------------------------------------------------ #
#                                          MPS                                         #
# ------------------------------------------------------------------------------------ #
def add_mps_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_var.MaxPlasmidScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
    previous_gc_score: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Add MPS constraints to MGC model."""
    _gc_score_lower_bound(
        m,
        var,
        network,
        intervals,
        coefficient,
        previous_gc_score,
        obj_fun_domain,
    )


def _gc_score_lower_bound(  # noqa: PLR0913
    m: gp.Model,
    var: lp_var.MaxPlasmidScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
    previous_gc_score: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """GC score lower bound."""
    # DOCU MPS: be carefull when previous_gc_score is < 0
    m.addConstr(
        previous_gc_score - (1 - coefficient) * abs(previous_gc_score)
        <= lp_obj.gc_score(network, intervals, var.frag_gc(), obj_fun_domain),
        name="gc_score_lower_bound",
    )


# ------------------------------------------------------------------------------------ #
#                                         MRCF                                         #
# ------------------------------------------------------------------------------------ #
def add_mrcf_constraints(
    m: gp.Model,
    var: lp_var.MaxPlasmidScore,
    network: net.Network,
    mps_obj_value: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Add MRCF constraints."""
    _fix_plasmidness_score(m, var, network, mps_obj_value, obj_fun_domain)


def _fix_plasmidness_score(
    m: gp.Model,
    var: lp_var.MaxPlasmidScore,
    network: net.Network,
    mps_obj_value: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    # FIXME numerical problem can occur
    # Warning: max constraint violation (9.5055e-06) exceeds tolerance
    m.addConstr(
        lp_obj.plasmidness_score(network, var.frag(), obj_fun_domain) >= mps_obj_value,
        name="fix_plasmidness_score",
    )
