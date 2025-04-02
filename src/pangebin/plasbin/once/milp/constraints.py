"""PangeBin-flow once MILP constraints."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.constraints as lp_cst
import pangebin.plasbin.network as net
import pangebin.plasbin.once.milp.variables as lp_vars


# ------------------------------------------------------------------------------------ #
#                                         MGCLB                                        #
# ------------------------------------------------------------------------------------ #
def set_mgclb_constraints(  # noqa: PLR0913
    m: gp.Model,
    var: lp_vars.MaxGCLabelBinScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    min_flow: float,
    min_cumulative_len: int,
    circular: bool,  # noqa: FBT001
) -> None:
    """Set MGCLB constraints."""
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
    # GC
    #
    lp_cst.exactly_one_interval_is_active(m, var.gc(), intervals)
    lp_cst.define_inflow_gc(
        m,
        var.frag(),
        var.gc(),
        var.flow(),
        var.inflow_gc(),
        network,
        intervals,
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
