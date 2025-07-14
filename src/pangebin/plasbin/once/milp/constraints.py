"""PangeBin-flow once MILP constraints."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.connected_component.constraints as ccomp_cst
import pangebin.plasbin.milp.constraints as cmn_lp_cst
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
) -> list[gp.Constr]:
    """Set MGCLB constraints."""
    constraints: list[gp.Constr] = []
    constraints += cmn_lp_cst.active_fragments_active_one_of_their_orientations(
        m,
        var.sub_vertices(),
        var.sub_frags(),
        network,
    )
    # XXX PAER TESTS
    # DOCU needed if source arc involves all vertices (and not only the seeds)
    constraints.append(
        cmn_lp_cst.active_seeds_lower_bound(m, var.sub_frags(), network, 1),
    )
    # DOCU UB for number of source and sink arcs
    # DOCU needed if source arc involves all vertices (and not only the seeds)
    constraints.append(
        cmn_lp_cst.source_arcs_upper_bound(m, var.sub_arcs(), network, 1),
    )
    constraints.append(cmn_lp_cst.sink_arcs_upper_bound(m, var.sub_arcs(), network, 1))

    # DOCU needed if source arc involves all vertices (and not only the seeds)
    # DOCU if source_arc, no other incoming link
    constraints += cmn_lp_cst.s_connected_orfrag_incoming_arcs_ub(
        m,
        var.sub_arcs(),
        network,
    )[0]
    # DOCU needed if source arc involves all vertices (and not only the seeds)
    # DOCU if sink_arc, no other outgoing link
    constraints += cmn_lp_cst.t_connected_orfrag_outgoing_arcs_ub(
        m,
        var.sub_arcs(),
        network,
    )[0]

    constraints += cmn_lp_cst.arc_capacities_limit_arc_flows(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
    )

    # DOCU modif vs decomp and binlab
    constraints += cmn_lp_cst.active_arcs_imply_strict_positive_flow(
        m,
        var.flows(),
        var.sub_arcs(),
        network,
        min_flow,
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
    constraints += cmn_lp_cst.active_arcs_imply_vertices_in_same_component(
        m,
        var.sub_vertices(),
        var.sub_arcs(),
        network,
    )
    #
    # Connectivity
    #
    # DOCU use of a beta rev

    constraints += ccomp_cst.arcs_in_tree_are_active(
        m,
        var.sub_arcs(),
        var.tree_edges().dtree(),
        network,
    )

    # DOCU needed if can consider circular bin
    # DOCU C: (opti) force each rev beta to be 0
    beta_rev_up_csts = ccomp_cst.beta_rev_upper_bound(m, var.tree_edges(), network)
    constraints += beta_rev_up_csts
    if not circular:
        for beta_rev_cst in beta_rev_up_csts:
            beta_rev_cst.RHS = 1

    constraints += ccomp_cst.rev_link_arcs_in_tree_are_active(
        m,
        var.sub_arcs(),
        var.tree_edges(),
        network,
    )

    constraints += ccomp_cst.subtree_size_from_source(
        m,
        var.sub_vertices(),
        var.tree_edges().dtree(),
        network,
    )
    constraints += ccomp_cst.subtree_size_for_sink_rev_version(
        m,
        var.tree_edges(),
        network,
    )

    constraints += ccomp_cst.subtree_size_fragment_rev_version(
        m,
        var.sub_vertices(),
        var.tree_edges(),
        network,
    )

    # DOCU needed if several source arcs
    constraints.append(
        ccomp_cst.only_one_oriented_fragment_connects_the_source(
            m,
            var.root(),
            network,
        ),
    )
    # DOCU needed if several source arcs
    constraints += ccomp_cst.at_most_one_source_arc_in_tree(
        m,
        var.tree_edges(),
        var.root(),
        network,
    )
    # DOCU needed if several source arcs
    constraints += ccomp_cst.root_priority(m, var.root(), var.sub_arcs(), network)
    #
    # GC
    #
    constraints.append(
        cmn_lp_cst.exactly_one_interval_is_active(m, var.gc(), intervals),
    )
    constraints += cmn_lp_cst.define_inflow_gc(
        m,
        var.sub_frags(),
        var.gc(),
        var.flows(),
        var.inflow_gc(),
        network,
        intervals,
    )
    #
    # Plasmid property
    #
    # DOCU MCF: + minimum flow value constraint
    constraints += cmn_lp_cst.total_flow_is_strictly_positive(m, var.flows(), min_flow)
    # DOCU MCF: + min cumulative len constraint
    constraints.append(
        cmn_lp_cst.minimum_cumulative_length(
            m,
            var.sub_frags(),
            network,
            min_cumulative_len,
        ),
    )
    # FIXME use other constraint
    if circular:
        # DOCU circularity remove a little flow to the first seed but no problems(?)
        # DOCU new consraints for circularity
        constraints += cmn_lp_cst.cycle_before_out(
            m,
            var.sub_vertices(),
            var.sub_arcs(),
            network,
        )
        constraints += cmn_lp_cst.cycle_before_in(
            m,
            var.sub_vertices(),
            var.sub_arcs(),
            network,
        )

        constraints += cmn_lp_cst.same_oriented_fragment_connects_s_and_t_eq(
            m,
            var.sub_arcs(),
            network,
        )
        # DOCU C: [opti] only the seeds can be connected to the source and the sink
        constraints += ccomp_cst.non_seed_root_upper_bound(m, var.root(), network)

    return constraints
