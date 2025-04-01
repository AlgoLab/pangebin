"""PangeBin-flow MILP constraints library.

The following constraint function serve as base bricks to compose MILP models.
"""

from __future__ import annotations

from itertools import chain, product

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net


def active_fragments_active_one_of_their_orientations(
    m: gurobipy.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_f_vars: pb_lp_var.SubFragments,
    network: net.Network,
) -> None:
    """Active fragment equivalent to at least one active extremity constraint."""
    # DOCU MCF: + active fragment equivalent to at least one active extremity
    m.addConstrs(
        (
            sub_v_vars.x(frag) <= sub_f_vars.x(frag_id)
            for frag_id in network.fragment_ids()
            for frag in network.to_oriented(frag_id)
        ),
        name="active_extremity_implies_active_fragment",
    )
    m.addConstrs(
        (
            sub_f_vars.x(frag_id)
            <= gurobipy.quicksum(
                sub_v_vars.x(frag) for frag in network.to_oriented(frag_id)
            )
            for frag_id in network.fragment_ids()
        ),
        name="active_fragment_implies_active_extremity",
    )


def exactly_one_active_source_arc(
    m: gurobipy.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> None:
    """Exactly one out source arc constraint."""
    m.addConstr(
        gurobipy.quicksum(sub_arc_vars.s(a) for a in network.source_arcs()) == 1,
        name="exactly_one_out_source_link",
    )


def circularity(
    m: gurobipy.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> None:
    """Circularity constraint."""
    # DOCU MCF: + circularity
    source = network.SOURCE_VERTEX
    sink = network.SINK_VERTEX
    m.addConstrs(
        (
            sub_arc_vars.s((source, frag)) == sub_arc_vars.t((frag, sink))
            for seed in network.seeds()
            for frag in network.to_oriented(seed)
        ),
        name="same_source_and_sink",
    )


def arc_capacities_limit_arc_flows(
    m: gurobipy.Model,
    flow_vars: pb_lp_var.Flow,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> None:
    """Arc capacities limit arc flows."""
    m.addConstrs(
        (
            const
            for const in chain(
                (
                    flow_vars.s(a) <= network.cap_s(a) * sub_arc_vars.s(a)
                    for a in network.source_arcs()
                ),
                (
                    flow_vars.l(a) <= network.cap(a) * sub_arc_vars.l(a)
                    for a in network.link_arcs()
                ),
                (
                    flow_vars.t(a) <= network.cap_t(a) * sub_arc_vars.t(a)
                    for a in network.sink_arcs()
                ),
            )
        ),
        name="arc_capacities_limit_arc_flows",
    )


def fragment_coverages_limit_cumulative_flows(
    m: gurobipy.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
) -> None:
    """Fragment coverages limit cumulative flows."""
    m.addConstrs(
        (
            flow_vars.incoming_forward_reverse(network, frag_id)
            <= network.coverage(frag_id)
            for frag_id in network.gfa_graph().segment_names
        ),
        name="fragment_coverages_limit_cumulative_flows",
    )


def flow_conservation(
    m: gurobipy.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
) -> None:
    """Flow conservation."""
    m.addConstrs(
        (
            flow_vars.incoming(network, fragment)
            == flow_vars.outgoing(network, fragment)
            for fragment in network.oriented_fragments()
        ),
        name="flow_conservation",
    )


def total_flow_value(
    m: gurobipy.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
) -> None:
    """Set the total flow equals the source outgoing flow and the sink incoming flow."""
    m.addConstr(
        flow_vars.total()
        == gurobipy.quicksum(flow_vars.s(s_w) for s_w in network.source_arcs()),
        name="total_flow_equals_source_outgoing_flow",
    )
    m.addConstr(
        flow_vars.total()
        == gurobipy.quicksum(flow_vars.t(u_t) for u_t in network.sink_arcs()),
        name="total_flow_equals_sink_incoming_flow",
    )


def active_arcs_implies_active_fragments(
    m: gurobipy.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> None:
    """Active arcs implies active fragments."""
    for source_arc in network.source_arcs():
        m.addConstr(
            sub_arc_vars.s(source_arc) <= sub_v_vars.x(source_arc[1]),
            name=f"active_arc_{source_arc[0]}_{source_arc[1]}_implies_active_fragment_{source_arc[1]}",
        )
    for arc in network.link_arcs():
        m.addConstr(
            sub_arc_vars.l(arc) <= sub_v_vars.x(arc.predecessor()),
            name=f"active_arc_{arc}_implies_active_fragment_{arc.predecessor()}",
        )
        m.addConstr(
            sub_arc_vars.l(arc) <= sub_v_vars.x(arc.successor()),
            name=f"active_arc_{arc}_implies_active_fragment_{arc.successor()}",
        )
    for sink_arc in network.sink_arcs():
        m.addConstr(
            sub_arc_vars.t(sink_arc) <= sub_v_vars.x(sink_arc[0]),
            name=f"active_arc_{sink_arc[0]}_{sink_arc[1]}_implies_active_fragment_{sink_arc[0]}",
        )


# REFACTOR Potentially useless constraint
def active_fragments_imply_at_least_one_active_arc(
    m: gurobipy.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> None:
    """Active fragments imply at least one active arc."""
    m.addConstrs(
        (
            sub_v_vars.x(frag) <= sub_arc_vars.incoming(network, frag)
            for frag in network.oriented_fragments()
        ),
        name="active_fragments_imply_at_least_one_active_arc",
    )


def active_arcs_imply_vertices_in_same_component(
    m: gurobipy.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> None:
    """Active arcs imply vertices in the same connected component."""
    # REFACTOR potentially redondant with active_fragments_imply_at_least_one_active_arc
    for source_arc in network.source_arcs():
        m.addConstr(
            1 - sub_arc_vars.s(source_arc) >= 1 - sub_v_vars.x(source_arc[1]),
            name=f"active_arcs_imply_vertices_in_same_component_{source_arc[0]}_{source_arc[1]}",
        )
    for arc in network.link_arcs():
        m.addConstr(
            1 - sub_arc_vars.l(arc)
            >= sub_v_vars.x(arc.predecessor()) - sub_v_vars.x(arc.successor()),
            name=f"active_arcs_imply_vertices_in_same_component_{arc.predecessor()}_{arc.successor()}",
        )
    for sink_arc in network.sink_arcs():
        m.addConstr(
            1 - sub_arc_vars.t(sink_arc) >= sub_v_vars.x(sink_arc[0]) - 1,
            name=f"active_arcs_imply_vertices_in_same_component_{sink_arc[0]}_{sink_arc[1]}",
        )


def subtree_depth_min_distance(
    m: gurobipy.Model,
    ccomp_vars: pb_lp_var.ConnectedComponent,
    network: net.Network,
) -> None:
    """Subtree depth min distance."""
    #
    # For the source
    #
    m.addConstr(
        ccomp_vars.alpha()
        + gurobipy.quicksum(
            ccomp_vars.beta_s(source_arcs) for source_arcs in network.source_arcs()
        )
        <= 1,
        name="source_subtree_depth_min_distance",
    )
    #
    # For the fragments
    #
    m.addConstrs(
        (
            ccomp_vars.outgoing_beta(network, fragment)
            - ccomp_vars.incoming_beta(network, fragment)
            <= 1
            for fragment in network.oriented_fragments()
        ),
        name="subtree_depth_min_distance",
    )
    #
    # For the sink
    #
    m.addConstr(
        -gurobipy.quicksum(
            ccomp_vars.beta_t(sink_arcs) for sink_arcs in network.sink_arcs()
        )
        <= 1,
        name="sink_subtree_depth_min_distance",
    )


def arcs_in_tree_are_active(
    m: gurobipy.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    ccomp_vars: pb_lp_var.ConnectedComponent,
    network: net.Network,
) -> None:
    """Arcs in tree are active."""
    #
    # Source-arcs
    #
    m.addConstrs(
        (
            ccomp_vars.beta_s(source_arc)
            >= -sub_arc_vars.s(source_arc) * network.number_of_vertices()
            for source_arc in network.source_arcs()
        ),
        name="source_arcs_in_tree_are_active",
    )
    #
    # Link-arcs
    #
    m.addConstrs(
        (
            ccomp_vars.beta(arc) >= -sub_arc_vars.l(arc) * network.number_of_vertices()
            for arc in network.link_arcs()
        ),
        name="Link_arcs_in_tree_are_active",
    )
    #
    # Sink-arcs
    #
    m.addConstrs(
        (
            ccomp_vars.beta_t(sink_arc)
            >= -sub_arc_vars.t(sink_arc) * network.number_of_vertices()
            for sink_arc in network.sink_arcs()
        ),
        name="sink_arcs_in_tree_are_active",
    )


def alpha_is_the_number_of_vertices_in_the_solution(
    m: gurobipy.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    ccomp_vars: pb_lp_var.ConnectedComponent,
    network: net.Network,
) -> None:
    """Alpha is the number of vertices in the solution connected subgraph."""
    m.addConstr(
        ccomp_vars.alpha()
        == 2
        + gurobipy.quicksum(
            sub_v_vars.x(fragment) for fragment in network.oriented_fragments()
        ),
        name="alpha_is_the_number_of_vertices_in_the_solution",
    )


def active_arcs_have_flow_at_least_total_flow(
    m: gurobipy.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    flow_vars: pb_lp_var.Flow,
    pos_flow_vars: pb_lp_var.PositiveFlow,
    network: net.Network,
) -> None:
    """Active links have flow at least total flow."""
    max_seed_cov = max(
        network.coverage(frag_id) for frag_id in network.source_connected_fragment_ids()
    )
    m.addConstrs(
        (
            const
            for const in chain(
                (
                    pos_flow_vars.s(source_arc)
                    >= flow_vars.total()
                    - (1 - sub_arc_vars.s(source_arc)) * max_seed_cov
                    for source_arc in network.source_arcs()
                ),
                (
                    pos_flow_vars.l(arc)
                    >= flow_vars.total() - (1 - sub_arc_vars.l(arc)) * max_seed_cov
                    for arc in network.link_arcs()
                ),
                (
                    pos_flow_vars.t(sink_arc)
                    >= flow_vars.total() - (1 - sub_arc_vars.t(sink_arc)) * max_seed_cov
                    for sink_arc in network.sink_arcs()
                ),
            )
        ),
        name="active_links_have_flow_at_least_total_flow_1",
    )
    m.addConstrs(
        (
            const
            for const in chain(
                (
                    pos_flow_vars.s(source_arc) <= flow_vars.total()
                    for source_arc in network.source_arcs()
                ),
                (
                    pos_flow_vars.l(arc) <= flow_vars.total()
                    for arc in network.link_arcs()
                ),
                (
                    pos_flow_vars.t(sink_arc) <= flow_vars.total()
                    for sink_arc in network.sink_arcs()
                ),
            )
        ),
        name="active_links_have_flow_at_least_total_flow_2",
    )
    m.addConstrs(
        (
            const
            for const in chain(
                (
                    pos_flow_vars.s(source_arc) <= flow_vars.s(source_arc)
                    for source_arc in network.source_arcs()
                ),
                (
                    pos_flow_vars.l(arc) <= flow_vars.l(arc)
                    for arc in network.link_arcs()
                ),
                (
                    pos_flow_vars.t(sink_arc) <= flow_vars.t(sink_arc)
                    for sink_arc in network.sink_arcs()
                ),
            )
        ),
        name="active_links_have_flow_at_least_total_flow_3",
    )


def total_flow_is_strictly_positive(
    m: gurobipy.Model,
    flow_vars: pb_lp_var.Flow,
    epsilon_total_flow: float,
) -> None:
    """Total flow is strictly positive."""
    m.addConstr(
        flow_vars.total() >= epsilon_total_flow,
        name="total_flow_is_strictly_positive",
    )


def minimum_cumulative_length(
    m: gurobipy.Model,
    frag_vars: pb_lp_var.SubFragments,
    network: net.Network,
    minimum_cumulative_length: int,
) -> None:
    """Minimum cumulative length."""
    # DOCU note that the repetition is not considered
    m.addConstr(
        minimum_cumulative_length
        <= gurobipy.quicksum(
            frag_vars.x(frag_id)
            * gfa_segment.length(network.gfa_graph().segment(frag_id))
            for frag_id in network.fragment_ids()
        ),
        name="minimum_cumulative_length",
    )


# ------------------------------------------------------------------------------------ #
#                                          GC                                          #
# ------------------------------------------------------------------------------------ #
def exactly_one_interval_is_active(
    m: gurobipy.Model,
    gc_vars: pb_lp_var.GCIntervals,
    intervals: gc_items.Intervals,
) -> None:
    """Exactly one interval is active."""
    m.addConstr(
        gurobipy.quicksum(gc_vars.gc(interval) for interval in intervals) == 1,
        name="exactly_one_interval_is_active",
    )


def define_frag_gc(
    m: gurobipy.Model,
    frag_vars: pb_lp_var.SubFragments,
    gc_vars: pb_lp_var.GCIntervals,
    network: net.Network,
    intervals: gc_items.Intervals,
) -> None:
    """Define frag gc variables."""
    m.addConstrs(
        (
            gc_vars.frag_gc(frag_id, interval) <= frag_vars.x(frag_id)
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_1",
    )
    m.addConstrs(
        (
            gc_vars.frag_gc(frag_id, interval) <= gc_vars.gc(interval)
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_2",
    )
    m.addConstrs(
        (
            gc_vars.frag_gc(frag_id, interval)
            >= frag_vars.x(frag_id) + gc_vars.gc(interval) - 1
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_3",
    )
