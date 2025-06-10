"""PangeBin-flow MILP constraints library.

The following constraint function serve as base bricks to compose MILP models.
"""

from __future__ import annotations

from itertools import chain, product

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net


def active_fragments_active_one_of_their_orientations(
    m: gp.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_f_vars: pb_lp_var.SubFragments,
    network: net.Network,
) -> list[gp.Constr]:
    """Active fragment equivalent to at least one active extremity constraint."""
    constraints: list[gp.Constr] = []
    # DOCU MCF: + active fragment equivalent to at least one active extremity
    constraints.extend(
        m.addConstrs(
            (
                sub_v_vars.x(frag) <= sub_f_vars.frag(frag_id)
                for frag_id in network.fragment_ids()
                for frag in network.to_oriented(frag_id)
            ),
            name="active_extremity_implies_active_fragment",
        ).values(),
    )
    constraints.extend(
        m.addConstrs(
            (
                sub_f_vars.frag(frag_id)
                <= gp.quicksum(
                    sub_v_vars.x(frag) for frag in network.to_oriented(frag_id)
                )
                for frag_id in network.fragment_ids()
            ),
            name="active_fragment_implies_active_extremity",
        ).values(),
    )
    return constraints


def exactly_one_active_source_arc(
    m: gp.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> gp.Constr:
    """Exactly one out source arc constraint."""
    return m.addConstr(
        gp.quicksum(sub_arc_vars.s(a) for a in network.source_arcs()) == 1,
        name="exactly_one_out_source_link",
    )


def circularity(
    m: gp.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Circularity constraint."""
    constraints: list[gp.Constr] = []
    # DOCU MCF: + circularity
    source = network.SOURCE_VERTEX
    sink = network.SINK_VERTEX
    constraints.extend(
        m.addConstrs(
            (
                sub_arc_vars.s((source, v)) == sub_arc_vars.t((v, sink))
                for frag_id in network.source_connected_fragment_ids()
                for v in network.to_oriented(frag_id)
            ),
            name="same_source_and_sink",
        ).values(),
    )
    return constraints


def arc_capacities_limit_arc_flows(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Arc capacities limit arc flows."""
    # DOCU change interpretation
    constraints: list[gp.Constr] = []
    constraints.extend(
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
        ).values(),
    )
    return constraints


# ------------------------------------------------------------------------------------ #
#                          Active Arc Has Strict Positive Flow                         #
# ------------------------------------------------------------------------------------ #
def active_arcs_imply_strict_positive_flow(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
    flow_lower_bound: float,
) -> list[gp.Constr]:
    """Active arcs imply strict positive flow."""
    # DOCU new required constraint for all the models(? see rmk below)
    # FIXME add to models binlab and decomp?
    # See if Pos flow csts was not sufficient before
    constraints: list[gp.Constr] = []
    invert_epsilon_f = 1 / flow_lower_bound
    constraints.extend(
        m.addConstrs(
            (
                const
                for const in chain(
                    (
                        sub_arc_vars.s(a) <= invert_epsilon_f * flow_vars.s(a)
                        for a in network.source_arcs()
                    ),
                    (
                        sub_arc_vars.l(a) <= invert_epsilon_f * flow_vars.l(a)
                        for a in network.link_arcs()
                    ),
                    (
                        sub_arc_vars.t(a) <= invert_epsilon_f * flow_vars.t(a)
                        for a in network.sink_arcs()
                    ),
                )
            ),
            name="active_arcs_imply_strict_positive_flow",
        ).values(),
    )
    return constraints


def active_source_arc_has_strict_positive_flow(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
    flow_lower_bound: float,
) -> list[gp.Constr]:
    """Active source arc has a strict positive flow.

    Note
    ----
    To authorize source arc to be active with a null flow,
    change the RHS to `1` (default: `0`)
    """
    invert_epsilon_f = 1 / flow_lower_bound
    return list(
        m.addConstrs(
            (
                sub_arc_vars.s(a) - invert_epsilon_f * flow_vars.s(a) <= 0
                for a in network.source_arcs()
            ),
            name="active_source_arc_has_strict_positive_flow",
        ).values(),
    )


def active_link_arc_has_strict_positive_flow(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
    flow_lower_bound: float,
) -> list[gp.Constr]:
    """Active link arc has a strict positive flow.

    Note
    ----
    To authorize link arc to be active with a null flow,
    change the RHS to `1` (default: `0`)
    """
    invert_epsilon_f = 1 / flow_lower_bound
    return list(
        m.addConstrs(
            (
                sub_arc_vars.l(a) - invert_epsilon_f * flow_vars.l(a) <= 0
                for a in network.link_arcs()
            ),
            name="active_link_arc_has_strict_positive_flow",
        ).values(),
    )


def active_sink_arc_has_strict_positive_flow(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
    flow_lower_bound: float,
) -> list[gp.Constr]:
    """Active sink arc has a strict positive flow.

    Note
    ----
    To authorize sink arc to be active with a null flow,
    change the RHS to `1` (default: `0`)
    """
    invert_epsilon_f = 1 / flow_lower_bound
    return list(
        m.addConstrs(
            (
                sub_arc_vars.t(a) - invert_epsilon_f * flow_vars.t(a) <= 0
                for a in network.sink_arcs()
            ),
            name="active_sink_arc_has_strict_positive_flow",
        ).values(),
    )


# ------------------------------------------------------------------------------------ #
def fragment_coverages_limit_cumulative_flows(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
) -> list[gp.Constr]:
    """Fragment coverages limit cumulative flows."""
    constraints: list[gp.Constr] = []
    constraints.extend(
        m.addConstrs(
            (
                flow_vars.incoming_forward_reverse(network, frag_id)
                <= network.coverage(frag_id)
                for frag_id in network.gfa_graph().segment_names
            ),
            name="fragment_coverages_limit_cumulative_flows",
        ).values(),
    )
    return constraints


def flow_conservation(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
) -> list[gp.Constr]:
    """Flow conservation."""
    constraints: list[gp.Constr] = []
    constraints.extend(
        m.addConstrs(
            (
                flow_vars.incoming(network, fragment)
                == flow_vars.outgoing(network, fragment)
                for fragment in network.oriented_fragments()
            ),
            name="flow_conservation",
        ).values(),
    )
    return constraints


def total_flow_value(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
) -> list[gp.Constr]:
    """Set the total flow equals the source outgoing flow and the sink incoming flow."""
    constraints: list[gp.Constr] = []
    constraints.append(
        m.addConstr(
            flow_vars.total()
            == gp.quicksum(flow_vars.s(s_w) for s_w in network.source_arcs()),
            name="total_flow_equals_source_outgoing_flow",
        ),
    )
    constraints.append(
        m.addConstr(
            flow_vars.total()
            == gp.quicksum(flow_vars.t(u_t) for u_t in network.sink_arcs()),
            name="total_flow_equals_sink_incoming_flow",
        ),
    )
    return constraints


def active_arcs_implies_active_fragments(
    m: gp.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Active arcs implies active fragments."""
    constraints: list[gp.Constr] = []
    constraints.extend(
        m.addConstrs(
            (
                sub_arc_vars.s(source_arc) <= sub_v_vars.x(source_arc[1])
                for source_arc in network.source_arcs()
            ),
            name="active_source_arc_implies_active_fragment",
        ).values(),
    )
    constraints.extend(
        m.addConstrs(
            (
                sub_arc_vars.l(arc) <= sub_v_vars.x(arc.predecessor())
                for arc in network.link_arcs()
            ),
            name="active_link_arc_implies_active_predecessor_fragment",
        ).values(),
    )
    constraints.extend(
        m.addConstrs(
            (
                sub_arc_vars.l(arc) <= sub_v_vars.x(arc.successor())
                for arc in network.link_arcs()
            ),
            name="active_link_arc_implies_active_successor_fragment",
        ).values(),
    )
    constraints.extend(
        m.addConstrs(
            (
                sub_arc_vars.t(sink_arc) <= sub_v_vars.x(sink_arc[0])
                for sink_arc in network.sink_arcs()
            ),
            name="active_sink_arc_implies_active_fragment",
        ).values(),
    )
    return constraints


# REFACTOR Potentially useless constraint
def active_fragments_imply_at_least_one_active_arc(
    m: gp.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Active fragments imply at least one active arc."""
    constraints: list[gp.Constr] = []
    constraints.extend(
        m.addConstrs(
            (
                sub_v_vars.x(frag) <= sub_arc_vars.incoming(network, frag)
                for frag in network.oriented_fragments()
            ),
            name="active_fragments_imply_at_least_one_active_arc",
        ).values(),
    )
    return constraints


def active_arcs_imply_vertices_in_same_component(
    m: gp.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    sub_arc_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Active arcs imply vertices in the same connected component."""
    constraints: list[gp.Constr] = []
    # REFACTOR potentially redondant with active_fragments_imply_at_least_one_active_arc
    constraints.extend(
        m.addConstrs(
            (
                1 - sub_arc_vars.s(source_arc) >= 1 - sub_v_vars.x(source_arc[1])
                for source_arc in network.source_arcs()
            ),
            name="active_source_arcs_imply_vertices_in_same_component",
        ).values(),
    )
    constraints.extend(
        m.addConstrs(
            (
                1 - sub_arc_vars.l(arc)
                >= sub_v_vars.x(arc.predecessor()) - sub_v_vars.x(arc.successor())
                for arc in network.link_arcs()
            ),
            name="active_link_arcs_imply_vertices_in_same_component",
        ).values(),
    )
    constraints.extend(
        m.addConstrs(
            (
                1 - sub_arc_vars.t(sink_arc) >= sub_v_vars.x(sink_arc[0]) - 1
                for sink_arc in network.sink_arcs()
            ),
            name="active_sink_arcs_imply_vertices_in_same_component",
        ).values(),
    )
    return constraints


def subtree_depth_min_distance(
    m: gp.Model,
    ccomp_vars: pb_lp_var.ConnectedComponent,
    network: net.Network,
) -> list[gp.Constr]:
    """Subtree depth min distance."""
    constraints: list[gp.Constr] = []
    #
    # For the source
    #
    constraints.append(
        m.addConstr(
            ccomp_vars.alpha()
            + gp.quicksum(
                ccomp_vars.beta_s(source_arcs) for source_arcs in network.source_arcs()
            )
            <= 1,
            name="source_subtree_depth_min_distance",
        ),
    )
    #
    # For the fragments
    #
    constraints.extend(
        m.addConstrs(
            (
                ccomp_vars.outgoing_beta(network, fragment)
                - ccomp_vars.incoming_beta(network, fragment)
                <= 1
                for fragment in network.oriented_fragments()
            ),
            name="subtree_depth_min_distance",
        ).values(),
    )
    #
    # For the sink
    #
    constraints.append(
        m.addConstr(
            -gp.quicksum(
                ccomp_vars.beta_t(sink_arcs) for sink_arcs in network.sink_arcs()
            )
            <= 1,
            name="sink_subtree_depth_min_distance",
        ),
    )
    return constraints


def arcs_in_tree_are_active(
    m: gp.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    ccomp_vars: pb_lp_var.ConnectedComponent,
    network: net.Network,
) -> list[gp.Constr]:
    """Arcs in tree are active."""
    constraints: list[gp.Constr] = []
    #
    # Source-arcs
    #
    constraints.extend(
        m.addConstrs(
            (
                ccomp_vars.beta_s(source_arc)
                >= -sub_arc_vars.s(source_arc) * network.number_of_vertices()
                for source_arc in network.source_arcs()
            ),
            name="source_arcs_in_tree_are_active",
        ).values(),
    )
    #
    # Link-arcs
    #
    constraints.extend(
        m.addConstrs(
            (
                ccomp_vars.beta(arc)
                >= -sub_arc_vars.l(arc) * network.number_of_vertices()
                for arc in network.link_arcs()
            ),
            name="Link_arcs_in_tree_are_active",
        ).values(),
    )
    #
    # Sink-arcs
    #
    constraints.extend(
        m.addConstrs(
            (
                ccomp_vars.beta_t(sink_arc)
                >= -sub_arc_vars.t(sink_arc) * network.number_of_vertices()
                for sink_arc in network.sink_arcs()
            ),
            name="sink_arcs_in_tree_are_active",
        ).values(),
    )
    return constraints


def alpha_is_the_number_of_vertices_in_the_solution(
    m: gp.Model,
    sub_v_vars: pb_lp_var.SubVertices,
    ccomp_vars: pb_lp_var.ConnectedComponent,
    network: net.Network,
) -> list[gp.Constr]:
    """Alpha is the number of vertices in the solution connected subgraph."""
    constraints: list[gp.Constr] = []
    constraints.append(
        m.addConstr(
            ccomp_vars.alpha()
            == 2
            + gp.quicksum(
                sub_v_vars.x(fragment) for fragment in network.oriented_fragments()
            ),
            name="alpha_is_the_number_of_vertices_in_the_solution",
        ),
    )
    return constraints


def active_arcs_have_flow_at_least_total_flow(
    m: gp.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    flow_vars: pb_lp_var.Flow,
    pos_flow_vars: pb_lp_var.PositiveFlow,
    network: net.Network,
) -> list[gp.Constr]:
    """Active links have flow at least total flow."""
    constraints: list[gp.Constr] = []
    max_seed_cov = max(
        network.coverage(frag_id) for frag_id in network.source_connected_fragment_ids()
    )
    constraints.extend(
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
                        >= flow_vars.total()
                        - (1 - sub_arc_vars.t(sink_arc)) * max_seed_cov
                        for sink_arc in network.sink_arcs()
                    ),
                )
            ),
            name="active_links_have_flow_at_least_total_flow_1",
        ).values(),
    )
    constraints.extend(
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
        ).values(),
    )
    constraints.extend(
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
        ).values(),
    )
    return constraints


def total_flow_is_strictly_positive(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    epsilon_total_flow: float,
) -> list[gp.Constr]:
    """Total flow is strictly positive."""
    # DOCU total flow can be the coverage of plasmid repeats
    # FIXME now if want a unique region of a plasmid to have a lb flow
    # And the total flow does not represent unique region coverages
    constraints: list[gp.Constr] = []
    constraints.append(
        m.addConstr(
            flow_vars.total() >= epsilon_total_flow,
            name="total_flow_is_strictly_positive",
        ),
    )
    return constraints


def minimum_cumulative_length(
    m: gp.Model,
    frag_vars: pb_lp_var.SubFragments,
    network: net.Network,
    minimum_cumulative_length: int,
) -> gp.Constr:
    """Minimum cumulative length."""
    # DOCU note that the repetition is not considered
    return m.addConstr(
        gp.quicksum(
            frag_vars.frag(frag_id)
            * gfa_segment.length(network.gfa_graph().segment(frag_id))
            for frag_id in network.fragment_ids()
        )
        >= minimum_cumulative_length,
        name="minimum_cumulative_length",
    )


# ------------------------------------------------------------------------------------ #
#                                          GC                                          #
# ------------------------------------------------------------------------------------ #
def exactly_one_interval_is_active(
    m: gp.Model,
    gc_vars: pb_lp_var.GCIntervals,
    intervals: gc_items.Intervals,
) -> gp.Constr:
    """Exactly one interval is active."""
    return m.addConstr(
        gp.quicksum(gc_vars.x(interval) for interval in intervals) == 1,
        name="exactly_one_interval_is_active",
    )


def define_frag_gc(  # noqa: PLR0913
    m: gp.Model,
    frag_vars: pb_lp_var.SubFragments,
    gc_vars: pb_lp_var.GCIntervals,
    frag_gc_vars: pb_lp_var.FragmentGC,
    network: net.Network,
    intervals: gc_items.Intervals,
) -> list[gp.Constr]:
    """Define frag gc variables."""
    return list(
        chain(
            m.addConstrs(
                (
                    frag_gc_vars.x(frag_id, interval) <= frag_vars.frag(frag_id)
                    for (frag_id, interval) in product(
                        network.fragment_ids(),
                        intervals,
                    )
                ),
                name="define_frag_gc_1",
            ).values(),
            m.addConstrs(
                (
                    frag_gc_vars.x(frag_id, interval) <= gc_vars.x(interval)
                    for (frag_id, interval) in product(
                        network.fragment_ids(),
                        intervals,
                    )
                ),
                name="define_frag_gc_2",
            ).values(),
            m.addConstrs(
                (
                    frag_gc_vars.x(frag_id, interval)
                    >= frag_vars.frag(frag_id) + gc_vars.x(interval) - 1
                    for (frag_id, interval) in product(
                        network.fragment_ids(),
                        intervals,
                    )
                ),
                name="define_frag_gc_3",
            ).values(),
        ),
    )


def define_inflow_gc(  # noqa: PLR0913
    m: gp.Model,
    frag_vars: pb_lp_var.SubFragments,
    gc_vars: pb_lp_var.GCIntervals,
    flow_vars: pb_lp_var.Flow,
    inflow_gc_vars: pb_lp_var.InflowGC,
    network: net.Network,
    intervals: gc_items.Intervals,
) -> list[gp.Constr]:
    """Define inflow_gc variables."""
    return list(
        chain(
            m.addConstrs(
                (
                    inflow_gc_vars.x(frag_id, interval)
                    <= frag_vars.frag(frag_id) * network.coverage(frag_id)
                    for (frag_id, interval) in product(
                        network.fragment_ids(),
                        intervals,
                    )
                ),
                name="define_inflow_gc_1",
            ).values(),
            m.addConstrs(
                (
                    inflow_gc_vars.x(frag_id, interval)
                    <= gc_vars.x(interval) * network.coverage(frag_id)
                    for (frag_id, interval) in product(
                        network.fragment_ids(),
                        intervals,
                    )
                ),
                name="define_inflow_gc_2",
            ).values(),
            m.addConstrs(
                (
                    inflow_gc_vars.x(frag_id, interval)
                    >= flow_vars.incoming_forward_reverse(network, frag_id)
                    - (2 - frag_vars.frag(frag_id) - gc_vars.x(interval))
                    * network.coverage(frag_id)
                    for (frag_id, interval) in product(
                        network.fragment_ids(),
                        intervals,
                    )
                ),
                name="define_inflow_gc_3",
            ).values(),
        ),
    )


# ------------------------------------------------------------------------------------ #
#                                 Seed Tree Constraints                                #
# ------------------------------------------------------------------------------------ #
def arc_in_seed_tree_is_active(
    m: gp.Model,
    seed_tree_vars: pb_lp_var.SeedTreeArcs,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set arc in seed tree is active."""
    number_of_seeds = len(network.seeds())
    return list(
        chain(
            m.addConstrs(
                (
                    seed_tree_vars.zeta_s(source_arc)
                    <= sub_arcs_vars.s(source_arc) * number_of_seeds
                    for source_arc in network.source_arcs()
                ),
                name="source_link_in_seed_tree_is_active",
            ).values(),
            m.addConstrs(
                (
                    seed_tree_vars.zeta(link) <= sub_arcs_vars.l(link) * number_of_seeds
                    for link in network.link_arcs()
                ),
                name="link_in_seed_tree_is_active",
            ).values(),
            m.addConstrs(
                (
                    seed_tree_vars.zeta_t(sink_link)
                    <= sub_arcs_vars.t(sink_link) * number_of_seeds
                    for sink_link in network.sink_arcs()
                ),
                name="sink_link_in_seed_tree_is_active",
            ).values(),
        ),
    )


def size_of_seed_tree_under_the_source(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    seed_tree_vars: pb_lp_var.SeedTreeArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set size of seed tree under the source."""
    return [
        m.addConstr(
            gp.quicksum(
                sub_vertices_vars.x(oriented_seed)
                for seed in network.seeds()
                for oriented_seed in network.to_oriented(seed)
            )
            - gp.quicksum(
                seed_tree_vars.zeta_s(source_arc)
                for source_arc in network.source_arcs()
            )
            == 0,
        ),
    ]


def size_of_seed_subtree(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    seed_tree_vars: pb_lp_var.SeedTreeArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set size of seed subtree."""
    return list(
        chain(
            m.addConstrs(
                seed_tree_vars.incoming_zeta(network, fragment)
                - seed_tree_vars.outgoing_zeta(network, fragment)
                == (
                    sub_vertices_vars.x(fragment)
                    if fragment.identifier() in network.seeds()
                    else 0
                )
                for fragment in network.oriented_fragments()
            ).values(),
            m.addConstrs(
                gp.quicksum(
                    seed_tree_vars.zeta_t(sink_arc) for sink_arc in network.sink_arcs()
                )
                == 0
                for _ in (network.SINK_VERTEX,)
            ).values(),
        ),
    )


def active_source_arc_implies_subgraph_contains_seed(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    seed_tree_vars: pb_lp_var.SeedTreeArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set active source arc implies subgraph contains seeds."""
    return list(
        m.addConstrs(
            sub_arcs_vars.s(source_arc) <= seed_tree_vars.zeta_s(source_arc)
            for source_arc in network.source_arcs()
        ).values(),
    )


# ------------------------------------------------------------------------------------ #
#                                  Network Properties                                  #
# ------------------------------------------------------------------------------------ #
def no_more_starter_than_active_seeds(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    sub_frag_vars: pb_lp_var.SubFragments,
    network: net.Network,
) -> gp.Constr:
    """Set no more starter than seeds."""
    return m.addConstr(
        gp.quicksum(sub_arcs_vars.s(a) for a in network.source_arcs())
        <= gp.quicksum(sub_frag_vars.frag(seed_id) for seed_id in network.seeds()),
        name="No_more_starter_than_active_seeds",
    )


def active_seeds_lower_bound(
    m: gp.Model,
    sub_frag_vars: pb_lp_var.SubFragments,
    network: net.Network,
    min_number_of_seeds: int,
) -> gp.Constr:
    """Set at least one seed in ccomp."""
    return m.addConstr(
        gp.quicksum(sub_frag_vars.frag(seed_id) for seed_id in network.seeds())
        >= min_number_of_seeds,
        name="active_seeds_lb",
    )


def source_arc_flow_upper_bound(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
    max_source_arc_flow: float,
) -> list[gp.Constr]:
    """Set source arc flow upper bound."""
    return list(
        m.addConstrs(
            (
                flow_vars.s(source_arc) <= max_source_arc_flow
                for source_arc in network.source_arcs()
            ),
            name="source_arc_flow_ub",
        ).values(),
    )


def sink_arc_flow_upper_bound(
    m: gp.Model,
    flow_vars: pb_lp_var.Flow,
    network: net.Network,
    max_sink_arc_flow: float,
) -> list[gp.Constr]:
    """Set sink arc flow upper bound."""
    return list(
        m.addConstrs(
            (
                flow_vars.t(sink_arc) <= max_sink_arc_flow
                for sink_arc in network.sink_arcs()
            ),
            name="sink_arc_flow_ub",
        ).values(),
    )


def same_oriented_fragment_connects_s_and_t_eq(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set same oriented fragment connects s and t equality."""
    return list(
        m.addConstrs(
            (
                sub_arcs_vars.s((network.SOURCE_VERTEX, or_frag))
                - sub_arcs_vars.t((or_frag, network.SINK_VERTEX))
                == 0
                for frag_id in network.source_connected_fragment_ids()
                for or_frag in network.to_oriented(frag_id)
            ),
            name="same_oriented_fragment_connects_s_and_t_eq",
        ).values(),
    )


def same_seed_connects_s_and_t_lb(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set same seed connects s and t lower bound.

    Notes
    -----
    (Default) activate the constraint: set RHS to `0`.
    Deactivate the constraint: set RHS to `-1`.
    """
    return list(
        m.addConstrs(
            (
                sub_arcs_vars.s((network.SOURCE_VERTEX, orient_seed))
                - sub_arcs_vars.t((orient_seed, network.SINK_VERTEX))
                >= 0
                for seed_id in network.seeds()
                for orient_seed in network.to_oriented(seed_id)
            ),
            name="same_seed_connects_s_and_t_lb",
        ).values(),
    )


def same_seed_connects_s_and_t_ub(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set same seed connects s and t upper bound.

    Notes
    -----
    (Default) activate the constraint: set RHS to `0`.
    Deactivate the constraint: set RHS to `1`.
    """
    return list(
        m.addConstrs(
            (
                sub_arcs_vars.s((network.SOURCE_VERTEX, orient_seed))
                - sub_arcs_vars.t((orient_seed, network.SINK_VERTEX))
                <= 0
                for seed_id in network.seeds()
                for orient_seed in network.to_oriented(seed_id)
            ),
            name="same_seed_connects_s_and_t_ub",
        ).values(),
    )


def s_connected_orfrag_incoming_arcs_ub(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> tuple[list[gp.Constr], list[int]]:
    """Set s connected orfrag incoming arcs upper bound.

    Returns
    -------
    list[gp.Constr]
        List of constraints.
    list[int]
        List of number of predecessors for the oriented fragments
        connected to the source. Internal data structure where the order
        of the items is the same as the constraints.

    Notes
    -----
    (Default) authorize only the source arc or all but source arc:
    set RHS to `1 + nb_preds[k]`.
    No constraint: set RHS to `2 * nb_preds[k]`.
    """
    nb_preds = [
        network.number_predecessors(orient_frag)
        for frag_id in network.source_connected_fragment_ids()
        for orient_frag in network.to_oriented(frag_id)
    ]
    return list(
        m.addConstrs(
            (
                sub_arcs_vars.incoming(network, orient_frag)
                + sub_arcs_vars.s((network.SOURCE_VERTEX, orient_frag)) * nb_preds[k]
                <= 1 + nb_preds[k]
                # DOCU modif to put cst on one side
                # sum_(y_uv) - (1 - y_sv) A-v <= 0
                for k, orient_frag in enumerate(
                    orient_frag
                    for frag_id in network.source_connected_fragment_ids()
                    for orient_frag in network.to_oriented(frag_id)
                )
            ),
            name="s_connected_orfrag_incoming_arcs_ub",
        ).values(),
    ), nb_preds


def t_connected_orfrag_outgoing_arcs_ub(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> tuple[list[gp.Constr], list[int]]:
    """Set t connected orfrag outgoing arcs upper bound.

    Returns
    -------
    list[gp.Constr]
        List of constraints.
    list[int]
        List of number of successors for the oriented fragments
        connected to the sink. Internal data structure where the order
        of the items is the same as the constraints.

    Notes
    -----
    (Default) authorize only the sink arc or all but sink arc:
    set RHS to `1 + nb_succs[k]`.
    No constraint: set RHS to `2 * nb_succs[k]`.
    """
    nb_succs = [
        network.number_successors(orient_frag)
        for frag_id in network.sink_connected_fragment_ids()
        for orient_frag in network.to_oriented(frag_id)
    ]
    return list(
        m.addConstrs(
            (
                sub_arcs_vars.outgoing(network, orient_frag)
                + sub_arcs_vars.t((orient_frag, network.SINK_VERTEX)) * nb_succs[k]
                <= 1 + nb_succs[k]
                for k, orient_frag in enumerate(
                    orient_frag
                    for frag_id in network.sink_connected_fragment_ids()
                    for orient_frag in network.to_oriented(frag_id)
                )
            ),
            name="t_connected_orfrag_outgoing_arcs_ub",
        ).values(),
    ), nb_succs


def cycle_before_out(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set cycle before out to the target.

    Notes
    -----
    (Default) activate the constraint: set RHS to `0`.
    Deactivate the constraint: set RHS to `1`.
    """
    return list(
        m.addConstrs(
            (
                sub_vertices_vars.x(or_frag)
                - gp.quicksum(
                    sub_arcs_vars.l(out_link_arc)
                    for out_link_arc in network.out_link_arcs(or_frag)
                )
                <= 0
                for frag_id in network.fragment_ids()
                for or_frag in network.to_oriented(frag_id)
            ),
            name="cycle_before_out",
        ).values(),
    )


def cycle_before_in(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set cycle before in from the source.

    Notes
    -----
    (Default) activate the constraint: set RHS to `0`.
    Deactivate the constraint: set RHS to `1`.
    """
    return list(
        m.addConstrs(
            (
                sub_vertices_vars.x(or_frag)
                - gp.quicksum(
                    sub_arcs_vars.l(in_link_arc)
                    for in_link_arc in network.in_link_arcs(or_frag)
                )
                <= 0
                for frag_id in network.fragment_ids()
                for or_frag in network.to_oriented(frag_id)
            ),
            name="cycle_before_in",
        ).values(),
    )


def source_arcs_upper_bound(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
    max_number_of_source_arcs: int,
) -> gp.Constr:
    """Set number of source arcs upper bound."""
    return m.addConstr(
        gp.quicksum(sub_arcs_vars.s(a) for a in network.source_arcs())
        <= max_number_of_source_arcs,
        name="source_arcs_ub",
    )


def sink_arcs_upper_bound(
    m: gp.Model,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
    max_number_of_sink_arcs: int,
) -> gp.Constr:
    """Set number of sink arcs upper bound."""
    return m.addConstr(
        gp.quicksum(sub_arcs_vars.t(a) for a in network.sink_arcs())
        <= max_number_of_sink_arcs,
        name="sink_arcs_ub",
    )
