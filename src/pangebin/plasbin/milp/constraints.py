"""PangeBin-flow MILP constraints."""

from __future__ import annotations

from itertools import chain, product

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.network as pb_network


def set_mcf_constraints(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """MCF constraints."""
    # FIXME active arc implies active arc (should not do this normally...)
    _exactly_one_active_source_arc(m, var, network)
    _arc_capacities_limit_arc_flows(m, var, network)
    _fragment_coverages_limit_cumulative_flows(m, var, network)
    _flow_conservation(m, var, network)
    _total_flow_value(m, var, network)
    _active_arcs_implies_active_fragments(m, var, network)
    _active_fragments_imply_at_least_one_active_arc(m, var, network)
    _pos_flow_implies_active_vertices(m, var, network)
    #
    # Connectivity
    #
    _subtree_depth_min_distance(m, var, network)
    _arcs_in_tree_are_active(m, var, network)
    _alpha_is_the_number_of_vertices_in_the_solution(m, var, network)
    #
    # Arc flow lower bound
    #
    _active_arcs_have_flow_at_least_total_flow(m, var, network)
    #
    # Plasmid property
    #
    # REFACTOR magic values should be in config
    _total_flow_is_strictly_positive(m, var, 0.0001)
    _minimum_cumulative_length(m, var, network, 1000)


def _exactly_one_active_source_arc(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Exactly one out source arc constraint."""
    m.addConstr(
        gurobipy.quicksum(var.y_s(a) for a in network.source_arcs()) == 1,
        name="exactly_one_out_source_link",
    )


def _arc_capacities_limit_arc_flows(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Arc capacities limit arc flows."""
    m.addConstrs(
        (
            const
            for const in chain(
                (
                    var.f_s(a) <= network.cap_s(a) * var.y_s(a)
                    for a in network.source_arcs()
                ),
                (var.f(a) <= network.cap(a) * var.y(a) for a in network.link_arcs()),
                (
                    var.f_t(a) <= network.cap_t(a) * var.y_t(a)
                    for a in network.sink_arcs()
                ),
            )
        ),
        name="arc_capacities_limit_arc_flows",
    )


def _fragment_coverages_limit_cumulative_flows(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Fragment coverages limit cumulative flows."""
    m.addConstrs(
        (
            milp_vars.incoming_flow_forward_reverse(frag_id, network, var)
            <= network.coverage(frag_id)
            for frag_id in network.gfa_graph().segment_names
        ),
        name="fragment_coverages_limit_cumulative_flows",
    )


def _flow_conservation(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Flow conservation."""
    m.addConstrs(
        (
            milp_vars.incoming_flow(fragment, network, var)
            == milp_vars.outgoing_flow(fragment, network, var)
            for fragment in network.oriented_fragments()
        ),
        name="flow_conservation",
    )


def _total_flow_value(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Set the total flow equals the source outgoing flow and the sink incoming flow."""
    m.addConstr(
        var.total_flow()
        == gurobipy.quicksum(var.f_s(s_w) for s_w in network.source_arcs()),
        name="total_flow_equals_source_outgoing_flow",
    )
    m.addConstr(
        var.total_flow()
        == gurobipy.quicksum(var.f_t(u_t) for u_t in network.sink_arcs()),
        name="total_flow_equals_sink_incoming_flow",
    )


def _active_arcs_implies_active_fragments(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Active arcs implies active fragments."""
    for source_arc in network.source_arcs():
        m.addConstr(
            var.y_s(source_arc) <= var.x(source_arc[1]),
            name=f"active_arc_{source_arc[0]}_{source_arc[1]}_implies_active_fragment_{source_arc[1]}",
        )
    for arc in network.link_arcs():
        m.addConstr(
            var.y(arc) <= var.x(arc.predecessor()),
            name=f"active_arc_{arc}_implies_active_fragment_{arc.predecessor()}",
        )
        m.addConstr(
            var.y(arc) <= var.x(arc.successor()),
            name=f"active_arc_{arc}_implies_active_fragment_{arc.successor()}",
        )
    for sink_arc in network.sink_arcs():
        m.addConstr(
            var.y_t(sink_arc) <= var.x(sink_arc[0]),
            name=f"active_arc_{sink_arc[0]}_{sink_arc[1]}_implies_active_fragment_{sink_arc[0]}",
        )


# REFACTOR Potentially useless constraint
def _active_fragments_imply_at_least_one_active_arc(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Active fragments imply at least one active arc."""
    m.addConstrs(
        (
            var.x(frag)
            <= milp_vars.incoming_arcs_y(
                frag,
                network,
                var,
            )
            for frag in network.oriented_fragments()
        ),
        name="active_fragments_imply_at_least_one_active_arc",
    )


def _pos_flow_implies_active_vertices(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Positive flow implies the two incident vertices are in the component."""
    for source_arc in network.source_arcs():
        m.addConstr(
            1 - var.y_s(source_arc) >= 1 - var.x(source_arc[1]),
            name=f"pos_flow_implies_active_vertices_{source_arc[0]}_{source_arc[1]}",
        )
    for arc in network.link_arcs():
        m.addConstr(
            1 - var.y(arc) >= var.x(arc.predecessor()) - var.x(arc.successor()),
            name=f"pos_flow_implies_active_vertices_{arc.predecessor()}_{arc.successor()}",
        )
    for sink_arc in network.sink_arcs():
        m.addConstr(
            1 - var.y_t(sink_arc) >= var.x(sink_arc[0]) - 1,
            name=f"pos_flow_implies_active_vertices_{sink_arc[0]}_{sink_arc[1]}",
        )


def _subtree_depth_min_distance(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Subtree depth min distance."""
    #
    # For the source
    #
    m.addConstr(
        var.alpha()
        + gurobipy.quicksum(
            var.beta_s(source_arcs) for source_arcs in network.source_arcs()
        )
        <= 1,
        name="source_subtree_depth_min_distance",
    )
    #
    # For the fragments
    #
    m.addConstrs(
        (
            milp_vars.outgoing_beta(fragment, network, var)
            - milp_vars.incoming_beta(fragment, network, var)
            <= 1
            for fragment in network.oriented_fragments()
        ),
        name="subtree_depth_min_distance",
    )
    #
    # For the sink
    #
    m.addConstr(
        -gurobipy.quicksum(var.beta_t(sink_arcs) for sink_arcs in network.sink_arcs())
        <= 1,
        name="sink_subtree_depth_min_distance",
    )


def _arcs_in_tree_are_active(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Arcs in tree are active."""
    #
    # Source-arcs
    #
    m.addConstrs(
        (
            var.beta_s(source_arc)
            >= -var.y_s(source_arc) * network.number_of_vertices()
            for source_arc in network.source_arcs()
        ),
        name="source_arcs_in_tree_are_active",
    )
    #
    # Link-arcs
    #
    m.addConstrs(
        (
            var.beta(arc) >= -var.y(arc) * network.number_of_vertices()
            for arc in network.link_arcs()
        ),
        name="Link_arcs_in_tree_are_active",
    )
    #
    # Sink-arcs
    #
    m.addConstrs(
        (
            var.beta_t(sink_arc) >= -var.y_t(sink_arc) * network.number_of_vertices()
            for sink_arc in network.sink_arcs()
        ),
        name="sink_arcs_in_tree_are_active",
    )


def _alpha_is_the_number_of_vertices_in_the_solution(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """Alpha is the number of vertices in the solution connected subgraph."""
    m.addConstr(
        var.alpha()
        == 2
        + gurobipy.quicksum(
            var.x(fragment) for fragment in network.oriented_fragments()
        ),
        name="alpha_is_the_number_of_vertices_in_the_solution",
    )


def _active_arcs_have_flow_at_least_total_flow(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
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
                    var.pos_f_s(source_arc)
                    >= var.total_flow() - (1 - var.y_s(source_arc)) * max_seed_cov
                    for source_arc in network.source_arcs()
                ),
                (
                    var.pos_f(arc) >= var.total_flow() - (1 - var.y(arc)) * max_seed_cov
                    for arc in network.link_arcs()
                ),
                (
                    var.pos_f_t(sink_arc)
                    >= var.total_flow() - (1 - var.y_t(sink_arc)) * max_seed_cov
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
                    var.pos_f_s(source_arc) <= var.total_flow()
                    for source_arc in network.source_arcs()
                ),
                (var.pos_f(arc) <= var.total_flow() for arc in network.link_arcs()),
                (
                    var.pos_f_t(sink_arc) <= var.total_flow()
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
                    var.pos_f_s(source_arc) <= var.f_s(source_arc)
                    for source_arc in network.source_arcs()
                ),
                (var.pos_f(arc) <= var.f(arc) for arc in network.link_arcs()),
                (
                    var.pos_f_t(sink_arc) <= var.f_t(sink_arc)
                    for sink_arc in network.sink_arcs()
                ),
            )
        ),
        name="active_links_have_flow_at_least_total_flow_3",
    )


def _total_flow_is_strictly_positive(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    epsilon_total_flow: float,
) -> None:
    """Total flow is strictly positive."""
    # FIXME new constraint
    m.addConstr(
        var.total_flow() >= epsilon_total_flow,
        name="total_flow_is_strictly_positive",
    )


def _minimum_cumulative_length(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
    minimum_cumulative_length: int,
) -> None:
    """Minimum cumulative length."""
    # FIXME new constraints
    for frag_id in network.fragment_ids():
        frag_f = gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.FORWARD)
        frag_r = gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.REVERSE)
        m.addConstr(
            var.x(frag_f) <= var.frag(frag_id),
            name=f"forward_extremity_actives_fragment_{frag_id}",
        )
        m.addConstr(
            var.x(frag_r) <= var.frag(frag_id),
            name=f"reverse_extremity_actives_fragment_{frag_id}",
        )
    m.addConstr(
        minimum_cumulative_length
        <= gurobipy.quicksum(
            var.frag(frag_id) * gfa_segment.length(network.gfa_graph().segment(frag_id))
            for frag_id in network.fragment_ids()
        ),
        name="minimum_cumulative_length",
    )


# ------------------------------------------------------------------------------------ #
#                                          MGC                                         #
# ------------------------------------------------------------------------------------ #
def add_mgc_constraints(  # noqa: PLR0913
    m: gurobipy.Model,
    var: milp_vars.MaxGC,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
    mcf_obj_value: float,
) -> None:
    """Add MGC constraints."""
    _coverage_score_lower_bound(m, var, network, coefficient, mcf_obj_value)
    _exactly_one_interval_is_active(m, var, intervals)
    _define_frag_gc(m, var, network, intervals)


def _coverage_score_lower_bound(
    m: gurobipy.Model,
    var: milp_vars.MaxGC,
    network: pb_network.Network,
    coefficient: float,
    previous_coverage_score: float,
) -> None:
    """Coverage score lower bound."""
    # FIXME be carefull when previous_coverage_score is < 0
    m.addConstr(
        previous_coverage_score - (1 - coefficient) * abs(previous_coverage_score)
        <= milp_vars.coverage_score(network, var.mcf_vars()),
        name="coverage_score_lower_bound",
    )


def _exactly_one_interval_is_active(
    m: gurobipy.Model,
    var: milp_vars.MaxGC,
    intervals: gc_items.Intervals,
) -> None:
    """Exactly one interval is active."""
    m.addConstr(
        gurobipy.quicksum(var.gc(interval) for interval in intervals) == 1,
        name="exactly_one_interval_is_active",
    )


def _define_frag_gc(
    m: gurobipy.Model,
    var: milp_vars.MaxGC,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
) -> None:
    m.addConstrs(
        (
            var.frag_gc(frag_id, interval)
            <= var.mcf_vars().x(
                gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.FORWARD),
            )
            + var.mcf_vars().x(
                gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.REVERSE),
            )
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_1",
    )
    m.addConstrs(
        (
            var.frag_gc(frag_id, interval) <= var.gc(interval)
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_2",
    )
    m.addConstrs(
        (
            var.frag_gc(frag_id, interval)
            >= var.mcf_vars().x(
                gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.FORWARD),
            )
            + var.gc(interval)
            - 1
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_3",
    )
    m.addConstrs(
        (
            var.frag_gc(frag_id, interval)
            >= var.mcf_vars().x(
                gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.REVERSE),
            )
            + var.gc(interval)
            - 1
            for (frag_id, interval) in product(network.fragment_ids(), intervals)
        ),
        name="define_frag_gc_4",
    )


# ------------------------------------------------------------------------------------ #
#                                          MPS                                         #
# ------------------------------------------------------------------------------------ #
def add_mps_constraints(  # noqa: PLR0913
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
    mgc_obj_value: float,
) -> None:
    """Add MPS constraints to MGC model."""
    _gc_probability_score_lower_bound(
        m,
        var,
        network,
        intervals,
        coefficient,
        mgc_obj_value,
    )


def _gc_probability_score_lower_bound(  # noqa: PLR0913
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
    previous_gc_probability_score: float,
) -> None:
    """GC probability score lower bound."""
    # FIXME be carefull when previous_gc_probability_score is < 0
    m.addConstr(
        previous_gc_probability_score
        - (1 - coefficient) * abs(previous_gc_probability_score)
        <= milp_vars.gc_probability_score(network, intervals, var.mgc_vars()),
        name="gc_probability_score_lower_bound",
    )


# ------------------------------------------------------------------------------------ #
#                                         MPS'                                         #
# ------------------------------------------------------------------------------------ #
def add_mps_prime_constraints_to_mps(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    mps_obj_value: float,
) -> None:
    """Add MPS' constraints to MPS model."""
    _fix_plasmidness_score(m, var, network, intervals, mps_obj_value)


def _fix_plasmidness_score(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    mps_obj_value: float,
) -> None:
    m.addConstr(
        milp_vars.plasmidness_score(network, intervals, var) == mps_obj_value,
        name="fix_plasmidness_score",
    )
