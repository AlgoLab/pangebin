"""PangeBin-flow once MILP objectives."""

from collections.abc import Callable, Iterable

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net
import pangebin.plasbin.once.milp.variables as lp_vars


def _max_frag_length(
    network: net.Network,
    frag_set_fn: Callable[[net.Network], Iterable[str]],
) -> int:
    return max(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in frag_set_fn(network)
    )


def zeta_i(network: net.Network, frag_id: str, max_frag_length: int) -> float:
    """Get zeta_i coefficient."""
    return gfa_segment.length(network.gfa_graph().segment(frag_id)) / max_frag_length


def coverage_score(
    network: net.Network,
    flow_vars: pb_lp_var.Flow,
    frag_vars: pb_lp_var.SubFragments,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = _max_frag_length(network, frag_set_fn)
    return gp.quicksum(
        zeta_i(network, frag_id, max_frag_length)
        * (
            flow_vars.incoming_forward_reverse(network, frag_id)
            - (
                network.coverage(frag_id) * frag_vars.x(frag_id)
                - flow_vars.incoming_forward_reverse(network, frag_id)
            )
        )
        for frag_id in frag_set_fn(network)
    )


def plasmidness_score(
    network: net.Network,
    flow_vars: pb_lp_var.Flow,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = _max_frag_length(network, frag_set_fn)
    return gp.quicksum(
        zeta_i(network, frag_id, max_frag_length)
        * network.plasmidness(frag_id)
        * flow_vars.incoming_forward_reverse(network, frag_id)
        for frag_id in frag_set_fn(network)
    )


def gc_score(
    network: net.Network,
    intervals: gc_items.Intervals,
    inflow_gc_vars: pb_lp_var.InflowGC,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for GC score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = _max_frag_length(network, frag_set_fn)
    return gp.quicksum(
        zeta_i(network, frag_id, max_frag_length)
        * network.gc_score(frag_id)[b]
        * inflow_gc_vars.x(frag_id, interval)
        for b, interval in enumerate(intervals)
        for frag_id in frag_set_fn(network)
    )


def set_mgclb_objective(
    m: gp.Model,
    var: lp_vars.MaxGCLabelBinScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MGCLB objective."""
    m.setObjective(
        (
            coverage_score(network, var.flow(), var.frag(), obj_fun_domain)
            + plasmidness_score(network, var.flow(), obj_fun_domain)
            + gc_score(network, intervals, var.inflow_gc(), obj_fun_domain)
        ),
        gp.GRB.MAXIMIZE,
    )
