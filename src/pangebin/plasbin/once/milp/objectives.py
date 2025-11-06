"""PangeBin-flow once MILP objectives."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net
import pangebin.plasbin.once.milp.variables as lp_vars


def coverage_score(
    network: net.Network,
    flow_vars: cmn_lp_vars.Flow,
    frag_vars: cmn_lp_vars.SubFragments,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        net.length(network, frag_id)
        * (
            flow_vars.incoming_forward_reverse(network, frag_id)
            - (
                network.coverage(frag_id) * frag_vars.frag(frag_id)
                - flow_vars.incoming_forward_reverse(network, frag_id)
            )
        )
        for frag_id in frag_set_fn(network)
    )


def coverage_penalty(
    network: net.Network,
    flow_vars: cmn_lp_vars.Flow,
    frag_vars: cmn_lp_vars.SubFragments,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        net.length(network, frag_id)
        * (
            flow_vars.incoming_forward_reverse(network, frag_id)
            - network.coverage(frag_id) * frag_vars.frag(frag_id)
        )
        for frag_id in frag_set_fn(network)
    )


def plasmidness_score(
    network: net.Network,
    flow_vars: cmn_lp_vars.Flow,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        net.length(network, frag_id)
        * network.plasmidness(frag_id)
        * flow_vars.incoming_forward_reverse(network, frag_id)
        for frag_id in frag_set_fn(network)
    )


def gc_score(
    network: net.Network,
    intervals: gc_items.Intervals,
    inflow_gc_vars: cmn_lp_vars.InflowGC,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for GC score."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        net.length(network, frag_id)
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
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> None:
    """Set MGCLB objective."""
    m.setObjective(
        (
            coverage_penalty(network, var.flows(), var.sub_frags(), obj_fun_domain)
            + plasmidness_score(network, var.flows(), obj_fun_domain)
            + gc_score(network, intervals, var.inflow_gc(), obj_fun_domain)
        ),
        gp.GRB.MAXIMIZE,
    )
