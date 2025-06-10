"""PangeBin-flow binlab MILP objectives."""

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.binlab.milp.variables as lp_vars
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net


def binning_score(
    network: net.Network,
    flow_vars: pb_lp_var.Flow,
    frag_vars: pb_lp_var.SubFragments,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for binning score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = net.max_frag_length(network, frag_set_fn)
    return gurobipy.quicksum(
        pb_lp_obj.zeta_i(network, frag_id, max_frag_length)
        * (
            flow_vars.incoming_forward_reverse(network, frag_id)
            - (
                network.coverage(frag_id) * frag_vars.frag(frag_id)
                - flow_vars.incoming_forward_reverse(network, frag_id)
            )
            + network.plasmidness(frag_id)
            * flow_vars.incoming_forward_reverse(network, frag_id)
        )
        for frag_id in frag_set_fn(network)
    )


def set_mbs_objective(
    m: gurobipy.Model,
    var: lp_vars.MaxBinScore,
    network: net.Network,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MBS objective."""
    m.setObjective(
        binning_score(network, var.flow(), var.frag(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def gc_score(
    network: net.Network,
    intervals: gc_items.Intervals,
    frag_gc_vars: pb_lp_var.FragmentGC,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for GC score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gurobipy.quicksum(
        network.gc_score(frag_id)[b] * frag_gc_vars.x(frag_id, interval)
        for b, interval in enumerate(intervals)
        for frag_id in frag_set_fn(network)
    )


def set_mls_objective(
    m: gurobipy.Model,
    var: lp_vars.MaxLabScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MLS objective."""
    m.setObjective(
        gc_score(network, intervals, var.frag_gc(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def set_mrbs_objective(
    m: gurobipy.Model,
    var: lp_vars.MaxRefBinScore,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    network: net.Network,
) -> None:
    """Set MRBS objective."""
    m.setObjective(
        binning_score(network, var.flow(), var.frag(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )
