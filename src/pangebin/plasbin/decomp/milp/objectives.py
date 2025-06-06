"""PangeBin-flow MILP objectives."""

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.decomp.milp.variables as lp_var
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net


def coverage_score(
    network: net.Network,
    flow_vars: pb_lp_var.Flow,
    frag_vars: pb_lp_var.SubFragments,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for coverage score."""
    # DOCU MCF: Tests on coverage score
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = pb_lp_obj.max_frag_length(network, frag_set_fn)
    return gurobipy.quicksum(
        pb_lp_obj.zeta_i(network, frag_id, max_frag_length)
        * (
            flow_vars.incoming_forward_reverse(network, frag_id)
            - (
                network.coverage(frag_id) * frag_vars.x(frag_id)
                - flow_vars.incoming_forward_reverse(network, frag_id)
            )
        )
        for frag_id in frag_set_fn(network)
    )


def set_mcf_objective(
    m: gurobipy.Model,
    var: lp_var.MaxCovFlow,
    network: net.Network,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MCF objective."""
    m.setObjective(
        coverage_score(network, var.flow(), var.frag(), obj_fun_domain),
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


def set_mgc_objective(
    m: gurobipy.Model,
    var: lp_var.MaxGC,
    network: net.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MGC objective."""
    m.setObjective(
        gc_score(network, intervals, var.frag_gc(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def plasmidness_score(
    network: net.Network,
    frag_vars: pb_lp_var.SubFragments,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for plasmidness score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = pb_lp_obj.max_frag_length(network, frag_set_fn)
    # DOCU MPS: coeff in obj
    return gurobipy.quicksum(
        pb_lp_obj.zeta_i(network, frag_id, max_frag_length)
        * network.plasmidness(frag_id)
        * frag_vars.x(frag_id)
        for frag_id in frag_set_fn(network)
    )


def set_mps_objective(
    m: gurobipy.Model,
    var: lp_var.MaxPlasmidScore,
    network: net.Network,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MPS objective."""
    m.setObjective(
        plasmidness_score(network, var.frag(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def set_mrcf_objective(
    m: gurobipy.Model,
    var: lp_var.MaxRefCovFlow,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    network: net.Network,
) -> None:
    """Set MRCF objective."""
    # DOCU MRCF: max coverage score
    m.setObjective(
        coverage_score(network, var.flow(), var.frag(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )
