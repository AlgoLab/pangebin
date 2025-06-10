"""MILP model."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net

from . import constraints as mfb_cst
from . import objectives as mfb_obj
from . import variables as mfb_var


def new(
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    min_flow: float,
    plasmidness_coefficient: float,
    min_cumulative_len: int,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[
    gp.Model,
    list[mfb_var.BinVariables],
    pb_lp_var.SubFragments,
    mfb_cst.MultiFlowStateConstraints,
]:
    """Create Classify model."""
    m = gp.Model("Maximum Classification Flow")
    bins_vars, flow_union_frag_vars = mfb_var.init_binning(network, m)
    mfb_obj.set_objective(
        m,
        bins_vars,
        flow_union_frag_vars,
        network,
        gc_intervals,
        obj_fun_domain,
    )
    mfb_state_csts = mfb_cst.set_constraints(
        m,
        bins_vars,
        flow_union_frag_vars,
        network,
        min_flow,
        plasmidness_coefficient,
        obj_fun_domain,
        min_cumulative_len,
    )
    return m, bins_vars, flow_union_frag_vars, mfb_state_csts
