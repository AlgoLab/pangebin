"""MILP model."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.milp.variables as cmn_lp_vars
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
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> tuple[
    gp.Model,
    list[mfb_var.BinVariables],
    cmn_lp_vars.SubFragments,
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
