"""PangeBin-flow once MILP models."""

# DOCU once models

from enum import StrEnum

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.network as net
import pangebin.plasbin.once.milp.constraints as lp_cst
import pangebin.plasbin.once.milp.objectives as lp_obj
import pangebin.plasbin.once.milp.variables as lp_var


class Names(StrEnum):
    """PangeBin-flow once MILP model names."""

    MGCLB = "MGCLB"


def mgclb(  # noqa: PLR0913
    network: net.Network,
    intervals: gc_items.Intervals,
    min_flow: float,
    min_cumulative_len: int,
    circular: bool,  # noqa: FBT001
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.MaxGCLabelBinScore]:
    """Create MGCLB model."""
    m = gp.Model("Maximum GC Label-Binning Score")
    var = lp_var.init_mgclb(network, intervals, m)
    lp_obj.set_mgclb_objective(m, var, network, intervals, obj_fun_domain)
    lp_cst.set_mgclb_constraints(
        m,
        var,
        network,
        intervals,
        min_flow,
        min_cumulative_len,
        circular,
    )
    return m, var
