"""PangeBin-flow classbin MILP models."""

# DOCU classbin models

from enum import StrEnum

import gurobipy as gp

import pangebin.plasbin.classbin.milp.constraints as lp_cst
import pangebin.plasbin.classbin.milp.objectives as lp_obj
import pangebin.plasbin.classbin.milp.variables as lp_var
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net


class Names(StrEnum):
    """PangeBin-flow classbin MILP model names."""

    CLASSIFY = "Classify"
    MFB = "MFB"
    # FIXME change model name


def classify(
    network: net.Network,
    min_flow: float,
    plasmidness_coefficient: float,
    min_cumulative_len: int,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.Classify]:
    """Create Classify model."""
    m = gp.Model("Maximum Classification Flow")
    var = lp_var.init_classify(network, m)
    lp_obj.set_classify_objective(m, var, network, obj_fun_domain)
    lp_cst.set_classify_constraints(
        m,
        var,
        network,
        min_flow,
        plasmidness_coefficient,
        obj_fun_domain,
        min_cumulative_len,
    )
    return m, var
