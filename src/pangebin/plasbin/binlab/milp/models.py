"""PangeBin-flow binlab MILP models."""

# DOCU binlab models

from enum import StrEnum

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.binlab.milp.constraints as lp_cst
import pangebin.plasbin.binlab.milp.objectives as lp_obj
import pangebin.plasbin.binlab.milp.variables as lp_var
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net


class Names(StrEnum):
    """PangeBin-flow binlab MILP model names."""

    MBS = "MBS"
    MLS = "MLS"
    MRBS = "MRBS"


def mbs(
    network: net.Network,
    min_flow: float,
    min_cumulative_len: int,
    circular: bool,  # noqa: FBT001
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.MaxBinScore]:
    """Create MBS model."""
    m = gp.Model("Maximum Binning Score")
    var = lp_var.init_mbs(network, m)
    lp_obj.set_mbs_objective(m, var, network, obj_fun_domain)
    lp_cst.set_mbs_constraints(
        m,
        var,
        network,
        min_flow,
        min_cumulative_len,
        circular,
    )
    return m, var


def mls_from_mbs(  # noqa: PLR0913
    mbs_model: gp.Model,
    mbs_var: lp_var.MaxBinScore,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    network: net.Network,
    intervals: gc_items.Intervals,
    gamma_mbs: float,
) -> tuple[gp.Model, lp_var.MaxLabScore]:
    """Create MLS model from MBS model."""
    previous_binning_score = mbs_model.ObjVal
    mls_model = mbs_model
    mls_model.ModelName = "Maximum Labelling Score"
    mls_vars = lp_var.mls_from_mbs(mls_model, mbs_var, network, intervals)
    lp_obj.set_mls_objective(mls_model, mls_vars, network, intervals, obj_fun_domain)
    lp_cst.add_mls_constraints(
        mls_model,
        mls_vars,
        obj_fun_domain,
        network,
        intervals,
        gamma_mbs,
        previous_binning_score,
    )
    return mls_model, mls_vars


def mrbs_from_mls(
    mls_model: gp.Model,
    mls_var: lp_var.MaxLabScore,
    network: net.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.MaxRefBinScore]:
    """Create MRBS model from MLS model."""
    previous_labelling_score = mls_model.ObjVal
    mrbs_model = mls_model
    mrbs_model.ModelName = "Maximum Refined Binnings Score"
    mrbs_vars = lp_var.mrbs_from_mls(mls_var)
    lp_obj.set_mrbs_objective(mrbs_model, mrbs_vars, obj_fun_domain, network)
    lp_cst.add_mrbs_constraints(
        mrbs_model,
        mrbs_vars,
        network,
        intervals,
        previous_labelling_score,
        obj_fun_domain,
    )
    return mrbs_model, mrbs_vars
