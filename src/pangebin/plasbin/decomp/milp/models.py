"""PangeBin-flow hierarchical decomposition MILP models."""

from enum import StrEnum

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.decomp.milp.constraints as lp_cst
import pangebin.plasbin.decomp.milp.objectives as lp_obj
import pangebin.plasbin.decomp.milp.variables as lp_var
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net


class Names(StrEnum):
    """PangeBin-flow MILP model names."""

    MCF = "MCF"
    MGC = "MGC"
    MPS = "MPS"
    MRCF = "MRCF"


def mcf(
    network: net.Network,
    min_flow: float,
    min_cumulative_len: int,
    circular: bool,  # noqa: FBT001
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.MaxCovFlow]:
    """Create MCF model."""
    mcf_model = gp.Model("Maximum Coverage likelihood Flow")
    mcf_vars = lp_var.init_mcf(network, mcf_model)
    lp_obj.set_mcf_objective(mcf_model, mcf_vars, network, obj_fun_domain)
    lp_cst.set_mcf_constraints(
        mcf_model,
        mcf_vars,
        network,
        min_flow,
        min_cumulative_len,
        circular,
    )
    return mcf_model, mcf_vars


def mgc_from_mcf(  # noqa: PLR0913
    mcf_model: gp.Model,
    mcf_vars: lp_var.MaxCovFlow,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    network: net.Network,
    intervals: gc_items.Intervals,
    gamma_mcf: float,
) -> tuple[gp.Model, lp_var.MaxGC]:
    """Create MGC model from MCF model."""
    previous_coverage_score = mcf_model.ObjVal
    mgc_model = mcf_model
    mgc_model.ModelName = "Maximum GC Score"
    mgc_vars = lp_var.mgc_from_mcf(mgc_model, mcf_vars, network, intervals)
    lp_obj.set_mgc_objective(mgc_model, mgc_vars, network, intervals, obj_fun_domain)
    lp_cst.add_mgc_constraints(
        mgc_model,
        mgc_vars,
        obj_fun_domain,
        network,
        intervals,
        gamma_mcf,
        previous_coverage_score,
    )
    return mgc_model, mgc_vars


def mps_from_mgc(  # noqa: PLR0913
    mgc_model: gp.Model,
    mgc_vars: lp_var.MaxGC,
    network: net.Network,
    intervals: gc_items.Intervals,
    gamma_mgc: float,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.MaxPlasmidScore]:
    """Create MPS model from MGC model."""
    previous_gc_score = mgc_model.ObjVal
    mps_model = mgc_model
    mps_model.ModelName = "Maximum Plasmidness Score"
    mps_vars = lp_var.mps_from_mgc(mgc_vars)
    lp_obj.set_mps_objective(mps_model, mps_vars, network, obj_fun_domain)
    lp_cst.add_mps_constraints(
        mps_model,
        mps_vars,
        network,
        intervals,
        gamma_mgc,
        previous_gc_score,
        obj_fun_domain,
    )
    return mps_model, mps_vars


def mrcf_from_mps(
    mps_model: gp.Model,
    mps_var: lp_var.MaxPlasmidScore,
    network: net.Network,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> tuple[gp.Model, lp_var.MaxRefCovFlow]:
    """Create MRCF model from MPS model."""
    # DOCU MPS' renamed to MRCF
    previous_plasmidness_score = mps_model.ObjVal
    mrcf_model = mps_model
    mrcf_model.ModelName = "Maximum Refined Coverage likelihood Flow"
    mrcf_vars = lp_var.mrcf_from_mps(mps_var)
    lp_obj.set_mrcf_objective(mrcf_model, mrcf_vars, obj_fun_domain, network)
    lp_cst.add_mrcf_constraints(
        mrcf_model,
        mrcf_vars,
        network,
        previous_plasmidness_score,  # BUG numerical value approx when .GetValue
        obj_fun_domain,
    )
    return mrcf_model, mrcf_vars
