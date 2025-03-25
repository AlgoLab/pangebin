"""PangeBin-flow MILP models."""

from enum import StrEnum

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.constraints as milp_consts
import pangebin.plasbin.milp.objectives as milp_objs
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.network as pb_network


class Names(StrEnum):
    """PangeBin-flow MILP model names."""

    MCF = "MCF"
    MGC = "MGC"
    MPS = "MPS"
    MPS_PRIME = "MPS_PRIME"


def mcf(
    network: pb_network.Network,
) -> tuple[gurobipy.Model, milp_vars.MaxCovFlow]:
    """Create MCF model."""
    m = gurobipy.Model("Maximum Coverage likelihood Flow")
    var = milp_vars.MaxCovFlow(m, network)
    milp_objs.set_mcf_objective(m, var, network)
    milp_consts.set_mcf_constraints(m, var, network)
    return m, var


def mgc_from_mcf(
    mcf_model: gurobipy.Model,
    mcf_var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
) -> tuple[gurobipy.Model, milp_vars.MaxGC]:
    """Create MGC model from MCF model."""
    mcf_model.ModelName = "Maximum GC Probability Score"
    var = milp_vars.MaxGC(mcf_model, mcf_var, network, intervals)
    milp_objs.set_mgc_objective(mcf_model, var, network, intervals)
    milp_consts.add_mgc_constraints(
        mcf_model,
        var,
        network,
        intervals,
        coefficient,
        mcf_model.ObjVal,
    )
    return mcf_model, var


def mps_from_mgc(
    mgc_model: gurobipy.Model,
    mgc_var: milp_vars.MaxGC,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    coefficient: float,
) -> tuple[gurobipy.Model, milp_vars.MaxPlasmidScore]:
    """Create MPS model from MGC model."""
    mgc_model.ModelName = "Maximum Plasmidness Score"
    var = milp_vars.MaxPlasmidScore(mgc_var)
    milp_objs.set_mps_objective(mgc_model, var, network, intervals)
    milp_consts.add_mps_constraints(
        mgc_model,
        var,
        network,
        intervals,
        coefficient,
        mgc_model.ObjVal,
    )
    return mgc_model, var


def mps_prime_from_mps(
    mps_model: gurobipy.Model,
    mps_var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
) -> tuple[gurobipy.Model, milp_vars.MaxPlasmidScore]:
    """Create MPS' model from MPS model."""
    mps_model.ModelName = "Maximum Plasmidness Score'"
    milp_objs.set_mps_prime_objective(mps_model, mps_var, network, intervals)
    milp_consts.add_mps_prime_constraints_to_mps(
        mps_model,
        mps_var,
        network,
        intervals,
        mps_model.ObjVal,  # BUG numerical value approx when .GetValue
    )
    return mps_model, mps_var


def multiobjective(
    network: pb_network.Network,
    intervals: gc_items.Intervals,
) -> tuple[gurobipy.Model, milp_vars.MaxPlasmidScore]:
    """Multiobjective."""
    m = gurobipy.Model("Multiobjective Flow")
    var = milp_vars.MaxPlasmidScore(
        milp_vars.MaxGC(
            m,
            milp_vars.MaxCovFlow(m, network),
            network,
            intervals,
        ),
    )
    milp_objs.set_multiobjective(m, var, network, intervals)
    milp_consts.set_multiobj_constraints(m, var, network, intervals)
    return m, var
