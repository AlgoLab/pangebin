"""PangeBin-flow MILP objectives."""

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.network as pb_network


def set_mcf_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> None:
    """MCF objective."""
    m.setObjective(milp_vars.coverage_score(network, var), gurobipy.GRB.MAXIMIZE)


def set_mgc_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxGC,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
) -> None:
    """MGC objective."""
    m.setObjective(
        milp_vars.gc_probability_score(network, intervals, var),
        gurobipy.GRB.MAXIMIZE,
    )


def set_mps_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
) -> None:
    """MPS objective."""
    m.setObjective(
        milp_vars.plasmidness_score(network, intervals, var),
        gurobipy.GRB.MAXIMIZE,
    )


def set_mps_prime_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
) -> None:
    """MPS' objective."""
    m.setObjective(
        network.number_of_fragments()
        * milp_vars.coverage_score(network, var.mcf_vars())
        + milp_vars.gc_probability_score(network, intervals, var.mgc_vars()),
        gurobipy.GRB.MAXIMIZE,
    )
