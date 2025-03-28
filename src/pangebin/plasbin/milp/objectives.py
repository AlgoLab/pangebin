"""PangeBin-flow MILP objectives."""

from collections.abc import Callable, Iterable, Iterator
from enum import StrEnum

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.network as pb_network


class ObjectiveFunctionDomain(StrEnum):
    """Objective function domain."""

    ALL = "all"
    SEEDS = "seeds"

    def to_fn(self) -> Callable[[pb_network.Network], Iterator[str] | Iterable[str]]:
        """Get objective function domain function."""
        match self:
            case ObjectiveFunctionDomain.ALL:
                return pb_network.Network.fragment_ids
            case ObjectiveFunctionDomain.SEEDS:
                return pb_network.Network.seeds


def coverage_score(
    network: pb_network.Network,
    var: milp_vars.MaxCovFlow,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for coverage score."""
    # DOCU MCF: Tests on coverage score
    frag_set_fn = ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_length = max(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in frag_set_fn(network)
    )
    return gurobipy.quicksum(
        (gfa_segment.length(network.gfa_graph().segment(frag_id)) / max_frag_length)
        * (
            milp_vars.incoming_flow_forward_reverse(frag_id, network, var)
            - (
                network.coverage(frag_id) * var.frag(frag_id)
                - milp_vars.incoming_flow_forward_reverse(frag_id, network, var)
            )
        )
        for frag_id in frag_set_fn(network)
    )


def set_mcf_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxCovFlow,
    network: pb_network.Network,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> None:
    """MCF objective."""
    m.setObjective(coverage_score(network, var, obj_fun_domain), gurobipy.GRB.MAXIMIZE)


def gc_score(
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    var: milp_vars.MaxGC,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for GC probability score."""
    frag_set_fn = ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gurobipy.quicksum(
        network.gc_score(frag_id)[b] * var.frag_gc(frag_id, interval)
        for b, interval in enumerate(intervals)
        for frag_id in frag_set_fn(network)
    )


def set_mgc_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxGC,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> None:
    """MGC objective."""
    m.setObjective(
        gc_score(network, intervals, var, obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def plasmidness_score(
    network: pb_network.Network,
    var: milp_vars.MaxPlasmidScore,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> gurobipy.LinExpr:
    """Get linear expression for plasmidness score."""
    frag_set_fn = ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    max_frag_len = max(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in frag_set_fn(network)
    )  # DOCU MPS: coeff in obj
    return gurobipy.quicksum(
        (gfa_segment.length(network.gfa_graph().segment(frag_id)) / max_frag_len)
        * network.plasmidness(frag_id)
        * var.mcf_vars().frag(frag_id)
        for frag_id in frag_set_fn(network)
    )


def set_mps_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> None:
    """MPS objective."""
    m.setObjective(
        plasmidness_score(network, var, obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def set_mps_prime_objective(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    obj_fun_domain: ObjectiveFunctionDomain,
    network: pb_network.Network,
) -> None:
    """MPS' objective."""
    # DOCU MPS': max total flow
    m.setObjective(
        coverage_score(network, var.mcf_vars(), obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )


def set_multiobjective(
    m: gurobipy.Model,
    var: milp_vars.MaxPlasmidScore,
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: ObjectiveFunctionDomain,
) -> None:
    """Multiobjective."""
    # TODO multiobj plasbin''
    m.setObjective(
        var.mcf_vars().total_flow()
        + gc_score(network, intervals, var.mgc_vars(), obj_fun_domain)
        + plasmidness_score(network, var, obj_fun_domain),
        gurobipy.GRB.MAXIMIZE,
    )
