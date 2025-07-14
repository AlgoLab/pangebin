"""MILP objectives."""

import gurobipy as gp

import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net

from . import variables as lp_vars


def frag_coeff(network: net.Network, frag_id: str) -> float:
    """Get fragment coefficient."""
    return net.length(network, frag_id)


def active_coverage_penalty(
    bins_vars: list[lp_vars.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage penalty."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        frag_coeff(network, frag_id)
        # BUG TMP PENALTY MULTIPLIED BY PLASMIDNESS
        * network.plasmidness(frag_id)
        * (
            gp.quicksum(
                bin_vars.flows().incoming_forward_reverse(network, frag_id)
                for bin_vars in bins_vars
            )
            - network.coverage(frag_id) * flow_union_frag_vars.frag(frag_id)
        )
        for frag_id in frag_set_fn(network)
        # BUG TMP only for positive plasmidness
        if network.plasmidness(frag_id) > 0
    )


def all_coverage_penalty(
    bins_vars: list[lp_vars.BinVariables],
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get coverage penalty for plasmidness positive fragments."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        frag_coeff(network, frag_id)
        # BUG TMP PENALTY MULTIPLIED BY PLASMIDNESS
        * network.plasmidness(frag_id)
        * (
            gp.quicksum(
                bin_vars.flows().incoming_forward_reverse(network, frag_id)
                for bin_vars in bins_vars
            )
            - network.coverage(frag_id)
        )
        for frag_id in frag_set_fn(network)
        # BUG TMP only for positive plasmidness
        if network.plasmidness(frag_id) > 0
    )


def repeat_penalty(
    bins_vars: list[lp_vars.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
) -> gp.LinExpr:
    """Get linear expression for repeat penalty."""
    return gp.quicksum(
        frag_coeff(network, frag_id)
        * network.plasmidness(frag_id)
        * (
            flow_union_frag_vars.frag(frag_id)
            - gp.quicksum(bin_var.sub_frag().frag(frag_id) for bin_var in bins_vars)
        )
        # BUG repeat penalty coeff on coverage
        * (network.coverage(frag_id) / 3)
        for frag_id in network.fragment_ids()
        # BUG TMP only for positive plasmidness
        if network.plasmidness(frag_id) > 0
    )


def plasmidness_score(
    network: net.Network,
    flow_vars: cmn_lp_vars.Flow,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        frag_coeff(network, frag_id)
        * network.plasmidness(frag_id)
        * flow_vars.incoming_forward_reverse(network, frag_id)
        for frag_id in frag_set_fn(network)
    )


def multi_flow_plasmidness_score(
    bins_vars: list[lp_vars.BinVariables],
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
    from_until: tuple[int, int],
) -> gp.LinExpr:
    """Get plasmidness score linear expression for a multi-flow interval."""
    return gp.quicksum(
        plasmidness_score(
            network,
            bin_vars.flows(),
            obj_fun_domain,
        )
        for bin_vars in bins_vars[from_until[0] : from_until[1]]
    )


def circular_objective(
    bins_vars: list[lp_vars.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get MFB circular objective linear expression."""
    return (
        # BUG TMP
        # active_coverage_penalty(
        #     bins_vars,
        #     flow_union_frag_vars,
        #     network,
        #     obj_fun_domain,
        # )
        # all_coverage_penalty(
        #     bins_vars,
        #     network,
        #     obj_fun_domain,
        # )
        0  # BUG TMP REMOVE COV PENALTY
        + gp.quicksum(
            plasmidness_score(network, var.flows(), obj_fun_domain) for var in bins_vars
        )
        # BUG repeat penalty
        + repeat_penalty(
            bins_vars,
            flow_union_frag_vars,
            network,
        )
    )


def partially_circular_objective(
    bins_vars: list[lp_vars.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get MFB partially circular objective linear expression."""
    return (
        all_coverage_penalty(
            bins_vars,
            network,
            obj_fun_domain,
        )
        + gp.quicksum(
            plasmidness_score(network, var.flows(), obj_fun_domain) for var in bins_vars
        )
        # BUG repeat penalty
        + repeat_penalty(
            bins_vars,
            flow_union_frag_vars,
            network,
        )
    )
