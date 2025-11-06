"""MILP objectives."""

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.classbin.multi_flow.milp.variables as mfb_var
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net


def frag_coeff(network: net.Network, frag_id: str) -> float:
    """Get fragment coefficient."""
    return net.length(network, frag_id)


def coverage_penalty_all(
    bins_vars: list[mfb_var.BinVariables],
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
            -(
                network.coverage(frag_id) * flow_union_frag_vars.frag(frag_id)
                - gp.quicksum(
                    bin_vars.flows().incoming_forward_reverse(network, frag_id)
                    for bin_vars in bins_vars
                )
            )
        )
        for frag_id in frag_set_fn(network)
        # BUG TMP only for positive plasmidness
        if network.plasmidness(frag_id) >= 0
    )


def bad_plasmidness_chain_penalty(
    bins_vars: list[mfb_var.BinVariables],
    network: net.Network,
) -> gp.LinExpr:
    """Get linear expression for coverage penalty."""
    return gp.quicksum(
        gp.quicksum(
            (
                frag_coeff(network, link_arc.predecessor().identifier())
                + frag_coeff(network, link_arc.successor().identifier())
            )
            * (
                network.plasmidness(link_arc.predecessor().identifier())
                + network.plasmidness(link_arc.successor().identifier())
            )
            * bin_var.flows().l(link_arc)
            for link_arc in network.link_arcs()
            if network.plasmidness(link_arc.predecessor().identifier()) < 0
            and network.plasmidness(link_arc.successor().identifier()) < 0
        )
        for bin_var in bins_vars
    )


def repeat_penalty(
    bins_vars: list[mfb_var.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
) -> gp.LinExpr:
    return -gp.quicksum(
        frag_coeff(network, frag_id)
        * (
            gp.quicksum(bin_var.sub_frag().frag(frag_id) for bin_var in bins_vars)
            - flow_union_frag_vars.frag(frag_id)
        )
        * network.plasmidness(frag_id)
        for frag_id in network.fragment_ids()
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


def circular_flows_plasmidness_score(
    bins_vars: list[mfb_var.BinVariables],
    network: net.Network,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
    number_of_circular_flows: int,
) -> gp.LinExpr:
    """Get circular flows plasmidness score linear expression."""
    return gp.quicksum(
        plasmidness_score(
            network,
            bin_vars.flows(),
            obj_fun_domain,
        )
        for bin_vars in bins_vars[:number_of_circular_flows]
    )


# REFACTOR potentially useless
def gc_score(
    network: net.Network,
    intervals: gc_items.Intervals,
    inflow_gc_vars: cmn_lp_vars.InflowGC,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for GC score."""
    frag_set_fn = cmn_lp_objs.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        frag_coeff(network, frag_id)
        * network.gc_score(frag_id)[b]
        * inflow_gc_vars.x(frag_id, interval)
        for b, interval in enumerate(intervals)
        for frag_id in frag_set_fn(network)
    )


def mfb_objective(
    bins_vars: list[mfb_var.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get MFB objective linear expression.

    Interval of bin close-open
    """
    return (
        coverage_penalty_all(
            bins_vars,
            flow_union_frag_vars,
            network,
            obj_fun_domain,
        )
        + gp.quicksum(
            plasmidness_score(network, var.flows(), obj_fun_domain) for var in bins_vars
        )
        # + bad_plasmidness_chain_penalty(
        #     bins_vars,
        #     network,
        # )  # BUG TMP bad plasmidness chain penalty
        # TODO CONTINUE HERE
        # test constraint no more than 2 source arc and 2 sink arc
        # BUG repeat penalty
        + repeat_penalty(
            bins_vars,
            flow_union_frag_vars,
            network,
        )
    )


def set_objective(
    m: gp.Model,
    bin_vars: list[mfb_var.BinVariables],
    flow_union_frag_vars: cmn_lp_vars.SubFragments,
    network: net.Network,
    intervals: gc_items.Intervals,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> None:
    """Set MGCLB objective."""
    m.setObjective(
        mfb_objective(
            bin_vars,
            flow_union_frag_vars,
            network,
            intervals,
            obj_fun_domain,
        ),
        gp.GRB.MAXIMIZE,
    )
