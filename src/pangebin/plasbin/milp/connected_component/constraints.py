"""PangeBin-flow MILP constraints library.

The following constraint function serve as base bricks to compose MILP models.
"""

from __future__ import annotations

import gurobipy as gp

import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net

from . import variables as ccomp_var


def arcs_in_tree_are_active(
    m: gp.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    dtree_vars: ccomp_var.TreeArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Arcs in tree are active."""
    constraints: list[gp.Constr] = []
    #
    # Source-arcs
    #
    constraints.extend(
        m.addConstrs(
            (
                dtree_vars.beta_s(source_arc)
                <= sub_arc_vars.s(source_arc) * network.number_of_vertices()
                for source_arc in network.source_arcs()
            ),
            name="source_arcs_in_tree_are_active",
        ).values(),
    )
    #
    # Link-arcs
    #
    constraints.extend(
        m.addConstrs(
            (
                dtree_vars.beta(arc)
                <= sub_arc_vars.l(arc) * network.number_of_vertices()
                for arc in network.link_arcs()
            ),
            name="link_arcs_in_tree_are_active",
        ).values(),
    )
    #
    # Sink-arcs
    #
    constraints.extend(
        m.addConstrs(
            (
                dtree_vars.beta_t(sink_arc)
                <= sub_arc_vars.t(sink_arc) * network.number_of_vertices()
                for sink_arc in network.sink_arcs()
            ),
            name="sink_arcs_in_tree_are_active",
        ).values(),
    )
    return constraints


def subtree_size_from_source(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    dtree_vars: ccomp_var.TreeArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Subtree size from the source.

    Warning
    -------
    The source `s` is assumed to be in the component.
    To describe an empty graph, set the right hand side to `0`.
    Otherwise, set the right hand side to `-1` (default).
    """
    return [
        m.addConstr(
            gp.quicksum(
                sub_vertices_vars.x(oriented_fragment)
                for oriented_fragment in network.oriented_fragments()
            )
            - gp.quicksum(
                dtree_vars.beta_s(source_arcs) for source_arcs in network.source_arcs()
            )
            == -1,
            # The `-1` corresponds to the `s` and `t` not counted in the sum of `x_v`
            # Here we consider the source in the ccomp.
            name="subtree_size_from_source",
        ),
    ]


def subtree_size_but_the_source(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    dtree_vars: ccomp_var.TreeArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Subtree size from all vertices but the source."""
    constraints: list[gp.Constr] = []
    #
    # For the fragments
    #
    constraints.extend(
        m.addConstrs(
            (
                dtree_vars.incoming_beta(network, fragment)
                - dtree_vars.outgoing_beta(network, fragment)
                == sub_vertices_vars.x(fragment)
                for fragment in network.oriented_fragments()
            ),
            name="subtree_size_from_oriented_fragment",
        ).values(),
    )
    #
    # For the sink
    #
    constraints.append(
        m.addConstr(
            gp.quicksum(
                dtree_vars.beta_t(sink_arcs) for sink_arcs in network.sink_arcs()
            )
            == 1,
            name="subtree_size_from_sink",
        ),
    )
    return constraints


# ------------------------------------------------------------------------------------ #
#                              Connected Intermediate Flow                             #
# ------------------------------------------------------------------------------------ #
def rev_link_arcs_in_tree_are_active(
    m: gp.Model,
    sub_arc_vars: pb_lp_var.SubArcs,
    tree_vars: ccomp_var.TreeEdges,
    network: net.Network,
) -> list[gp.Constr]:
    """Arcs in tree are active."""
    #
    # Link-arcs
    #
    return list(
        m.addConstrs(
            (
                tree_vars.beta_rev(arc)
                <= sub_arc_vars.l(arc) * network.number_of_vertices()
                for arc in network.link_arcs()
            ),
            name="link_arcs_in_tree_are_active_rev",
        ).values(),
    )


def subtree_size_fragment_rev_version(
    m: gp.Model,
    sub_vertices_vars: pb_lp_var.SubVertices,
    tree_vars: ccomp_var.TreeEdges,
    network: net.Network,
) -> list[gp.Constr]:
    """Subtree size from all vertices but the source."""
    #
    # For the fragments
    #
    return list(
        m.addConstrs(
            (
                tree_vars.dtree().incoming_beta(network, fragment)
                + tree_vars.outgoing_beta_rev(network, fragment)
                - tree_vars.dtree().outgoing_beta(network, fragment)
                - tree_vars.incoming_beta_rev(network, fragment)
                == sub_vertices_vars.x(fragment)
                for fragment in network.oriented_fragments()
            ),
            name="subtree_size_from_oriented_fragment_rev",
        ).values(),
    )


def subtree_size_for_sink_rev_version(
    m: gp.Model,
    tree_vars: ccomp_var.TreeEdges,
    network: net.Network,
) -> list[gp.Constr]:
    """Subtree size from all vertices but the source.

    Warning
    -------
    The sink `t` is assumed to be in the component.
    To describe an empty graph, set the right hand side to `0`.
    Otherwise, set the right hand side to `1` (default).
    """
    #
    # For the sink
    #
    return [
        m.addConstr(
            gp.quicksum(
                tree_vars.dtree().beta_t(sink_arcs) for sink_arcs in network.sink_arcs()
            )
            == 1,
            name="subtree_size_for_sink_rev_version",
        ),
    ]


def only_one_oriented_fragment_connects_the_source(
    m: gp.Model,
    root_vars: ccomp_var.Root,
    network: net.Network,
) -> gp.Constr:
    """Set only one oriented fragment connects the source."""
    return m.addConstr(
        gp.quicksum(
            root_vars.r(v)
            for v in network.oriented_fragments()
            if network.is_source_connected(v.identifier())
        )
        == 1,
        name="only_one_oriented_fragment_connects_the_source",
    )


def at_most_one_source_arc_in_tree(
    m: gp.Model,
    tree_vars: ccomp_var.TreeEdges,
    root_vars: ccomp_var.Root,
    network: net.Network,
) -> list[gp.Constr]:
    """Set exactly one arc in tree."""
    return list(
        m.addConstrs(
            (
                tree_vars.dtree().beta_s(source_arc)
                <= root_vars.r(source_arc[1]) * network.number_of_vertices()
                for source_arc in network.source_arcs()
            ),
            name="at_most_one_source_arc_in_tree",
        ).values(),
    )


def root_priority(
    m: gp.Model,
    root_vars: ccomp_var.Root,
    sub_arcs_vars: pb_lp_var.SubArcs,
    network: net.Network,
) -> list[gp.Constr]:
    """Set root priority."""
    root_list = sorted(network.oriented_fragments(), key=lambda v: v.identifier())
    return list(
        m.addConstrs(
            (
                root_vars.r(root_list[k + 1])
                <= root_vars.r(root_list[k])
                + (1 - sub_arcs_vars.s((network.SOURCE_VERTEX, root_list[k])))
                for k in range(len(root_list) - 1)
            ),
            name="root_priority",
        ).values(),
    )


def beta_rev_upper_bound(
    m: gp.Model,
    tree_edges_vars: ccomp_var.TreeEdges,
    network: net.Network,
) -> list[gp.Constr]:
    """Set beta_rev upper bound.

    Notes
    -----
    (Default) Refuse the use of beta_rev: set RHS to `0`.
    Allow the use of beta_rev: set RHS to `1`.
    """
    return list(
        m.addConstrs(
            (tree_edges_vars.beta_rev(arc) <= 0 for arc in network.link_arcs()),
            name="beta_rev_upper_bound",
        ).values(),
    )


def non_seed_root_upper_bound(
    m: gp.Model,
    root_vars: ccomp_var.Root,
    network: net.Network,
) -> list[gp.Constr]:
    """Set non-seed root upper bound.

    Notes
    -----
    (Default) Refuse a non-seed to be a root: set RHS to `0`.
    Allow a non-seed to be a root: set RHS to `1`.
    """
    return list(
        m.addConstrs(
            (
                root_vars.r(orient_frag) <= 0
                for frag_id in network.source_connected_fragment_ids()
                if frag_id not in network.seeds()
                for orient_frag in network.to_oriented(frag_id)
            ),
            name="non_seed_root_ub",
        ).values(),
    )
