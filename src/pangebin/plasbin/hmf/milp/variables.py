"""Plasbin classbin MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.plasbin.milp.connected_component.variables as ccomp_var
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class BinVariables:
    """Variables for one bin."""

    def __init__(
        self,
        frag: cmn_lp_vars.SubFragments,
        sub_v: cmn_lp_vars.SubVertices,
        sub_arc: cmn_lp_vars.SubArcs,
        flow: cmn_lp_vars.Flow,
        tree_edges_vars: ccomp_var.TreeEdges,
        root: ccomp_var.Root,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__tree_edges_vars = tree_edges_vars
        self.__root = root

    def sub_frag(self) -> cmn_lp_vars.SubFragments:
        """Get fragment variables."""
        return self.__frag

    def sub_vertices(self) -> cmn_lp_vars.SubVertices:
        """Get subvertex variables."""
        return self.__sub_v

    def sub_arcs(self) -> cmn_lp_vars.SubArcs:
        """Get subarc variables."""
        return self.__sub_arc

    def flows(self) -> cmn_lp_vars.Flow:
        """Get flow variables."""
        return self.__flow

    def tree_edges(self) -> ccomp_var.TreeEdges:
        """Get tree edges variables."""
        return self.__tree_edges_vars

    def root(self) -> ccomp_var.Root:
        """Get root variables."""
        return self.__root


def _new_variables_for_one_bin(
    network: net.Network,
    model: gp.Model,
    bin_number: int,
) -> BinVariables:
    """Create default Binning variables."""
    variable_prefix = f"bin_{bin_number}"
    frag = cmn_lp_vars.SubFragments(
        network,
        model,
        Domain.continuous(0, 1),
        variable_prefix,
    )
    sub_v = cmn_lp_vars.SubVertices(
        network,
        model,
        Domain.continuous(0, 1),
        variable_prefix,
    )
    sub_arc = cmn_lp_vars.SubArcs(
        network,
        model,
        Domain.binary(),
        variable_prefix,
    )
    flow = cmn_lp_vars.Flow(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
        Domain.continuous(0, gp.GRB.INFINITY),
        variable_prefix,
    )
    tree_edges_vars = ccomp_var.TreeEdges(
        network,
        model,
        Domain.continuous(0, network.number_of_vertices()),
        variable_prefix,
    )
    root = ccomp_var.Root(network, model, Domain.binary(), variable_prefix)
    return BinVariables(frag, sub_v, sub_arc, flow, tree_edges_vars, root)


def init_binning(
    network: net.Network,
    model: gp.Model,
    max_number_of_bins: int,
) -> tuple[list[BinVariables], cmn_lp_vars.SubFragments]:
    """Create default Binning variables."""
    return (
        [
            _new_variables_for_one_bin(network, model, k)
            for k in range(max_number_of_bins)
        ],
        cmn_lp_vars.SubFragments(
            network,
            model,
            Domain.continuous(0, 1),
            "flow_union",
        ),
    )
