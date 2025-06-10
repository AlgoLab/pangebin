"""Plasbin classbin MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class Classify:
    """Classify model variables."""

    def __init__(
        self,
        frag: pb_lp_var.SubFragments,
        sub_v: pb_lp_var.SubVertices,
        sub_arc: pb_lp_var.SubArcs,
        flow: pb_lp_var.Flow,
        ccomp: pb_lp_var.ConnectedComponent,
        seed_tree_arcs: pb_lp_var.SeedTreeArcs,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__ccomp = ccomp
        self.__seed_tree_arcs = seed_tree_arcs

    def frag(self) -> pb_lp_var.SubFragments:
        """Get fragment variables."""
        return self.__frag

    def sub_vertices(self) -> pb_lp_var.SubVertices:
        """Get subvertex variables."""
        return self.__sub_v

    def sub_arcs(self) -> pb_lp_var.SubArcs:
        """Get subarc variables."""
        return self.__sub_arc

    def flow(self) -> pb_lp_var.Flow:
        """Get flow variables."""
        return self.__flow

    def ccomp(self) -> pb_lp_var.ConnectedComponent:
        """Get connected component variables."""
        return self.__ccomp

    def seed_tree_arcs(self) -> pb_lp_var.SeedTreeArcs:
        """Get seed tree arc variables."""
        return self.__seed_tree_arcs


def init_classify(network: net.Network, model: gp.Model) -> Classify:
    """Create default Classify variables."""
    frag = pb_lp_var.SubFragments(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    sub_v = pb_lp_var.SubVertices(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    sub_arc = pb_lp_var.SubArcs(
        network,
        model,
        Domain.binary(),
    )
    flow = pb_lp_var.Flow(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    ccomp = pb_lp_var.ConnectedComponent(
        network,
        model,
        Domain.continuous(-gp.GRB.INFINITY, gp.GRB.INFINITY),
        Domain.continuous(-gp.GRB.INFINITY, 0),
    )
    seed_tree_arcs = pb_lp_var.SeedTreeArcs(
        network,
        model,
        Domain.continuous(0, len(network.seeds())),
    )
    return Classify(frag, sub_v, sub_arc, flow, ccomp, seed_tree_arcs)
