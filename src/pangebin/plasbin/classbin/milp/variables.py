"""Plasbin classbin MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class Classify:
    """Classify model variables."""

    def __init__(
        self,
        frag: cmn_lp_vars.SubFragments,
        sub_v: cmn_lp_vars.SubVertices,
        sub_arc: cmn_lp_vars.SubArcs,
        flow: cmn_lp_vars.Flow,
        ccomp: cmn_lp_vars.ConnectedComponent,
        seed_tree_arcs: cmn_lp_vars.SeedTreeArcs,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__ccomp = ccomp
        self.__seed_tree_arcs = seed_tree_arcs

    def frag(self) -> cmn_lp_vars.SubFragments:
        """Get fragment variables."""
        return self.__frag

    def sub_vertices(self) -> cmn_lp_vars.SubVertices:
        """Get subvertex variables."""
        return self.__sub_v

    def sub_arcs(self) -> cmn_lp_vars.SubArcs:
        """Get subarc variables."""
        return self.__sub_arc

    def flow(self) -> cmn_lp_vars.Flow:
        """Get flow variables."""
        return self.__flow

    def ccomp(self) -> cmn_lp_vars.ConnectedComponent:
        """Get connected component variables."""
        return self.__ccomp

    def seed_tree_arcs(self) -> cmn_lp_vars.SeedTreeArcs:
        """Get seed tree arc variables."""
        return self.__seed_tree_arcs


def init_classify(network: net.Network, model: gp.Model) -> Classify:
    """Create default Classify variables."""
    frag = cmn_lp_vars.SubFragments(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    sub_v = cmn_lp_vars.SubVertices(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    sub_arc = cmn_lp_vars.SubArcs(
        network,
        model,
        Domain.binary(),
    )
    flow = cmn_lp_vars.Flow(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    ccomp = cmn_lp_vars.ConnectedComponent(
        network,
        model,
        Domain.continuous(-gp.GRB.INFINITY, gp.GRB.INFINITY),
        Domain.continuous(-gp.GRB.INFINITY, 0),
    )
    seed_tree_arcs = cmn_lp_vars.SeedTreeArcs(
        network,
        model,
        Domain.continuous(0, len(network.seeds())),
    )
    return Classify(frag, sub_v, sub_arc, flow, ccomp, seed_tree_arcs)
