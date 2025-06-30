"""Plasbin once MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.connected_component.variables as ccomp_var
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class MaxGCLabelBinScore:
    """MGCLB variables."""

    def __init__(  # noqa: PLR0913
        self,
        frag: cmn_lp_vars.SubFragments,
        sub_v: cmn_lp_vars.SubVertices,
        sub_arc: cmn_lp_vars.SubArcs,
        flow: cmn_lp_vars.Flow,
        pos_flow: cmn_lp_vars.PositiveFlow,
        tree_edges_vars: ccomp_var.TreeEdges,
        root: ccomp_var.Root,
        gc: cmn_lp_vars.GCIntervals,
        inflow_gc: cmn_lp_vars.InflowGC,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__pos_flow = pos_flow
        self.__tree_edges_vars = tree_edges_vars
        self.__root = root
        self.__gc = gc
        self.__inflow_gc = inflow_gc

    def sub_frags(self) -> cmn_lp_vars.SubFragments:
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

    def pos_flow(self) -> cmn_lp_vars.PositiveFlow:
        """Get positive flow variables."""
        return self.__pos_flow

    def gc(self) -> cmn_lp_vars.GCIntervals:
        """Get GC variables."""
        return self.__gc

    def inflow_gc(self) -> cmn_lp_vars.InflowGC:
        """Get inflow GC variables."""
        # REFACTOR inflowgc is completely useless: use inflow instead
        return self.__inflow_gc


def init_mgclb(
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    model: gp.Model,
) -> MaxGCLabelBinScore:
    """Create default MaxGCLabelBinScore variables."""
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
    pos_flow = cmn_lp_vars.PositiveFlow(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    tree_edges_vars = ccomp_var.TreeEdges(
        network,
        model,
        Domain.continuous(0, network.number_of_vertices()),
    )
    root = ccomp_var.Root(network, model, Domain.binary())
    gc = cmn_lp_vars.GCIntervals(
        gc_intervals,
        model,
        Domain.binary(),
    )
    inflow_gc = cmn_lp_vars.InflowGC(
        network,
        gc_intervals,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    return MaxGCLabelBinScore(
        frag,
        sub_v,
        sub_arc,
        flow,
        pos_flow,
        tree_edges_vars,
        root,
        gc,
        inflow_gc,
    )
