"""Plasbin once MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class MaxGCLabelBinScore:
    """MGCLB variables."""

    def __init__(  # noqa: PLR0913
        self,
        frag: pb_lp_var.SubFragments,
        sub_v: pb_lp_var.SubVertices,
        sub_arc: pb_lp_var.SubArcs,
        flow: pb_lp_var.Flow,
        pos_flow: pb_lp_var.PositiveFlow,
        ccomp: pb_lp_var.ConnectedComponent,
        gc: pb_lp_var.GCIntervals,
        inflow_gc: pb_lp_var.InflowGC,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__pos_flow = pos_flow
        self.__ccomp = ccomp
        self.__gc = gc
        self.__inflow_gc = inflow_gc

    def frag(self) -> pb_lp_var.SubFragments:
        """Get fragment variables."""
        return self.__frag

    def sub_v(self) -> pb_lp_var.SubVertices:
        """Get subvertex variables."""
        return self.__sub_v

    def sub_arc(self) -> pb_lp_var.SubArcs:
        """Get subarc variables."""
        return self.__sub_arc

    def flow(self) -> pb_lp_var.Flow:
        """Get flow variables."""
        return self.__flow

    def pos_flow(self) -> pb_lp_var.PositiveFlow:
        """Get positive flow variables."""
        return self.__pos_flow

    def ccomp(self) -> pb_lp_var.ConnectedComponent:
        """Get connected component variables."""
        return self.__ccomp

    def gc(self) -> pb_lp_var.GCIntervals:
        """Get GC variables."""
        return self.__gc

    def inflow_gc(self) -> pb_lp_var.InflowGC:
        """Get inflow GC variables."""
        return self.__inflow_gc


def init_mgclb(
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    model: gp.Model,
) -> MaxGCLabelBinScore:
    """Create default MaxGCLabelBinScore variables."""
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
    pos_flow = pb_lp_var.PositiveFlow(
        network,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    ccomp = pb_lp_var.ConnectedComponent(
        network,
        model,
        Domain.continuous(-gp.GRB.INFINITY, gp.GRB.INFINITY),
        Domain.continuous(-gp.GRB.INFINITY, 0),
    )
    gc = pb_lp_var.GCIntervals(
        gc_intervals,
        model,
        Domain.binary(),
    )
    inflow_gc = pb_lp_var.InflowGC(
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
        ccomp,
        gc,
        inflow_gc,
    )
