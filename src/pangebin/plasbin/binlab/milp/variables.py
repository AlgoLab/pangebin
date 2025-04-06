"""Plasbin binlab MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class MaxBinScore:
    """MBS variables."""

    def __init__(  # noqa: PLR0913
        self,
        frag: pb_lp_var.SubFragments,
        sub_v: pb_lp_var.SubVertices,
        sub_arc: pb_lp_var.SubArcs,
        flow: pb_lp_var.Flow,
        pos_flow: pb_lp_var.PositiveFlow,
        ccomp: pb_lp_var.ConnectedComponent,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__pos_flow = pos_flow
        self.__ccomp = ccomp

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


class MaxLabScore(MaxBinScore):
    """MLS variables."""

    def __init__(  # noqa: PLR0913
        self,
        frag: pb_lp_var.SubFragments,
        sub_v: pb_lp_var.SubVertices,
        sub_arc: pb_lp_var.SubArcs,
        flow: pb_lp_var.Flow,
        pos_flow: pb_lp_var.PositiveFlow,
        ccomp: pb_lp_var.ConnectedComponent,
        gc: pb_lp_var.GCIntervals,
        frag_gc: pb_lp_var.FragmentGC,
    ) -> None:
        """Create Labelling variables."""
        super().__init__(frag, sub_v, sub_arc, flow, pos_flow, ccomp)
        self.__gc = gc
        self.__frag_gc = frag_gc

    def gc(self) -> pb_lp_var.GCIntervals:
        """Get GC variables."""
        return self.__gc

    def frag_gc(self) -> pb_lp_var.FragmentGC:
        """Get fragment GC variables."""
        return self.__frag_gc


class MaxRefBinScore(MaxLabScore):
    """MRBS variables."""


def init_mbs(network: net.Network, model: gp.Model) -> MaxBinScore:
    """Create default Binning variables."""
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
    return MaxBinScore(frag, sub_v, sub_arc, flow, pos_flow, ccomp)


def mls_from_mbs(
    model: gp.Model,
    mbs_vars: MaxBinScore,
    network: net.Network,
    gc_intervals: gc_items.Intervals,
) -> MaxLabScore:
    """Create MLS variables from MBS variables."""
    gc = pb_lp_var.GCIntervals(
        gc_intervals,
        model,
        Domain.binary(),
    )
    frag_gc = pb_lp_var.FragmentGC(
        network,
        gc_intervals,
        model,
        Domain.continuous(0, gp.GRB.INFINITY),
    )
    return MaxLabScore(
        mbs_vars.frag(),
        mbs_vars.sub_v(),
        mbs_vars.sub_arc(),
        mbs_vars.flow(),
        mbs_vars.pos_flow(),
        mbs_vars.ccomp(),
        gc,
        frag_gc,
    )


def mrbs_from_mls(mls_vars: MaxLabScore) -> MaxRefBinScore:
    """Create MRBS variables from MLS variables."""
    return MaxRefBinScore(
        mls_vars.frag(),
        mls_vars.sub_v(),
        mls_vars.sub_arc(),
        mls_vars.flow(),
        mls_vars.pos_flow(),
        mls_vars.ccomp(),
        mls_vars.gc(),
        mls_vars.frag_gc(),
    )
