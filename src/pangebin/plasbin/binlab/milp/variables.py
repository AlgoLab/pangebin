"""Plasbin binlab MILP variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.variables as cmn_lp_vars
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain


class MaxBinScore:
    """MBS variables."""

    def __init__(  # noqa: PLR0913
        self,
        frag: cmn_lp_vars.SubFragments,
        sub_v: cmn_lp_vars.SubVertices,
        sub_arc: cmn_lp_vars.SubArcs,
        flow: cmn_lp_vars.Flow,
        pos_flow: cmn_lp_vars.PositiveFlow,
        ccomp: cmn_lp_vars.ConnectedComponent,
    ) -> None:
        self.__frag = frag
        self.__sub_v = sub_v
        self.__sub_arc = sub_arc
        self.__flow = flow
        self.__pos_flow = pos_flow
        self.__ccomp = ccomp

    def frag(self) -> cmn_lp_vars.SubFragments:
        """Get fragment variables."""
        return self.__frag

    def sub_v(self) -> cmn_lp_vars.SubVertices:
        """Get subvertex variables."""
        return self.__sub_v

    def sub_arc(self) -> cmn_lp_vars.SubArcs:
        """Get subarc variables."""
        return self.__sub_arc

    def flow(self) -> cmn_lp_vars.Flow:
        """Get flow variables."""
        return self.__flow

    def pos_flow(self) -> cmn_lp_vars.PositiveFlow:
        """Get positive flow variables."""
        return self.__pos_flow

    def ccomp(self) -> cmn_lp_vars.ConnectedComponent:
        """Get connected component variables."""
        return self.__ccomp


class MaxLabScore(MaxBinScore):
    """MLS variables."""

    def __init__(  # noqa: PLR0913
        self,
        frag: cmn_lp_vars.SubFragments,
        sub_v: cmn_lp_vars.SubVertices,
        sub_arc: cmn_lp_vars.SubArcs,
        flow: cmn_lp_vars.Flow,
        pos_flow: cmn_lp_vars.PositiveFlow,
        ccomp: cmn_lp_vars.ConnectedComponent,
        gc: cmn_lp_vars.GCIntervals,
        frag_gc: cmn_lp_vars.FragmentGC,
    ) -> None:
        """Create Labelling variables."""
        super().__init__(frag, sub_v, sub_arc, flow, pos_flow, ccomp)
        self.__gc = gc
        self.__frag_gc = frag_gc

    def gc(self) -> cmn_lp_vars.GCIntervals:
        """Get GC variables."""
        return self.__gc

    def frag_gc(self) -> cmn_lp_vars.FragmentGC:
        """Get fragment GC variables."""
        return self.__frag_gc


class MaxRefBinScore(MaxLabScore):
    """MRBS variables."""


def init_mbs(network: net.Network, model: gp.Model) -> MaxBinScore:
    """Create default Binning variables."""
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
    ccomp = cmn_lp_vars.ConnectedComponent(
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
    gc = cmn_lp_vars.GCIntervals(
        gc_intervals,
        model,
        Domain.binary(),
    )
    frag_gc = cmn_lp_vars.FragmentGC(
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
