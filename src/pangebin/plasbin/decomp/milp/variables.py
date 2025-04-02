"""Hiearchical decomposition variables."""

from __future__ import annotations

import gurobipy as gp

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain

# TODO use new variable frag in MGC and MPS models


class MaxCovFlow:
    """Max coverage flow variables."""

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


class MaxGC(MaxCovFlow):
    """MGC variables."""

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
        """Create MaxGC variables."""
        super().__init__(frag, sub_v, sub_arc, flow, pos_flow, ccomp)
        self.__gc = gc
        self.__frag_gc = frag_gc

    def gc(self) -> pb_lp_var.GCIntervals:
        """Get GC variables."""
        return self.__gc

    def frag_gc(self) -> pb_lp_var.FragmentGC:
        """Get fragment GC variables."""
        return self.__frag_gc


class MaxPlasmidScore(MaxGC):
    """MPS variables."""


class MaxRefCovFlow(MaxPlasmidScore):
    """MRCF variables."""


def init_mcf(network: net.Network, model: gp.Model) -> MaxCovFlow:
    """Create MCF variables."""
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
    return MaxCovFlow(frag, sub_v, sub_arc, flow, pos_flow, ccomp)


def mgc_from_mcf(
    model: gp.Model,
    mcf_vars: MaxCovFlow,
    network: net.Network,
    gc_intervals: gc_items.Intervals,
) -> MaxGC:
    """Create MGC variables from MCF variables."""
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
    return MaxGC(
        mcf_vars.frag(),
        mcf_vars.sub_v(),
        mcf_vars.sub_arc(),
        mcf_vars.flow(),
        mcf_vars.pos_flow(),
        mcf_vars.ccomp(),
        gc,
        frag_gc,
    )


def mps_from_mgc(mgc_vars: MaxGC) -> MaxPlasmidScore:
    """Create MPS variables from MGC variables."""
    return MaxPlasmidScore(
        mgc_vars.frag(),
        mgc_vars.sub_v(),
        mgc_vars.sub_arc(),
        mgc_vars.flow(),
        mgc_vars.pos_flow(),
        mgc_vars.ccomp(),
        mgc_vars.gc(),
        mgc_vars.frag_gc(),
    )


def mrcf_from_mps(mps_vars: MaxPlasmidScore) -> MaxRefCovFlow:
    """Create MRCF variables from MPS variables."""
    return MaxRefCovFlow(
        mps_vars.frag(),
        mps_vars.sub_v(),
        mps_vars.sub_arc(),
        mps_vars.flow(),
        mps_vars.pos_flow(),
        mps_vars.ccomp(),
        mps_vars.gc(),
        mps_vars.frag_gc(),
    )
