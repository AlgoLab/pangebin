"""PangeBin-flow MILP variables."""

from __future__ import annotations

import logging
from itertools import chain
from typing import TYPE_CHECKING

import gurobipy

import pangebin.gc_content.items as gc_items
import pangebin.gfa.link as gfa_link
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.network as pb_network

if TYPE_CHECKING:
    from collections.abc import Iterator

_LOGGER = logging.getLogger(__name__)


class MaxCovFlow:
    """Max coverage flow variables."""

    def __init__(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> None:
        self.__x: dict[str, gurobipy.Var] = dict(self.__init_x(model, network))
        self.__y: dict[str, gurobipy.Var] = dict(self.__init_y(model, network))
        self.__f: dict[str, gurobipy.Var] = dict(self.__init_f(model, network))
        self.__total_flow = model.addVar(
            name="F",
            vtype=gurobipy.GRB.CONTINUOUS,
            lb=0.0,
        )
        # F_a for all arcs a
        self.__pos_f: dict[str, gurobipy.Var] = dict(
            self.__init_pos_f(model, network),
        )
        # alpha, be the number of vertices in the connected component
        self.__alpha: gurobipy.Var = model.addVar(
            name="alpha",
            vtype=gurobipy.GRB.CONTINUOUS,
        )
        self.__beta: dict[str, gurobipy.Var] = dict(
            self.__init_beta(model, network),
        )
        # TODO use new variable frag in MGC and MPS models
        self.__frag: dict[str, gurobipy.Var] = dict(
            self.__init_frag(model, network),
        )

    def x(self, fragment: gfa_segment.OrientedFragment) -> gurobipy.Var:
        """Get x_i variable, i in Fragments."""
        return self.__x[str(fragment)]

    def y_s(self, source_arc: tuple[str, gfa_segment.OrientedFragment]) -> gurobipy.Var:
        """Get y_a variable, a in source Arcs."""
        return self.__y[self.__fmt_s_arcs(source_arc)]

    def y(self, arc: gfa_link.Link) -> gurobipy.Var:
        """Get y_a variable, a in Arcs."""
        return self.__y[str(arc)]

    def y_t(self, sink_arc: tuple[gfa_segment.OrientedFragment, str]) -> gurobipy.Var:
        """Get y_a variable, a in sink Arcs."""
        return self.__y[self.__fmt_t_arcs(sink_arc)]

    def f_s(self, source_arc: tuple[str, gfa_segment.OrientedFragment]) -> gurobipy.Var:
        """Get f_a variable, a in source Arcs."""
        return self.__f[self.__fmt_s_arcs(source_arc)]

    def f(self, arc: gfa_link.Link) -> gurobipy.Var:
        """Get f_a variable, a in Arcs."""
        return self.__f[str(arc)]

    def f_t(self, sink_arc: tuple[gfa_segment.OrientedFragment, str]) -> gurobipy.Var:
        """Get f_a variable, a in sink Arcs."""
        return self.__f[self.__fmt_t_arcs(sink_arc)]

    def total_flow(self) -> gurobipy.Var:
        """Get F variable."""
        return self.__total_flow

    def pos_f_s(
        self,
        source_arc: tuple[str, gfa_segment.OrientedFragment],
    ) -> gurobipy.Var:
        """Get F_a variable, a in source Arcs."""
        return self.__pos_f[self.__fmt_s_arcs(source_arc)]

    def pos_f(self, arc: gfa_link.Link) -> gurobipy.Var:
        """Get F_a variable, a in Arcs."""
        return self.__pos_f[str(arc)]

    def pos_f_t(
        self,
        sink_arc: tuple[gfa_segment.OrientedFragment, str],
    ) -> gurobipy.Var:
        """Get F_a variable, a in sink Arcs."""
        return self.__pos_f[self.__fmt_t_arcs(sink_arc)]

    def alpha(self) -> gurobipy.Var:
        """Get alpha variable."""
        return self.__alpha

    def beta_s(
        self,
        source_arc: tuple[str, gfa_segment.OrientedFragment],
    ) -> gurobipy.Var:
        """Get beta_a variable, a in source Arcs."""
        return self.__beta[self.__fmt_s_arcs(source_arc)]

    def beta(self, arc: gfa_link.Link) -> gurobipy.Var:
        """Get beta_a variable, a in Arcs."""
        return self.__beta[str(arc)]

    def beta_t(
        self,
        sink_arc: tuple[gfa_segment.OrientedFragment, str],
    ) -> gurobipy.Var:
        """Get beta_a variable, a in sink Arcs."""
        return self.__beta[self.__fmt_t_arcs(sink_arc)]

    def frag(self, frag_id: str) -> gurobipy.Var:
        """Get frag variable."""
        return self.__frag[frag_id]

    def start_with_previous_values(self, variables: MaxCovFlow) -> None:
        """Set start variables values with previous result."""
        for key, var in self.__x.items():
            var.Start = variables.__x[key].X
        for key, var in self.__y.items():
            var.Start = variables.__y[key].X
        for key, var in self.__f.items():
            var.Start = variables.__f[key].X
        self.__total_flow.Start = variables.__total_flow.X
        for key, var in self.__pos_f.items():
            var.Start = variables.__pos_f[key].X
        self.__alpha.Start = variables.__alpha.X
        for key, var in self.__beta.items():
            var.Start = variables.__beta[key].X

    def __init_x(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                str(oriented_fragment),
                model.addVar(
                    name=f"x_{oriented_fragment!s}",
                    vtype=gurobipy.GRB.CONTINUOUS,
                    lb=0.0,
                    ub=gurobipy.GRB.INFINITY,
                ),
            )
            for oriented_fragment in network.oriented_fragments()
        )

    def __init_y(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (link_id, model.addVar(name=f"y_{link_id}", vtype=gurobipy.GRB.BINARY))
            for link_id in self.__network_arc_ids(network)
        )

    def __init_f(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                link_id,
                model.addVar(
                    name=f"f_{link_id}",
                    vtype=gurobipy.GRB.CONTINUOUS,
                    lb=0.0,
                    ub=gurobipy.GRB.INFINITY,
                ),
            )
            for link_id in self.__network_arc_ids(network)
        )

    def __init_pos_f(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                link_id,
                model.addVar(
                    name=f"F_{link_id}",
                    vtype=gurobipy.GRB.CONTINUOUS,
                    lb=0.0,
                    ub=gurobipy.GRB.INFINITY,
                ),
            )
            for link_id in self.__network_arc_ids(network)
        )

    def __init_beta(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                link_id,
                model.addVar(
                    name=f"beta_{link_id}",
                    vtype=gurobipy.GRB.CONTINUOUS,
                    lb=-gurobipy.GRB.INFINITY,
                    ub=0.0,
                ),
            )
            for link_id in self.__network_arc_ids(network)
        )

    def __init_frag(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                frag_id,
                model.addVar(
                    name=f"frag_{frag_id}",
                    vtype=gurobipy.GRB.CONTINUOUS,
                    ub=1.0,
                ),
            )
            for frag_id in network.fragment_ids()
        )

    def __network_arc_ids(
        self,
        network: pb_network.Network,
    ) -> Iterator[str]:
        return chain(
            (self.__fmt_s_arcs(s_link) for s_link in network.source_arcs()),
            (self.__fmt_t_arcs(t_link) for t_link in network.sink_arcs()),
            (str(arc) for arc in network.link_arcs()),
        )

    def __fmt_s_arcs(self, link: tuple[str, gfa_segment.OrientedFragment]) -> str:
        return f"{link[0]}_{link[1]}"

    def __fmt_t_arcs(self, link: tuple[gfa_segment.OrientedFragment, str]) -> str:
        return f"{link[0]}_{link[1]}"


def incoming_flow(
    fragment: gfa_segment.OrientedFragment,
    network: pb_network.Network,
    variables: MaxCovFlow,
) -> gurobipy.LinExpr:
    """Get linear expression for incoming flow."""
    in_flow_arcs = chain(
        (
            variables.f(link)
            for link in gfa_link.incoming_links(network.gfa_graph(), fragment)
        ),
    )
    if network.is_source_connected(fragment.identifier()):
        in_flow_arcs = chain(
            in_flow_arcs,
            iter((variables.f_s((network.SOURCE_VERTEX, fragment)),)),
        )
    return gurobipy.quicksum(in_flow_arcs)


def outgoing_flow(
    fragment: gfa_segment.OrientedFragment,
    network: pb_network.Network,
    variables: MaxCovFlow,
) -> gurobipy.LinExpr:
    """Get linear expression for outgoing flow."""
    out_flow_arcs = chain(
        (
            variables.f(link)
            for link in gfa_link.outgoing_links(network.gfa_graph(), fragment)
        ),
    )
    if network.is_sink_connected(fragment.identifier()):
        out_flow_arcs = chain(
            out_flow_arcs,
            iter((variables.f_t((fragment, network.SINK_VERTEX)),)),
        )
    return gurobipy.quicksum(out_flow_arcs)


def incoming_arcs_y(
    fragment: gfa_segment.OrientedFragment,
    network: pb_network.Network,
    variables: MaxCovFlow,
) -> gurobipy.LinExpr:
    """Get linear expression for incoming flow."""
    in_y_arcs = chain(
        (
            variables.y(link)
            for link in gfa_link.incoming_links(network.gfa_graph(), fragment)
        ),
    )
    if network.is_source_connected(fragment.identifier()):
        in_y_arcs = chain(
            in_y_arcs,
            iter((variables.y_s((network.SOURCE_VERTEX, fragment)),)),
        )
    return gurobipy.quicksum(in_y_arcs)


def incoming_beta(
    fragment: gfa_segment.OrientedFragment,
    network: pb_network.Network,
    variables: MaxCovFlow,
) -> gurobipy.LinExpr:
    """Get linear expression for incoming beta."""
    in_beta_arcs = chain(
        (
            variables.beta(link)
            for link in gfa_link.incoming_links(network.gfa_graph(), fragment)
        ),
    )
    if network.is_source_connected(fragment.identifier()):
        in_beta_arcs = chain(
            in_beta_arcs,
            iter((variables.beta_s((network.SOURCE_VERTEX, fragment)),)),
        )
    return gurobipy.quicksum(in_beta_arcs)


def outgoing_beta(
    fragment: gfa_segment.OrientedFragment,
    network: pb_network.Network,
    variables: MaxCovFlow,
) -> gurobipy.LinExpr:
    """Get linear expression for outgoing beta."""
    out_beta_arcs = chain(
        (
            variables.beta(link)
            for link in gfa_link.outgoing_links(network.gfa_graph(), fragment)
        ),
    )
    if network.is_sink_connected(fragment.identifier()):
        out_beta_arcs = chain(
            out_beta_arcs,
            iter((variables.beta_t((fragment, network.SINK_VERTEX)),)),
        )
    return gurobipy.quicksum(out_beta_arcs)


def incoming_flow_forward_reverse(
    frag_id: str,
    network: pb_network.Network,
    variables: MaxCovFlow,
) -> gurobipy.LinExpr:
    """Get linear expression for incoming flow both of the forward and the reverse."""
    frag_f = gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.FORWARD)
    frag_r = gfa_segment.OrientedFragment(frag_id, gfa_segment.Orientation.REVERSE)
    return incoming_flow(frag_f, network, variables) + incoming_flow(
        frag_r,
        network,
        variables,
    )


def coverage_score(network: pb_network.Network, var: MaxCovFlow) -> gurobipy.LinExpr:
    """Get linear expression for coverage score."""
    # HACK MCF: Tests on coverage score
    # max_seed_cov = max(network.coverage(frag_id) for frag_id in network.seeds())
    # return gurobipy.quicksum(
    #     (network.coverage(frag_id) / max_seed_cov)
    #     * 2
    #     * (
    #         incoming_flow_forward_reverse(frag_id, network, var)
    #         / network.coverage(frag_id)
    #         - 0.5
    #     )
    #     for frag_id in network.fragment_ids()
    #     if frag_id in network.seeds()
    # )
    # DOCU super rapide

    max_frag_len = max(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in network.fragment_ids()
        if frag_id in network.seeds()
    )
    return gurobipy.quicksum(
        (gfa_segment.length(network.gfa_graph().segment(frag_id)) / max_frag_len)
        * incoming_flow_forward_reverse(frag_id, network, var)
        - network.coverage(frag_id) * var.frag(frag_id)
        # -math.exp(
        #     (med_frag_len - gfa_segment.length(network.gfa_graph().segment(frag_id)))
        #     / med_frag_len,
        # )
        # * (
        #     network.coverage(frag_id) * var.frag(frag_id)
        #     - incoming_flow_forward_reverse(frag_id, network, var)
        # )
        for frag_id in network.fragment_ids()
        if frag_id in network.seeds()
    )


class MaxGC:
    """MGC variables."""

    def __init__(
        self,
        model: gurobipy.Model,
        mcf_vars: MaxCovFlow,
        network: pb_network.Network,
        gc_intervals: gc_items.Intervals,
    ) -> None:
        """Create MaxGC variables."""
        self.__mcf_vars = mcf_vars
        self.__gc = dict(self.__init_gc(model, gc_intervals))
        self.__frag_gc = dict(self.__init_frag_gc(model, network, gc_intervals))

    def mcf_vars(self) -> MaxCovFlow:
        """Get MCF variables."""
        return self.__mcf_vars

    def gc(self, interval: tuple[float, float]) -> gurobipy.Var:
        """Get GC variable."""
        return self.__gc[self.__fmt_interval(interval)]

    def frag_gc(self, frag_id: str, interval: tuple[float, float]) -> gurobipy.Var:
        """Get fragment GC variable."""
        return self.__frag_gc[self.__fmt_frag_gc(frag_id, interval)]

    def start_with_previous_values(self, variables: MaxGC) -> None:
        """Set start variables values with previous result."""
        variables.mcf_vars().start_with_previous_values(self.mcf_vars())
        for key, var in self.__gc.items():
            var.Start = variables.__gc[key].X
        for key, var in self.__frag_gc.items():
            var.Start = variables.__frag_gc[key].X

    def __init_gc(
        self,
        model: gurobipy.Model,
        gc_intervals: gc_items.Intervals,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                self.__fmt_interval(interval),
                model.addVar(
                    name=f"gc_{self.__fmt_interval(interval)}",
                    vtype=gurobipy.GRB.BINARY,
                ),
            )
            for interval in gc_intervals
        )

    def __init_frag_gc(
        self,
        model: gurobipy.Model,
        network: pb_network.Network,
        gc_intervals: gc_items.Intervals,
    ) -> Iterator[tuple[str, gurobipy.Var]]:
        return (
            (
                self.__fmt_frag_gc(frag_id, interval),
                model.addVar(
                    name=f"frag_gc_{self.__fmt_frag_gc(frag_id, interval)}",
                    vtype=gurobipy.GRB.CONTINUOUS,
                    lb=0.0,
                    ub=gurobipy.GRB.INFINITY,
                ),
            )
            for frag_id in network.fragment_ids()
            for interval in gc_intervals
        )

    def __fmt_interval(self, interval: tuple[float, float]) -> str:
        return f"{interval[0]}_{interval[1]}"

    def __fmt_frag_gc(self, frag_id: str, interval: tuple[float, float]) -> str:
        return f"{frag_id}_{self.__fmt_interval(interval)}"


def gc_probability_score(
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    var: MaxGC,
) -> gurobipy.LinExpr:
    """Get linear expression for GC probability score."""
    return gurobipy.quicksum(
        network.gc_score(frag_id)[b] * var.frag_gc(frag_id, interval)
        for b, interval in enumerate(intervals)
        for frag_id in network.seeds()
        if gfa_segment.length(
            gfa_segment.get_segment_line_by_name(network.gfa_graph(), frag_id),
        )
        > 1000  # HACK MGC: obj modifications (seeds and parameter)
    )


class MaxPlasmidScore:
    """Max plasmid score variables."""

    def __init__(self, mgc_vars: MaxGC) -> None:
        """Create MaxPlasmidScore variables."""
        self.__mgc_vars = mgc_vars

    def mcf_vars(self) -> MaxCovFlow:
        """Get MCF variables."""
        return self.__mgc_vars.mcf_vars()

    def mgc_vars(self) -> MaxGC:
        """Get MGC variables."""
        return self.__mgc_vars

    def start_with_previous_values(self, variables: MaxPlasmidScore) -> None:
        """Set start variables values with previous result."""
        self.mcf_vars().start_with_previous_values(variables.mcf_vars())
        self.mgc_vars().start_with_previous_values(variables.mgc_vars())


def fragment_gc_sum(
    frag_id: str,
    intervals: gc_items.Intervals,
    var: MaxPlasmidScore,
) -> gurobipy.LinExpr:
    """Get linear expression for fragment GC sum."""
    mgc_vars = var.mgc_vars()
    return gurobipy.quicksum(
        mgc_vars.frag_gc(frag_id, interval) for interval in intervals
    )


def plasmidness_score(
    network: pb_network.Network,
    intervals: gc_items.Intervals,
    var: MaxPlasmidScore,
) -> gurobipy.LinExpr:
    """Get linear expression for plasmidness score."""
    max_frag_len = max(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in network.fragment_ids()
        if frag_id in network.seeds()
    )  # HACK MPS: coeff in obj
    return gurobipy.quicksum(
        (gfa_segment.length(network.gfa_graph().segment(frag_id)) / max_frag_len)
        * network.plasmidness(frag_id)
        * var.mcf_vars().frag(frag_id)
        for frag_id in network.fragment_ids()
        if frag_id in network.seeds()  # HACK MPS: only seeds in plm score
    )
