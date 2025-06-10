"""PangeBin-flow MILP variables library.

Here is base bricks to describe composition of variables.
"""

from __future__ import annotations

from itertools import chain

import gurobipy as gp

import pangebin.gfa.link as gfa_link
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.network as net
from pangebin.plasbin.milp.variables import Domain, gen_vars


class TreeArcs:
    """Directed tree variables.

    Defines:

    * `beta_uv` is the number of vertices in tree defined by root `v`
      for all `(u, v)` in link, source and sink arcs.
    """

    def __init__(
        self,
        network: net.Network,
        model: gp.Model,
        beta_domain: Domain,
        prefix_name: str | None = None,
    ) -> None:
        if prefix_name is None:
            prefix_name = ""
        else:
            prefix_name += "_"
        self.__beta: dict[str, gp.Var] = dict(
            gen_vars(
                model,
                beta_domain,
                f"{prefix_name}beta",
                net.StrFormatter.arc_ids(network),
            ),
        )

    def beta_s(
        self,
        source_arc: tuple[str, gfa_segment.OrientedFragment],
    ) -> gp.Var:
        """Get beta_a variable, a in source-arcs."""
        return self.__beta[net.StrFormatter.s_arc(source_arc)]

    def beta(self, arc: gfa_link.Link) -> gp.Var:
        """Get beta_a variable, a in link-arcs."""
        return self.__beta[net.StrFormatter.arc(arc)]

    def beta_t(
        self,
        sink_arc: tuple[gfa_segment.OrientedFragment, str],
    ) -> gp.Var:
        """Get beta_a variable, a in sink-arcs."""
        return self.__beta[net.StrFormatter.t_arc(sink_arc)]

    def incoming_beta(
        self,
        network: net.Network,
        fragment: gfa_segment.OrientedFragment,
    ) -> gp.LinExpr:
        """Get linear expression for incoming beta."""
        in_beta_arcs = chain(
            (
                self.beta(link)
                for link in gfa_link.incoming_links(network.gfa_graph(), fragment)
            ),
        )
        if network.is_source_connected(fragment.identifier()):
            in_beta_arcs = chain(
                in_beta_arcs,
                iter((self.beta_s((network.SOURCE_VERTEX, fragment)),)),
            )
        return gp.quicksum(in_beta_arcs)

    def outgoing_beta(
        self,
        network: net.Network,
        fragment: gfa_segment.OrientedFragment,
    ) -> gp.LinExpr:
        """Get linear expression for outgoing beta."""
        out_beta_arcs = chain(
            (
                self.beta(link)
                for link in gfa_link.outgoing_links(network.gfa_graph(), fragment)
            ),
        )
        if network.is_sink_connected(fragment.identifier()):
            out_beta_arcs = chain(
                out_beta_arcs,
                iter((self.beta_t((fragment, network.SINK_VERTEX)),)),
            )
        return gp.quicksum(out_beta_arcs)


class TreeEdges:
    """Undirected tree variables.

    Defines:

    * `beta_uv` is the number of vertices in tree defined by root `v`
      for all `(u, v)` in link, source and sink arcs.
      `beta_vu` is also defined for all link except when the predecessor is `s`
      or when the successor is `t`
    """

    def __init__(
        self,
        network: net.Network,
        model: gp.Model,
        beta_domain: Domain,
        prefix_name: str | None = None,
    ) -> None:
        if prefix_name is None:
            prefix_name = ""
        else:
            prefix_name += "_"
        self.__dir_tree_vars = TreeArcs(network, model, beta_domain, prefix_name)
        self.__beta_rev: dict[str, gp.Var] = dict(
            gen_vars(
                model,
                beta_domain,
                f"{prefix_name}beta_rev",
                # FIXME is there a pb because of loop/double arcs?
                (
                    net.StrFormatter.arc(
                        gfa_link.Link(link_arc.successor(), link_arc.predecessor()),
                    )
                    for link_arc in network.link_arcs()
                ),
            ),
        )

    def dtree(self) -> TreeArcs:
        """Get directed tree variables."""
        return self.__dir_tree_vars

    def beta_rev(self, link_arc: gfa_link.Link) -> gp.Var:
        """Get beta_a variable, a in link-arcs."""
        rev_link_arc = gfa_link.Link(link_arc.successor(), link_arc.predecessor())
        return self.__beta_rev[net.StrFormatter.arc(rev_link_arc)]

    def incoming_beta_rev(
        self,
        network: net.Network,
        fragment: gfa_segment.OrientedFragment,
    ) -> gp.LinExpr:
        """Get sum over incoming `beta_vu` where `u != s`."""
        return gp.quicksum(
            (
                self.beta_rev(link)
                for link in gfa_link.incoming_links(network.gfa_graph(), fragment)
            ),
        )

    def outgoing_beta_rev(
        self,
        network: net.Network,
        fragment: gfa_segment.OrientedFragment,
    ) -> gp.LinExpr:
        """Get sum over outgoing `beta_wv` where `w != t`."""
        return gp.quicksum(
            (
                self.beta_rev(link)
                for link in gfa_link.outgoing_links(network.gfa_graph(), fragment)
            ),
        )


class Root:
    """Source vertex variables."""

    def __init__(
        self,
        network: net.Network,
        model: gp.Model,
        r_domain: Domain,
        prefix_name: str | None = None,
    ) -> None:
        if prefix_name is None:
            prefix_name = ""
        else:
            prefix_name += "_"
        self.__r: dict[str, gp.Var] = dict(
            gen_vars(
                model,
                r_domain,
                f"{prefix_name}r",
                map(str, network.oriented_fragments()),
            ),
        )

    def r(self, orient_frag: gfa_segment.OrientedFragment) -> gp.Var:
        """Get r variable for vertex v."""
        return self.__r[str(orient_frag)]
