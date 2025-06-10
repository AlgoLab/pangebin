"""Plasbin Classbin MILP views."""

from __future__ import annotations

import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net
from pangebin.yaml_interface import YAMLInterface

from . import objectives as mfb_obj
from . import variables as mfb_var


class Stats(YAMLInterface):
    """Classify stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_MIN_NON_ZERO_ARC_FLOW = "min_non_zero_arc_flow"
    KEY_MIN_NON_ZERO_INFLOW = "min_non_zero_inflow"
    KEY_PLASMIDNESS_SCORE = "plasmidness_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> Stats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_MIN_NON_ZERO_ARC_FLOW],
            obj_dict[cls.KEY_MIN_NON_ZERO_INFLOW],
            obj_dict[cls.KEY_PLASMIDNESS_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        min_non_zero_arc_flow: float,
        min_non_zero_inflow: float,
        plasmidness_score: float,
    ) -> None:
        """Classify stats."""
        self.__total_flow = total_flow
        self.__min_non_zero_arc_flow = min_non_zero_arc_flow
        self.__min_non_zero_inflow = min_non_zero_inflow
        self.__plasmidness_score = plasmidness_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def min_non_zero_arc_flow(self) -> float:
        """Min non-zero arc flow."""
        return self.__min_non_zero_arc_flow

    def min_non_zero_inflow(self) -> float:
        """Min non-zero inflow."""
        return self.__min_non_zero_inflow

    def plasmidness_score(self) -> float:
        """Plasmidness score."""
        return self.__plasmidness_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_MIN_NON_ZERO_ARC_FLOW: self.__min_non_zero_arc_flow,
            self.KEY_MIN_NON_ZERO_INFLOW: self.__min_non_zero_inflow,
            self.KEY_PLASMIDNESS_SCORE: self.__plasmidness_score,
        }


def stats_from_opt_bin_vars(
    network: net.Network,
    bin_vars: mfb_var.BinVariables,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
    partially_circular: bool,  # noqa: FBT001
) -> Stats:
    """Create MGCLB stats from optimal variables."""
    # -------------------------------------------------------------------------------- #
    # # FIXME the error is constraits are missing: force circular to have link arcs !
    # # TODO CONTINUE HERE
    # # BUG TMP
    # for link_arc in network.source_arcs():
    #     print(
    #         str((link_arc[0], str(link_arc[1]))),
    #         bin_vars.sub_arcs().s(link_arc).X,
    #         bin_vars.flow().s(link_arc).X,
    #     )
    #     print("--")
    # print(bin_vars.flow().total().X)
    # print("--")
    # -------------------------------------------------------------------------------- #

    return Stats(
        bin_vars.flows().total().X,
        (
            min(
                (
                    bin_vars.flows().l(link_arc).X
                    for link_arc in network.link_arcs()
                    if bin_vars.flows().l(link_arc).X > 0
                ),
                default=0,
            )
            if partially_circular
            else min(
                bin_vars.flows().l(link_arc).X
                for link_arc in network.link_arcs()
                if bin_vars.flows().l(link_arc).X > 0
            )
        ),
        min(
            bin_vars.flows().incoming_forward_reverse(network, frag_id).getValue()
            for frag_id in network.fragment_ids()
            if bin_vars.flows().incoming_forward_reverse(network, frag_id).getValue()
            > 0
        ),
        mfb_obj.plasmidness_score(
            network,
            bin_vars.flows(),
            obj_fun_domain,
        ).getValue(),
    )
