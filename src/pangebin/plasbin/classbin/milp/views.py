"""Plasbin Classbin MILP views."""

from __future__ import annotations

import pangebin.plasbin.classbin.milp.objectives as lp_obj
import pangebin.plasbin.classbin.milp.variables as lp_vars
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.network as net
from pangebin.yaml_interface import YAMLInterface


class ClassifyStats(YAMLInterface):
    """Classify stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_PLASMIDNESS_SCORE = "plasmidness_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> ClassifyStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_PLASMIDNESS_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        plasmidness_score: float,
    ) -> None:
        """Classify stats."""
        self.__total_flow = total_flow
        self.__plasmidness_score = plasmidness_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def plasmidness_score(self) -> float:
        """Plasmidness score."""
        return self.__plasmidness_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_PLASMIDNESS_SCORE: self.__plasmidness_score,
        }


def classify_stats_from_opt_vars(
    network: net.Network,
    mgclb_vars: lp_vars.Classify,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> ClassifyStats:
    """Create MGCLB stats from optimal variables."""
    return ClassifyStats(
        mgclb_vars.flow().total().X,
        lp_obj.plasmidness_score(
            network,
            mgclb_vars.flow(),
            obj_fun_domain,
        ).getValue(),
    )
