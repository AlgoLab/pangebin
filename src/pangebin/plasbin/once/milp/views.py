"""Plasbin once MILP views."""

from __future__ import annotations

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net
import pangebin.plasbin.once.milp.objectives as lp_obj
import pangebin.plasbin.once.milp.variables as lp_vars
from pangebin.yaml_interface import YAMLInterface


class MGCLBStats(YAMLInterface):
    """MGCLB stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_COVERAGE_SCORE = "coverage_score"
    KEY_PLASMIDNESS_SCORE = "plasmidness_score"
    KEY_GC_SCORE = "gc_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MGCLBStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_COVERAGE_SCORE],
            obj_dict[cls.KEY_PLASMIDNESS_SCORE],
            obj_dict[cls.KEY_GC_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        coverage_score: float,
        plasmidness_score: float,
        gc_score: float,
    ) -> None:
        """MGCLB stats."""
        self.__total_flow = total_flow
        self.__coverage_score = coverage_score
        self.__plasmidness_score = plasmidness_score
        self.__gc_score = gc_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def coverage_score(self) -> float:
        """Binning score."""
        return self.__coverage_score

    def plasmidness_score(self) -> float:
        """Plasmidness score."""
        return self.__plasmidness_score

    def gc_score(self) -> float:
        """GC score."""
        return self.__gc_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_COVERAGE_SCORE: self.__coverage_score,
            self.KEY_PLASMIDNESS_SCORE: self.__plasmidness_score,
            self.KEY_GC_SCORE: self.__gc_score,
        }


def mgclb_stats_from_opt_vars(
    network: net.Network,
    intervals: gc_items.Intervals,
    mgclb_vars: lp_vars.MaxGCLabelBinScore,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> MGCLBStats:
    """Create MGCLB stats from optimal variables."""
    return MGCLBStats(
        mgclb_vars.flow().total().X,
        lp_obj.coverage_score(
            network,
            mgclb_vars.flow(),
            mgclb_vars.frag(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.plasmidness_score(
            network,
            mgclb_vars.flow(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.gc_score(
            network,
            intervals,
            mgclb_vars.inflow_gc(),
            obj_fun_domain,
        ).getValue(),
    )
