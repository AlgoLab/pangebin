"""Plasbin binlab MILP views."""

from __future__ import annotations

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.binlab.milp.models as lp_mod
import pangebin.plasbin.binlab.milp.objectives as lp_obj
import pangebin.plasbin.binlab.milp.variables as lp_vars
import pangebin.plasbin.milp.objectives as cmn_lp_objs
import pangebin.plasbin.network as net
from pangebin.yaml_interface import YAMLInterface


class MBSStats:
    """MBS stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_BINNING_SCORE = "binning_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MBSStats:
        """Convert dict to object."""
        return cls(obj_dict[cls.KEY_TOTAL_FLOW], obj_dict[cls.KEY_BINNING_SCORE])

    def __init__(self, total_flow: float, binning_score: float) -> None:
        """MBS stats."""
        self.__total_flow = total_flow
        self.__binning_score = binning_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def binning_score(self) -> float:
        """Binning score."""
        return self.__binning_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_BINNING_SCORE: self.__binning_score,
        }


def mbs_stats_from_opt_vars(
    network: net.Network,
    mbs_vars: lp_vars.MaxBinScore,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> MBSStats:
    """Create MBS stats from optimal variables."""
    return MBSStats(
        mbs_vars.flow().total().X,
        lp_obj.binning_score(
            network,
            mbs_vars.flow(),
            mbs_vars.frag(),
            obj_fun_domain,
        ).getValue(),
    )


class MLSStats(MBSStats):
    """MLS stats."""

    KEY_GC_SCORE = "gc_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MLSStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_BINNING_SCORE],
            obj_dict[cls.KEY_GC_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        coverage_score: float,
        gc_score: float,
    ) -> None:
        """MLS stats."""
        super().__init__(total_flow, coverage_score)
        self.__gc_score = gc_score

    def gc_score(self) -> float:
        """GC score."""
        return self.__gc_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return super().to_dict() | {
            self.KEY_GC_SCORE: self.__gc_score,
        }


def mls_stats_from_opt_vars(
    network: net.Network,
    intervals: gc_items.Intervals,
    mls_vars: lp_vars.MaxLabScore,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> MLSStats:
    """Create MLS stats from optimal variables."""
    return MLSStats(
        mls_vars.flow().total().X,
        lp_obj.binning_score(
            network,
            mls_vars.flow(),
            mls_vars.frag(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.gc_score(
            network,
            intervals,
            mls_vars.frag_gc(),
            obj_fun_domain,
        ).getValue(),
    )


class MRBSStats(MLSStats):
    """MRBS stats."""

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MRBSStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_BINNING_SCORE],
            obj_dict[cls.KEY_GC_SCORE],
        )


def mrbs_stats_from_opt_vars(
    network: net.Network,
    intervals: gc_items.Intervals,
    mls_vars: lp_vars.MaxLabScore,
    obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
) -> MRBSStats:
    """Create MLS stats from optimal variables."""
    return MRBSStats(
        mls_vars.flow().total().X,
        lp_obj.binning_score(
            network,
            mls_vars.flow(),
            mls_vars.frag(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.gc_score(
            network,
            intervals,
            mls_vars.frag_gc(),
            obj_fun_domain,
        ).getValue(),
    )


class StatsContainer(YAMLInterface):
    """Stats container."""

    KEY_MBS_STATS = str(lp_mod.Names.MBS)
    KEY_MLS_STATS = str(lp_mod.Names.MLS)
    KEY_MRBS_STATS = str(lp_mod.Names.MRBS)

    @classmethod
    def from_dict(cls, obj_dict: dict) -> StatsContainer:
        """Convert dict to object."""
        return cls(
            MBSStats.from_dict(obj_dict[cls.KEY_MBS_STATS]),
            MLSStats.from_dict(obj_dict[cls.KEY_MLS_STATS]),
            MRBSStats.from_dict(obj_dict[cls.KEY_MRBS_STATS]),
        )

    def __init__(
        self,
        mbs_stats: MBSStats,
        mls_stats: MLSStats,
        mrbs_stats: MRBSStats,
    ) -> None:
        """Initialize object."""
        self.__mbs_stats = mbs_stats
        self.__mls_stats = mls_stats
        self.__mrbs_stats = mrbs_stats

    def mbs(self) -> MBSStats:
        """MBS stats."""
        return self.__mbs_stats

    def mls(self) -> MLSStats:
        """MLS stats."""
        return self.__mls_stats

    def mrbs(self) -> MRBSStats:
        """MRBS stats."""
        return self.__mrbs_stats

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_MBS_STATS: self.__mbs_stats.to_dict(),
            self.KEY_MLS_STATS: self.__mls_stats.to_dict(),
            self.KEY_MRBS_STATS: self.__mrbs_stats.to_dict(),
        }
