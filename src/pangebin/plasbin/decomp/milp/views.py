"""MILP views."""

from __future__ import annotations

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.decomp.milp.models as lp_mod
import pangebin.plasbin.decomp.milp.objectives as lp_obj
import pangebin.plasbin.decomp.milp.variables as lp_vars
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net
from pangebin.yaml_interface import YAMLInterface


class MCFStats:
    """MCF stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_COVERAGE_SCORE = "coverage_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MCFStats:
        """Convert dict to object."""
        return cls(obj_dict[cls.KEY_TOTAL_FLOW], obj_dict[cls.KEY_COVERAGE_SCORE])

    def __init__(self, total_flow: float, coverage_score: float) -> None:
        """MCF stats."""
        self.__total_flow = total_flow
        self.__coverage_score = coverage_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def coverage_score(self) -> float:
        """Coverage score."""
        return self.__coverage_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_COVERAGE_SCORE: self.__coverage_score,
        }


def mcf_stats_from_opt_vars(
    network: net.Network,
    mcf_vars: lp_vars.MaxCovFlow,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> MCFStats:
    """Create MCF stats from optimal variables."""
    return MCFStats(
        mcf_vars.flow().total().X,
        lp_obj.coverage_score(
            network,
            mcf_vars.flow(),
            mcf_vars.frag(),
            obj_fun_domain,
        ).getValue(),
    )


class MGCStats(MCFStats):
    """MGC stats."""

    KEY_GC_SCORE = "gc_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MGCStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_COVERAGE_SCORE],
            obj_dict[cls.KEY_GC_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        coverage_score: float,
        gc_score: float,
    ) -> None:
        """MGC stats."""
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


def mgc_stats_from_opt_vars(
    network: net.Network,
    intervals: gc_items.Intervals,
    mgc_vars: lp_vars.MaxGC,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> MGCStats:
    """Create MGC stats from optimal variables."""
    return MGCStats(
        mgc_vars.flow().total().X,
        lp_obj.coverage_score(
            network,
            mgc_vars.flow(),
            mgc_vars.frag(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.gc_score(
            network,
            intervals,
            mgc_vars.frag_gc(),
            obj_fun_domain,
        ).getValue(),
    )


class MPSStats(MGCStats):
    """MPS stats."""

    KEY_PLASMIDNESS_SCORE = "plasmidness_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MPSStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_COVERAGE_SCORE],
            obj_dict[cls.KEY_GC_SCORE],
            obj_dict[cls.KEY_PLASMIDNESS_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        coverage_score: float,
        gc_score: float,
        plasmidness_score: float,
    ) -> None:
        """MPS stats."""
        super().__init__(total_flow, coverage_score, gc_score)
        self.__plasmidness_score = plasmidness_score

    def plasmidness_score(self) -> float:
        """Plasmidness score."""
        return self.__plasmidness_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return super().to_dict() | {
            self.KEY_PLASMIDNESS_SCORE: self.__plasmidness_score,
        }


def mps_stats_from_opt_vars(
    network: net.Network,
    intervals: gc_items.Intervals,
    mps_vars: lp_vars.MaxPlasmidScore,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> MPSStats:
    """Create MPS stats from optimal variables."""
    return MPSStats(
        mps_vars.flow().total().X,
        lp_obj.coverage_score(
            network,
            mps_vars.flow(),
            mps_vars.frag(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.gc_score(
            network,
            intervals,
            mps_vars.frag_gc(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.plasmidness_score(
            network,
            mps_vars.frag(),
            obj_fun_domain,
        ).getValue(),
    )


class MRCFStats(MPSStats):
    """MRCF stats."""

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MRCFStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_COVERAGE_SCORE],
            obj_dict[cls.KEY_GC_SCORE],
            obj_dict[cls.KEY_PLASMIDNESS_SCORE],
        )


def mrcf_stats_from_opt_vars(
    network: net.Network,
    intervals: gc_items.Intervals,
    mrcf_vars: lp_vars.MaxRefCovFlow,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> MRCFStats:
    """Create MRCF stats from optimal variables."""
    return MRCFStats(
        mrcf_vars.flow().total().X,
        lp_obj.coverage_score(
            network,
            mrcf_vars.flow(),
            mrcf_vars.frag(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.gc_score(
            network,
            intervals,
            mrcf_vars.frag_gc(),
            obj_fun_domain,
        ).getValue(),
        lp_obj.plasmidness_score(
            network,
            mrcf_vars.frag(),
            obj_fun_domain,
        ).getValue(),
    )


class StatsContainer(YAMLInterface):
    """Stats container."""

    KEY_MCF_STATS = str(lp_mod.Names.MCF)
    KEY_MGC_STATS = str(lp_mod.Names.MGC)
    KEY_MPS_STATS = str(lp_mod.Names.MPS)
    KEY_MRCF_STATS = str(lp_mod.Names.MRCF)

    @classmethod
    def from_dict(cls, obj_dict: dict) -> StatsContainer:
        """Convert dict to object."""
        return cls(
            MCFStats.from_dict(obj_dict[cls.KEY_MCF_STATS]),
            MGCStats.from_dict(obj_dict[cls.KEY_MGC_STATS]),
            MPSStats.from_dict(obj_dict[cls.KEY_MPS_STATS]),
            MRCFStats.from_dict(obj_dict[cls.KEY_MRCF_STATS]),
        )

    def __init__(
        self,
        mcf_stats: MCFStats,
        mgc_stats: MGCStats,
        mps_stats: MPSStats,
        mrcf_stats: MRCFStats,
    ) -> None:
        """Initialize object."""
        self.__mcf_stats = mcf_stats
        self.__mgc_stats = mgc_stats
        self.__mps_stats = mps_stats
        self.__mrcf_stats = mrcf_stats

    def mcf(self) -> MCFStats:
        """MCF stats."""
        return self.__mcf_stats

    def mgc(self) -> MGCStats:
        """MGC stats."""
        return self.__mgc_stats

    def mps(self) -> MPSStats:
        """MPS stats."""
        return self.__mps_stats

    def mrcf(self) -> MRCFStats:
        """MRCF stats."""
        return self.__mrcf_stats

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_MCF_STATS: self.__mcf_stats.to_dict(),
            self.KEY_MGC_STATS: self.__mgc_stats.to_dict(),
            self.KEY_MPS_STATS: self.__mps_stats.to_dict(),
            self.KEY_MRCF_STATS: self.__mrcf_stats.to_dict(),
        }
