"""Bins item."""

from __future__ import annotations

import pangebin.plasbin.milp.models as milp_models
import pangebin.plasbin.milp.views as milp_views
from pangebin.yaml import YAMLInterface


class Stats(YAMLInterface):
    """Bin stats."""

    KEY_CUMULATIVE_SEQUENCE_LENGTH = "cumulative_sequence_length"
    KEY_COVERAGE_FLOW = "coverage_flow"
    KEY_GC_CONTENT_INTERVAL = "GC_content_interval"
    KEY_MILP_STATS = "milp_stats"
    SUBKEY_MCF_STATS = str(milp_models.Names.MCF)
    SUBKEY_MGC_STATS = str(milp_models.Names.MGC)
    SUBKEY_MPS_STATS = str(milp_models.Names.MPS)

    @classmethod
    def from_dict(cls, obj_dict: dict) -> Stats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_CUMULATIVE_SEQUENCE_LENGTH],
            obj_dict[cls.KEY_COVERAGE_FLOW],
            tuple(obj_dict[cls.KEY_GC_CONTENT_INTERVAL]),
            milp_views.MCFStats.from_dict(
                obj_dict[cls.KEY_MILP_STATS][cls.SUBKEY_MCF_STATS],
            ),
            milp_views.MGCStats.from_dict(
                obj_dict[cls.KEY_MILP_STATS][cls.SUBKEY_MGC_STATS],
            ),
            milp_views.MPSStats.from_dict(
                obj_dict[cls.KEY_MILP_STATS][cls.SUBKEY_MPS_STATS],
            ),
        )

    def __init__(  # noqa: PLR0913
        self,
        cumultative_sequence_length: int,
        coverage_flow: float,
        gc_content_interval: tuple[float, float],
        mcf_stats: milp_views.MCFStats,
        mgc_stats: milp_views.MGCStats,
        mps_stats: milp_views.MPSStats,
    ) -> None:
        """Initialize object."""
        self.__cumultative_sequence_length = cumultative_sequence_length
        self.__coverage_flow = coverage_flow
        self.__gc_content_interval = gc_content_interval
        self.__mcf_stats = mcf_stats
        self.__mgc_stats = mgc_stats
        self.__mps_stats = mps_stats

    def cumultative_sequence_length(self) -> int:
        """Cumulative sequence length."""
        return self.__cumultative_sequence_length

    def coverage_flow(self) -> float:
        """Coverage flow."""
        return self.__coverage_flow

    def gc_content_interval(self) -> tuple[float, float]:
        """GC content interval."""
        return self.__gc_content_interval

    def mcf_stats(self) -> milp_views.MCFStats:
        """MCF stats."""
        return self.__mcf_stats

    def mgc_stats(self) -> milp_views.MGCStats:
        """MGC stats."""
        return self.__mgc_stats

    def mps_stats(self) -> milp_views.MPSStats:
        """MPS stats."""
        return self.__mps_stats

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_CUMULATIVE_SEQUENCE_LENGTH: self.__cumultative_sequence_length,
            self.KEY_COVERAGE_FLOW: self.__coverage_flow,
            self.KEY_GC_CONTENT_INTERVAL: list(self.__gc_content_interval),
            self.KEY_MILP_STATS: {
                self.SUBKEY_MCF_STATS: self.__mcf_stats.to_dict(),
                self.SUBKEY_MGC_STATS: self.__mgc_stats.to_dict(),
                self.SUBKEY_MPS_STATS: self.__mps_stats.to_dict(),
            },
        }


class FragmentNormCoverage:
    """Fragment normalized coverage."""

    def __init__(self, identifier: str, normalized_coverage: float) -> None:
        """Fragment normalized coverage."""
        self.__identifier = identifier
        self.__normalized_coverage = normalized_coverage

    def identifier(self) -> str:
        """Fragment identifier."""
        return self.__identifier

    def normalized_coverage(self) -> float:
        """Fragment normalized coverage."""
        return self.__normalized_coverage
