"""MILP views."""

from __future__ import annotations


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


class MGCStats:
    """MGC stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_COVERAGE_SCORE = MCFStats.KEY_COVERAGE_SCORE
    KEY_GC_PROBABILITY_SCORE = "gc_probability_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MGCStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_COVERAGE_SCORE],
            obj_dict[cls.KEY_GC_PROBABILITY_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        coverage_score: float,
        gc_probability_score: float,
    ) -> None:
        """MGC stats."""
        self.__total_flow = total_flow
        self.__coverage_score = coverage_score
        self.__gc_probability_score = gc_probability_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def coverage_score(self) -> float:
        """Coverage score."""
        return self.__coverage_score

    def gc_probability_score(self) -> float:
        """GC probability score."""
        return self.__gc_probability_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_COVERAGE_SCORE: self.__coverage_score,
            self.KEY_GC_PROBABILITY_SCORE: self.__gc_probability_score,
        }


class MPSStats:
    """MPS stats."""

    KEY_TOTAL_FLOW = "total_flow"
    KEY_COVERAGE_SCORE = MCFStats.KEY_COVERAGE_SCORE
    KEY_GC_PROBABILITY_SCORE = MGCStats.KEY_GC_PROBABILITY_SCORE
    KEY_PLASMIDNESS_SCORE = "plasmidness_score"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> MPSStats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_COVERAGE_SCORE],
            obj_dict[cls.KEY_GC_PROBABILITY_SCORE],
            obj_dict[cls.KEY_PLASMIDNESS_SCORE],
        )

    def __init__(
        self,
        total_flow: float,
        coverage_score: float,
        gc_probability_score: float,
        plasmidness_score: float,
    ) -> None:
        """MPS stats."""
        self.__total_flow = total_flow
        self.__coverage_score = coverage_score
        self.__gc_probability_score = gc_probability_score
        self.__plasmidness_score = plasmidness_score

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def coverage_score(self) -> float:
        """Coverage score."""
        return self.__coverage_score

    def gc_probability_score(self) -> float:
        """GC probability score."""
        return self.__gc_probability_score

    def plasmidness_score(self) -> float:
        """Plasmidness score."""
        return self.__plasmidness_score

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_COVERAGE_SCORE: self.__coverage_score,
            self.KEY_GC_PROBABILITY_SCORE: self.__gc_probability_score,
            self.KEY_PLASMIDNESS_SCORE: self.__plasmidness_score,
        }
