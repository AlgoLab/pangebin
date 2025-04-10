"""Bins item."""

from __future__ import annotations

from pangebin.yaml_interface import YAMLInterface


class Stats(YAMLInterface):
    """Bin stats."""

    KEY_CUMULATIVE_SEQUENCE_LENGTH = "cumulative_sequence_length"
    KEY_TOTAL_FLOW = "total_flow"
    KEY_NORMALIZING_COVERAGE = "normalizing_coverage"
    KEY_GC_CONTENT_INTERVAL = "GC_content_interval"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> Stats:
        """Convert dict to object."""
        return cls(
            obj_dict[cls.KEY_CUMULATIVE_SEQUENCE_LENGTH],
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_NORMALIZING_COVERAGE],
            tuple(obj_dict[cls.KEY_GC_CONTENT_INTERVAL]),
        )

    def __init__(
        self,
        cumultative_sequence_length: int,
        total_flow: float,
        normalizing_coverage: float,
        gc_content_interval: tuple[float, float],
    ) -> None:
        """Initialize object."""
        self.__cumultative_sequence_length = cumultative_sequence_length
        self.__total_flow = total_flow
        self.__normalizing_coverage = normalizing_coverage
        self.__gc_content_interval = gc_content_interval

    def cumultative_sequence_length(self) -> int:
        """Cumulative sequence length."""
        return self.__cumultative_sequence_length

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def normalizing_coverage(self) -> float:
        """Get Normalized flow."""
        return self.__normalizing_coverage

    def gc_content_interval(self) -> tuple[float, float]:
        """GC content interval."""
        return self.__gc_content_interval

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_CUMULATIVE_SEQUENCE_LENGTH: self.__cumultative_sequence_length,
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_NORMALIZING_COVERAGE: self.__normalizing_coverage,
            self.KEY_GC_CONTENT_INTERVAL: list(self.__gc_content_interval),
        }


class SequenceNormCoverage:
    """Sequence normalized coverage."""

    def __init__(self, identifier: str, normalized_coverage: float) -> None:
        """Sequence normalized coverage."""
        self.__identifier = identifier
        self.__normalized_coverage = normalized_coverage

    def identifier(self) -> str:
        """Sequence identifier."""
        return self.__identifier

    def normalized_coverage(self) -> float:
        """Sequence normalized coverage."""
        return self.__normalized_coverage
