"""HMF bins module."""

from __future__ import annotations

from enum import StrEnum

from pangebin.yaml_interface import YAMLInterface


class Topology(StrEnum):
    """Bin type."""

    CIRCULAR = "circular"
    PARTIALLY_CIRCULAR = "partially_circular"


class SeedConstraint(StrEnum):
    """Seed constraint."""

    REQUIRED = "must_have_seed"
    NOT_REQUIRED = "can_be_free_of_seed"


class Stats(YAMLInterface):
    """Bin stats."""

    KEY_BIN_TYPE = "bin_type"
    KEY_SEED_STATUS = "seed_status"
    KEY_CUMULATIVE_SEQUENCE_LENGTH = "cumulative_sequence_length"
    KEY_TOTAL_FLOW = "total_flow"
    KEY_NORMALIZING_COVERAGE = "normalizing_coverage"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> Stats:
        """Convert dict to object."""
        return cls(
            Topology(obj_dict[cls.KEY_BIN_TYPE]),
            SeedConstraint(obj_dict[cls.KEY_SEED_STATUS]),
            obj_dict[cls.KEY_CUMULATIVE_SEQUENCE_LENGTH],
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_NORMALIZING_COVERAGE],
        )

    def __init__(
        self,
        bin_type: Topology,
        seed_status: SeedConstraint,
        cumultative_sequence_length: int,
        total_flow: float,
        normalizing_coverage: float,
    ) -> None:
        """Initialize object."""
        self.__bin_type = bin_type
        self.__seed_status = seed_status
        self.__cumultative_sequence_length = cumultative_sequence_length
        self.__total_flow = total_flow
        self.__normalizing_coverage = normalizing_coverage

    def bin_type(self) -> Topology:
        """Bin type."""
        return self.__bin_type

    def seed_status(self) -> SeedConstraint:
        """Seed status."""
        return self.__seed_status

    def cumultative_sequence_length(self) -> int:
        """Cumulative sequence length."""
        return self.__cumultative_sequence_length

    def total_flow(self) -> float:
        """Total flow."""
        return self.__total_flow

    def normalizing_coverage(self) -> float:
        """Get Normalized flow."""
        return self.__normalizing_coverage

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_BIN_TYPE: self.__bin_type.value,
            self.KEY_CUMULATIVE_SEQUENCE_LENGTH: self.__cumultative_sequence_length,
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_NORMALIZING_COVERAGE: self.__normalizing_coverage,
        }
