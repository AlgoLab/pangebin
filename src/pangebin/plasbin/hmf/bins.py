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

    KEY_TOPOLOGY = "topology"
    KEY_SEED_CONSTRAINT = "seed_constraint"
    KEY_CUMULATIVE_SEQUENCE_LENGTH = "cumulative_sequence_length"
    KEY_TOTAL_FLOW = "total_flow"
    KEY_NORMALIZING_COVERAGE = "normalizing_coverage"

    @classmethod
    def from_dict(cls, obj_dict: dict) -> Stats:
        """Convert dict to object."""
        return cls(
            Topology(obj_dict[cls.KEY_TOPOLOGY]),
            SeedConstraint(obj_dict[cls.KEY_SEED_CONSTRAINT]),
            obj_dict[cls.KEY_CUMULATIVE_SEQUENCE_LENGTH],
            obj_dict[cls.KEY_TOTAL_FLOW],
            obj_dict[cls.KEY_NORMALIZING_COVERAGE],
        )

    def __init__(
        self,
        topology: Topology,
        seed_constraint: SeedConstraint,
        cumultative_sequence_length: int,
        total_flow: float,
        normalizing_coverage: float,
    ) -> None:
        """Initialize object."""
        self.__topology = topology
        self.__seed_constraint = seed_constraint
        self.__cumultative_sequence_length = cumultative_sequence_length
        self.__total_flow = total_flow
        self.__normalizing_coverage = normalizing_coverage

    def topology(self) -> Topology:
        """Bin topology."""
        return self.__topology

    def seed_constraint(self) -> SeedConstraint:
        """Seed constraint."""
        return self.__seed_constraint

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
            self.KEY_TOPOLOGY: self.__topology.value,
            self.KEY_SEED_CONSTRAINT: self.__seed_constraint.value,
            self.KEY_CUMULATIVE_SEQUENCE_LENGTH: self.__cumultative_sequence_length,
            self.KEY_TOTAL_FLOW: self.__total_flow,
            self.KEY_NORMALIZING_COVERAGE: self.__normalizing_coverage,
        }
