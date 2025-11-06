"""HMF configuration module."""

from __future__ import annotations

from typing import Any

import pangebin.plasbin.network as net
from pangebin.yaml_interface import YAMLInterface

from . import bins
from .milp import config as milp_cfg


class _SeedConstraint(YAMLInterface):
    """Seed constraint configuration."""

    KEY_WITH_SEEDS = "with_seeds"
    KEY_MAX_NUMBER_FREE_OF_SEEDS = "max_number_free_of_seeds"

    DEFAULT_WITH_SEEDS = True
    DEFAULT_MAX_NUMBER_FREE_OF_SEEDS = 0

    @classmethod
    def default(cls) -> _SeedConstraint:
        return cls(cls.DEFAULT_WITH_SEEDS, cls.DEFAULT_MAX_NUMBER_FREE_OF_SEEDS)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> _SeedConstraint:
        return cls(
            config_dict[cls.KEY_WITH_SEEDS],
            config_dict[cls.KEY_MAX_NUMBER_FREE_OF_SEEDS],
        )

    def __init__(self, with_seeds: bool, max_number_free_of_seeds: int) -> None:  # noqa: FBT001
        self.__with_seeds = with_seeds
        self.__max_number_free_of_seeds = max_number_free_of_seeds

    def with_seeds(self) -> bool:
        return self.__with_seeds

    def max_number_free_of_seeds(self) -> int:
        return self.__max_number_free_of_seeds

    def to_dict(self) -> dict[str, Any]:
        return {
            self.KEY_WITH_SEEDS: self.__with_seeds,
            self.KEY_MAX_NUMBER_FREE_OF_SEEDS: self.__max_number_free_of_seeds,
        }


class BinProperties(YAMLInterface):
    """HMF bin properties class."""

    KEY_CIRCULAR = "circular"
    KEY_PARTIALLY_CIRCULAR = "partially_circular"

    @classmethod
    def default(cls) -> BinProperties:
        """Get default config."""
        return cls(_SeedConstraint.default(), _SeedConstraint.default())

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> BinProperties:
        """Convert dict to object."""
        return cls(
            _SeedConstraint.from_dict(config_dict[cls.KEY_CIRCULAR])
            if cls.KEY_CIRCULAR in config_dict
            else _SeedConstraint.default(),
            _SeedConstraint.from_dict(config_dict[cls.KEY_PARTIALLY_CIRCULAR])
            if cls.KEY_PARTIALLY_CIRCULAR in config_dict
            else _SeedConstraint.default(),
        )

    def __init__(
        self,
        circular_seed_constraint: _SeedConstraint,
        partially_circular_seed_constraint: _SeedConstraint,
    ) -> None:
        """Initialize object."""
        self.__circular_seed_constraint = circular_seed_constraint
        self.__partially_circular_seed_constraint = partially_circular_seed_constraint

    def circular_seed_constraint(self) -> _SeedConstraint:
        """Get circular seed constraint."""
        return self.__circular_seed_constraint

    def partially_circular_seed_constraint(self) -> _SeedConstraint:
        """Get partially circular seed constraint."""
        return self.__partially_circular_seed_constraint

    def seed_constraint_for_topology(self, topology: bins.Topology) -> _SeedConstraint:
        """Get seed constraint for topology."""
        match topology:
            case bins.Topology.CIRCULAR:
                return self.circular_seed_constraint()
            case bins.Topology.PARTIALLY_CIRCULAR:
                return self.partially_circular_seed_constraint()

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_CIRCULAR: self.__circular_seed_constraint.to_dict(),
            self.KEY_PARTIALLY_CIRCULAR: (
                self.__partially_circular_seed_constraint.to_dict()
            ),
        }


class Network(YAMLInterface):
    """PangeBin-flow network config class."""

    DEFAULT_SINK_ARCS_DOMAIN = net.SinkArcsDomain.ALL

    KEY_SINK_ARC_DEFINITION = "sink_arc_definition"

    @classmethod
    def default(cls) -> Network:
        """Get default config."""
        return cls(cls.DEFAULT_SINK_ARCS_DOMAIN)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Network:
        """Convert dict to object."""
        return cls(
            net.SinkArcsDomain(
                config_dict.get(
                    cls.KEY_SINK_ARC_DEFINITION,
                    cls.DEFAULT_SINK_ARCS_DOMAIN,
                ),
            ),
        )

    def __init__(self, sink_arcs_domain: net.SinkArcsDomain) -> None:
        """Initialize object."""
        self.__sink_arcs_domain = sink_arcs_domain

    def sink_arcs_domain(self) -> net.SinkArcsDomain:
        """Get sink-arcs domain."""
        return self.__sink_arcs_domain

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_SINK_ARC_DEFINITION: str(self.__sink_arcs_domain),
        }


class Config(YAMLInterface):
    """PangeBin-flow binning config class."""

    KEY_BIN_PROPERTIES = "bin_properties"
    KEY_NETWORK = "network"
    KEY_MILP = "milp"

    @classmethod
    def default(cls) -> Config:
        """Get default config."""
        return cls(
            BinProperties.default(),
            Network.default(),
            milp_cfg.Config.default(),
        )

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            BinProperties.from_dict(config_dict[cls.KEY_BIN_PROPERTIES]),
            Network.from_dict(config_dict[cls.KEY_NETWORK]),
            milp_cfg.Config.from_dict(config_dict[cls.KEY_MILP]),
        )

    def __init__(
        self,
        bin_properties: BinProperties,
        network_config: Network,
        milp_config: milp_cfg.Config,
    ) -> None:
        """Initialize object."""
        self.__bin_properties = bin_properties
        self.__network_config = network_config
        self.__milp_config = milp_config

    def bin_properties(self) -> BinProperties:
        """Get bin properties."""
        return self.__bin_properties

    def network(self) -> Network:
        """Get network config."""
        return self.__network_config

    def model(self) -> milp_cfg.Config:
        """Get MILP config."""
        return self.__milp_config

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_BIN_PROPERTIES: self.__bin_properties.to_dict(),
            self.KEY_NETWORK: self.__network_config.to_dict(),
            self.KEY_MILP: self.__milp_config.to_dict(),
        }


if __name__ == "__main__":
    from pathlib import Path

    config_yaml = Path("./hmf_config.yaml")
    Config.default().to_yaml(config_yaml)
