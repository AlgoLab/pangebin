"""Ground truth config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class Config:
    """Ground truth config class."""

    KEY_MIN_PIDENT = "min_pident"
    KEY_MIN_CONTIG_COVERAGE = "min_contig_coverage"

    DEFAULT_MIN_PIDENT = 95
    DEFAULT_MIN_CONTIG_COVERAGE = 0.95

    DEFAULT_YAML_FILE = Path("ground_truth_config.yaml")

    NAME = "Ground truth config"

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> Config:
        """Create config instance from a YAML file."""
        with Path(yaml_filepath).open("r") as file:
            config_data = yaml.safe_load(file)
        return cls.from_dict(config_data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_MIN_PIDENT, cls.DEFAULT_MIN_PIDENT),
            config_dict.get(
                cls.KEY_MIN_CONTIG_COVERAGE,
                cls.DEFAULT_MIN_CONTIG_COVERAGE,
            ),
        )

    def __init__(
        self,
        min_pident: float = DEFAULT_MIN_PIDENT,
        min_contig_coverage: float = DEFAULT_MIN_CONTIG_COVERAGE,
    ) -> None:
        """Initialize object.

        Parameters
        ----------
        min_pident : float
            Pourcentage of identity threshold (between 0 and 100)
        min_contig_coverage : float
            Contig coverage threshold (between 0 and 1)
        """
        self.__min_pident = min_pident
        self.__min_contig_coverage = min_contig_coverage

    def min_pident(self) -> float:
        """Get min pident."""
        return self.__min_pident

    def min_contig_coverage(self) -> float:
        """Get min contig coverage."""
        return self.__min_contig_coverage

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MIN_PIDENT: self.__min_pident,
            self.KEY_MIN_CONTIG_COVERAGE: self.__min_contig_coverage,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        with yaml_filepath.open("w") as yaml_file:
            yaml.dump(self.to_dict(), yaml_file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
