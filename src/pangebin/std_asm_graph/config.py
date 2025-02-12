"""Standardize config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class Config:
    """Standardize config class."""

    DEFAULT_MIN_CONTIG_LENGTH = 1

    KEY_MIN_CONTIG_LENGTH = "min_contig_length"

    DEFAULT_YAML_FILE = Path("standardize_config.yaml")

    NAME = "Standardize config"

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
            config_dict.get(cls.KEY_MIN_CONTIG_LENGTH, cls.DEFAULT_MIN_CONTIG_LENGTH),
        )

    def __init__(
        self,
        min_contig_length: int = DEFAULT_MIN_CONTIG_LENGTH,
    ) -> None:
        """Initialize object."""
        self.__min_contig_length = min_contig_length

    def min_contig_length(self) -> int:
        """Get min contig length."""
        return self.__min_contig_length

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MIN_CONTIG_LENGTH: self.__min_contig_length,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
