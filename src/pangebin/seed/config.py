"""Seed sequences configuration module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class ThresholdRanges:
    """Threshold ranges."""

    KEY_LENGTH = "length"
    KEY_GENE_DENSITY = "gene_density"

    KEY_MIN = "min"
    KEY_MAX = "max"
    KEY_STEP = "step"

    DEFAULT_MIN_LENGTH = 50
    DEFAULT_MAX_LENGTH = 5000
    DEFAULT_STEP_LENGTH = 50
    DEFAULT_MIN_GENE_DENSITY = 0.01
    DEFAULT_MAX_GENE_DENSITY = 1.0
    DEFAULT_STEP_GENE_DENSITY = 0.01

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> ThresholdRanges:
        """Create config instance from a YAML file."""
        with Path(yaml_filepath).open("r") as file:
            config_data = yaml.safe_load(file)
        return cls.from_dict(config_data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, dict[str, Any]]) -> ThresholdRanges:
        """Convert dict to object."""
        return cls(
            (
                config_dict[cls.KEY_LENGTH].get(cls.KEY_MIN, cls.DEFAULT_MIN_LENGTH)
                if cls.KEY_LENGTH in config_dict
                else cls.DEFAULT_MIN_LENGTH
            ),
            (
                config_dict[cls.KEY_LENGTH].get(cls.KEY_MAX, cls.DEFAULT_MAX_LENGTH)
                if cls.KEY_LENGTH in config_dict
                else cls.DEFAULT_MAX_LENGTH
            ),
            (
                config_dict[cls.KEY_LENGTH].get(cls.KEY_STEP, cls.DEFAULT_STEP_LENGTH)
                if cls.KEY_LENGTH in config_dict
                else cls.DEFAULT_STEP_LENGTH
            ),
            (
                config_dict[cls.KEY_GENE_DENSITY].get(
                    cls.KEY_MIN,
                    cls.DEFAULT_MIN_GENE_DENSITY,
                )
                if cls.KEY_GENE_DENSITY in config_dict
                else cls.DEFAULT_MIN_GENE_DENSITY
            ),
            (
                config_dict[cls.KEY_GENE_DENSITY].get(
                    cls.KEY_MAX,
                    cls.DEFAULT_MAX_GENE_DENSITY,
                )
                if cls.KEY_GENE_DENSITY in config_dict
                else cls.DEFAULT_MAX_GENE_DENSITY
            ),
            (
                config_dict[cls.KEY_GENE_DENSITY].get(
                    cls.KEY_STEP,
                    cls.DEFAULT_STEP_GENE_DENSITY,
                )
                if cls.KEY_GENE_DENSITY in config_dict
                else cls.DEFAULT_STEP_GENE_DENSITY
            ),
        )

    def __init__(  # noqa: PLR0913
        self,
        min_length: int = DEFAULT_MIN_LENGTH,
        max_length: int = DEFAULT_MAX_LENGTH,
        step_length: int = DEFAULT_STEP_LENGTH,
        min_gene_density: float = DEFAULT_MIN_GENE_DENSITY,
        max_gene_density: float = DEFAULT_MAX_GENE_DENSITY,
        step_gene_density: float = DEFAULT_STEP_GENE_DENSITY,
    ) -> None:
        self.__min_length = min_length
        self.__max_length = max_length
        self.__step_length = step_length
        self.__min_gene_density = min_gene_density
        self.__max_gene_density = max_gene_density
        self.__step_gene_density = step_gene_density

    def min_length(self) -> int:
        """Minimum length."""
        return self.__min_length

    def max_length(self) -> int:
        """Maximum length."""
        return self.__max_length

    def step_length(self) -> int:
        """Step length."""
        return self.__step_length

    def min_gene_density(self) -> float:
        """Minimum gene density."""
        return self.__min_gene_density

    def max_gene_density(self) -> float:
        """Maximum gene density."""
        return self.__max_gene_density

    def step_gene_density(self) -> float:
        """Step gene density."""
        return self.__step_gene_density

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_LENGTH: {
                self.KEY_MIN: self.__min_length,
                self.KEY_MAX: self.__max_length,
                self.KEY_STEP: self.__step_length,
            },
            self.KEY_GENE_DENSITY: {
                self.KEY_MIN: self.__min_gene_density,
                self.KEY_MAX: self.__max_gene_density,
                self.KEY_STEP: self.__step_gene_density,
            },
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        with Path(yaml_filepath).open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper)
        return yaml_filepath
