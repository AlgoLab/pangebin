"""Seed sequence input-output module."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper
import logging

if TYPE_CHECKING:
    from collections.abc import Iterator


import pangebin.seed.items as seed_items

_LOGGER = logging.getLogger(__name__)


def test_items_from_file(
    file: Path,
) -> Iterator[seed_items.SeedContigThresholdTestItem]:
    """Get seed contig threshold test items from file.

    Parameters
    ----------
    file : Path
        Path to the file

    Yield
    -----
    SeedContigThresholdTestItem
        Seed contig threshold test item
    """
    with file.open() as f_in:
        for line in f_in:
            yield seed_items.SeedContigThresholdTestItem.from_datatest_line(line)


class ThresholdsManager:
    """Seed thresholds IO manager."""

    SEED_THRESHOLD_YAMLNAME = Path("seed_thresholds.yaml")

    def __init__(self, config: ThresholdsConfig) -> None:
        """Initialize object."""
        self.__config = config

    def seed_threshold_yaml(self) -> Path:
        """Get threshold YAML path."""
        return self.__config.output_directory() / self.SEED_THRESHOLD_YAMLNAME

    def config(self) -> ThresholdsConfig:
        """Get config."""
        return self.__config


class ThresholdsConfig:
    """Seed thresholds IO configuration."""

    DEFAULT_OUTPUT_DIR = Path("./seed_thresholds")

    KEY_OUTPUT_DIR = "output_directory"

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> ThresholdsConfig:
        """Create config instance from a YAML file."""
        with Path(yaml_filepath).open("r") as file:
            config_data = yaml.safe_load(file)
        return cls.from_dict(config_data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> ThresholdsConfig:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_OUTPUT_DIR, cls.DEFAULT_OUTPUT_DIR),
        )

    def __init__(self, output_directory: Path = DEFAULT_OUTPUT_DIR) -> None:
        """Initialize object."""
        self.__output_directory = output_directory

    def output_directory(self) -> Path:
        """Get output directory."""
        return self.__output_directory

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_OUTPUT_DIR: self.__output_directory,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        yaml_filepath = Path(yaml_filepath)
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
