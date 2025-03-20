"""Seed contig thresholds input-output module."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pangebin.yaml import YAMLInterface

if TYPE_CHECKING:
    from collections.abc import Iterator


import pangebin.seed.thresholds.items as seed_thr_items

_LOGGER = logging.getLogger(__name__)


# REFACTOR use Reader
def test_items_from_file(
    file: Path,
) -> Iterator[seed_thr_items.TestItem]:
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
            yield seed_thr_items.TestItem.from_datatest_line(line)


class Manager:
    """Seed thresholds IO manager."""

    SEED_THRESHOLD_YAMLNAME = Path("seed_thresholds.yaml")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def seed_threshold_yaml(self) -> Path:
        """Get threshold YAML path."""
        return self.__config.output_directory() / self.SEED_THRESHOLD_YAMLNAME

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config(YAMLInterface):
    """Seed thresholds IO configuration."""

    DEFAULT_OUTPUT_DIR = Path("./seed_thresholds")

    KEY_OUTPUT_DIR = "output_directory"

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
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
