"""Standardize config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pangebin.yaml import YAMLInterface


class Config(YAMLInterface):
    """Standardize config class."""

    DEFAULT_MIN_CONTIG_LENGTH = 1

    KEY_MIN_CONTIG_LENGTH = "min_contig_length"

    DEFAULT_YAML_FILE = Path("standardize_config.yaml")

    NAME = "Standardize config"

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
