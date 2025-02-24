"""Assembly config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class Unicycler:
    """Unicycler configuration."""

    KEY_NUMBER_THREADS = "number_threads"

    DEFAULT_NUMBER_THREADS = 4

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> Unicycler:
        """Create config instance from a YAML file."""
        with Path(yaml_filepath).open("r") as file:
            config_data = yaml.safe_load(file)
        return cls.from_dict(config_data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Unicycler:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_NUMBER_THREADS, cls.DEFAULT_NUMBER_THREADS),
        )

    def __init__(self, number_threads: int = DEFAULT_NUMBER_THREADS) -> None:
        self.__number_threads = number_threads

    def number_threads(self) -> int:
        """Get number of threads option."""
        return self.__number_threads

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_NUMBER_THREADS: self.__number_threads,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        with yaml_filepath.open("w") as yaml_file:
            yaml.dump(self.to_dict(), yaml_file, Dumper=Dumper)
        return yaml_filepath
