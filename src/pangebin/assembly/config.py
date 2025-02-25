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
    KEY_SPADES_MAX_MEMORY = "spades_max_memory"

    DEFAULT_NUMBER_THREADS = 4
    DEFAULT_SPADES_MAX_MEMORY = 250

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
            config_dict.get(cls.KEY_SPADES_MAX_MEMORY, cls.DEFAULT_SPADES_MAX_MEMORY),
        )

    def __init__(
        self,
        number_threads: int = DEFAULT_NUMBER_THREADS,
        spades_max_memory: int = DEFAULT_SPADES_MAX_MEMORY,
    ) -> None:
        """Initialize object.

        Parameters
        ----------
        number_threads : int, optional
            Number of threads
        spades_max_memory : int, optional
            Max memory usage (in GB)
        """
        self.__number_threads = number_threads
        self.__spades_max_memory = spades_max_memory

    def number_threads(self) -> int:
        """Get number of threads option."""
        return self.__number_threads

    def spades_max_memory(self) -> int:
        """Get spades max memory option."""
        return self.__spades_max_memory

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_NUMBER_THREADS: self.__number_threads,
            self.KEY_SPADES_MAX_MEMORY: self.__spades_max_memory,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        with yaml_filepath.open("w") as yaml_file:
            yaml.dump(self.to_dict(), yaml_file, Dumper=Dumper)
        return yaml_filepath
