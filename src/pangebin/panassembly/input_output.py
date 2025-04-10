"""Pan-assembly IO module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pangebin.yaml_interface import YAMLInterface


class Manager:
    """Pan-assembly input/output manager."""

    PANASSEMBLY_GFA_FILENAME = Path("panassembly.gfa")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def panassembly_gfa_path(self) -> Path:
        """Get pan-assembly GFA path."""
        return self.__config.output_directory() / self.PANASSEMBLY_GFA_FILENAME

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config(YAMLInterface):
    """Pan-assembly config class."""

    DEFAULT_DIR = Path("./panassembly")

    KEY_OUTPUT_DIR = "output_directory"

    NAME = "Pan-assembly IO config"

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_OUTPUT_DIR, cls.DEFAULT_DIR),
        )

    def __init__(self, output_directory: Path = DEFAULT_DIR) -> None:
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
