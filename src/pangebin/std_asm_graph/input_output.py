"""Standardize IO module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class Manager:
    """Standardize input/output manager."""

    UNICYCLER_GFA_FILENAME = Path("unicycler.gfa")
    SKESA_GFA_FILENAME = Path("skesa.gfa")

    UNICYCLER_FASTA_FILENAME = Path("unicycler.fasta")
    SKESA_FASTA_FILENAME = Path("skesa.fasta")
    MIXED_FASTA_FILENAME = Path("mixed.fasta")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def unicycler_gfa_path(self) -> Path:
        """Get unicycler GFA path."""
        return self.__config.output_directory() / self.UNICYCLER_GFA_FILENAME

    def skesa_gfa_path(self) -> Path:
        """Get skesa GFA path."""
        return self.__config.output_directory() / self.SKESA_GFA_FILENAME

    def unicycler_fasta_path(self) -> Path:
        """Get unicycler FASTA path."""
        return self.__config.output_directory() / self.UNICYCLER_FASTA_FILENAME

    def skesa_fasta_path(self) -> Path:
        """Get skesa FASTA path."""
        return self.__config.output_directory() / self.SKESA_FASTA_FILENAME

    def mixed_fasta_path(self) -> Path:
        """Get mixed FASTA path."""
        return self.__config.output_directory() / self.MIXED_FASTA_FILENAME

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config:
    """Standardize IO config class."""

    DEFAULT_DIR = Path("./standardize")

    KEY_OUTPUT_DIR = "output_directory"

    NAME = "Standardize IO config"

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

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        yaml_filepath = Path(yaml_filepath)
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
