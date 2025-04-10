"""Database input/output module."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml  # type: ignore[import-untyped]

import pangebin.database.items as db_items
from pangebin.yaml_interface import YAMLInterface

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


def accessions_from_file(path: Path) -> Iterator[str]:
    """Read accessions from file."""
    with path.open() as f_in:
        for line in f_in:
            yield line.strip()


def illumina_biosamples_from_file(
    path: Path,
) -> Iterator[db_items.IlluminaBioSamples]:
    """Read illumina biosamples from file."""
    with path.open() as f_in:
        for dict_from_yaml in yaml.safe_load_all(f_in):
            yield db_items.IlluminaBioSamples.from_dict(dict_from_yaml)


def illumina_biosamples_to_file(
    illumina_biosamples: Iterable[db_items.IlluminaBioSamples],
    path: Path,
) -> None:
    """Write illumina biosamples to file."""
    with path.open("w") as f_out:
        yaml.dump_all(
            (sample.to_dict() for sample in illumina_biosamples),
            f_out,
            Dumper=Dumper,
            sort_keys=False,
        )


def non_illumina_biosamples_from_file(
    path: Path,
) -> Iterator[db_items.NonIlluminaBioSamples]:
    """Read non illumina biosamples from file."""
    with path.open() as f_in:
        for dict_from_yaml in yaml.safe_load_all(f_in):
            yield db_items.NonIlluminaBioSamples.from_dict(dict_from_yaml)


def non_illumina_biosamples_to_file(
    non_illumina_biosamples: Iterable[db_items.NonIlluminaBioSamples],
    path: Path,
) -> None:
    """Write non illumina biosamples to file."""
    with path.open("w") as f_out:
        yaml.dump_all(
            (sample.to_dict() for sample in non_illumina_biosamples),
            f_out,
            Dumper=Dumper,
            sort_keys=False,
        )


class Manager:
    """Plasmid database IO manager."""

    ILLUMINA_BIOSAMPLES_YAMLNAME = Path("illumina_biosamples.yaml")
    NON_ILLUMINA_BIOSAMPLES_YAMLNAME = Path("non_illumina_biosamples.yaml")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def illumina_biosamples_yaml(self) -> Path:
        """Get illumina biosamples YAML path."""
        return self.__config.output_directory() / self.ILLUMINA_BIOSAMPLES_YAMLNAME

    def non_illumina_biosamples_yaml(self) -> Path:
        """Get not illumina biosamples YAML path."""
        return self.__config.output_directory() / self.NON_ILLUMINA_BIOSAMPLES_YAMLNAME

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config(YAMLInterface):
    """Plasmid database IO configuration."""

    DEFAULT_OUTPUT_DIR = Path("./database")

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
