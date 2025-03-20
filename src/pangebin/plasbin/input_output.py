"""PangeBin-Flow input-output module."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml  # type: ignore[import-untyped]

import pangebin.plasbin.milp.input_output as milp_io

if TYPE_CHECKING:
    import pangebin.plasbin.milp.models as milp_models

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class Manager:
    """PangeBin-Flow input/output manager."""

    __BIN_DIR_PREFIX = "bin"
    __BIN_STATS_FILENAME = Path("bin_stats.yaml")
    __BIN_SEQ_NORMCOV_FILENAME = Path("bin_seq_normcov.tsv")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def bin_outdir(self, iteration: int) -> Path:
        """Get bin output directory."""
        return self.__config.output_directory() / f"{self.__BIN_DIR_PREFIX}_{iteration}"

    def bin_stats_path(self, iteration: int) -> Path:
        """Get bin stats YAML file path."""
        return self.bin_outdir(iteration) / self.__BIN_STATS_FILENAME

    def bin_seq_normcov_path(self, iteration: int) -> Path:
        """Get bin stats YAML file path."""
        return self.bin_outdir(iteration) / self.__BIN_SEQ_NORMCOV_FILENAME

    def gurobi_log_path(self, iteration: int, model: milp_models.Names) -> Path:
        """Get Gurobi log file path."""
        return self.bin_outdir(iteration) / f"{model}.log"

    def move_gurobi_logs(self, log_files: list[Path]) -> None:
        """Move Gurobi log files to the bin directory."""
        for log_file in log_files:
            shutil.move(
                log_file,
                self.gurobi_log_path(
                    *milp_io.Manager.attributes_from_gurobi_log_path(log_file),
                ),
            )

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config:
    """PangeBin-Flow config class."""

    DEFAULT_OUTPUT_DIR = Path("./plasbin")

    KEY_OUTPUT_DIR = "output_directory"

    NAME = "PangeBin-Flow IO config"

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
