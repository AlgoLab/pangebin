"""PangeBin-flow input-output module."""

from __future__ import annotations

import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any

import typer

from pangebin.yaml import YAMLInterface


class Manager:
    """PangeBin-flow input/output manager."""

    __BIN_DIR_PREFIX = "bin"
    __BIN_STATS_FILENAME = Path("bin_stats.yaml")
    __MILP_STATS_FILENAME = Path("milp_stats.yaml")
    __BIN_SEQ_NORMCOV_FILENAME = Path("bin_seq_normcov.tsv")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def bin_directory(self, iteration: int) -> Path:
        """Get bin output directory."""
        return self.__config.output_directory() / f"{self.__BIN_DIR_PREFIX}_{iteration}"

    def bin_stats_path(self, iteration: int) -> Path:
        """Get bin stats YAML file path."""
        return self.bin_directory(iteration) / self.__BIN_STATS_FILENAME

    def milp_stats_path(self, iteration: int) -> Path:
        """Get MILP stats YAML file path."""
        return self.bin_directory(iteration) / self.__MILP_STATS_FILENAME

    def bin_seq_normcov_path(self, iteration: int) -> Path:
        """Get bin stats YAML file path."""
        return self.bin_directory(iteration) / self.__BIN_SEQ_NORMCOV_FILENAME

    def gurobi_log_path(self, iteration: int, model_name: str) -> Path:
        """Get Gurobi log file path."""
        return self.bin_directory(iteration) / f"{model_name}.log"

    def move_gurobi_logs(
        self,
        log_files: list[Path],
        fn_log_file_to_iteration_and_model_name: Callable[[Path], tuple[int, str]],
    ) -> None:
        """Move Gurobi log files to the bin directory."""
        for log_file in log_files:
            shutil.move(
                log_file,
                self.gurobi_log_path(
                    *fn_log_file_to_iteration_and_model_name(log_file),
                ),
            )

    def number_of_bins(self) -> int:
        """Get number of bins."""
        return len(
            list(self.__config.output_directory().glob(f"{self.__BIN_DIR_PREFIX}_*")),
        )

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config(YAMLInterface):
    """PangeBin-flow config class."""

    DEFAULT_OUTPUT_DIR = Path("./plasbin")

    KEY_OUTPUT_DIR = "output_directory"

    NAME = "PangeBin-flow IO config"

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


class IOOptions:
    """Input-output options."""

    _RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        help="Output directory",
        rich_help_panel=_RICH_HELP_PANEL,
    )
