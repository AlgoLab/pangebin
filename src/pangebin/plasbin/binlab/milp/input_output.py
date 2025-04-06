"""Plasbin binlab MILP input/output module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pangebin.plasbin.binlab.milp.models as lp_mod
from pangebin.yaml import YAMLInterface


class Manager(YAMLInterface):
    """MILP input/output manager."""

    DEFAULT_OUTPUT_DIR = Path("./milp")

    KEY_OUTPUT_DIR = "output_directory"

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Manager:
        """Convert dict to object."""
        return Manager(config_dict.get(cls.KEY_OUTPUT_DIR, Manager.DEFAULT_OUTPUT_DIR))

    def __init__(self, output_directory: Path = DEFAULT_OUTPUT_DIR) -> None:
        """Initialize object."""
        self.__output_directory = output_directory

    def output_directory(self) -> Path:
        """Get output directory."""
        return self.__output_directory

    def gurobi_log_path(
        self,
        iteration: int,
        model: lp_mod.Names,
    ) -> Path:
        """Get Gurobi log file path."""
        return self.__output_directory / f"{iteration}_{model}.log"

    @classmethod
    def attributes_from_gurobi_log_path(
        cls,
        gurobi_log_path: Path,
    ) -> tuple[int, lp_mod.Names]:
        """Get iteration and model from Gurobi log file path."""
        pos_first_sep = gurobi_log_path.stem.find("_")
        iteration = int(gurobi_log_path.stem[:pos_first_sep])
        model_name = gurobi_log_path.stem[pos_first_sep + 1 :]
        return iteration, lp_mod.Names(model_name)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {self.KEY_OUTPUT_DIR: self.__output_directory}
