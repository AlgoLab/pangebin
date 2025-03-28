"""MILP config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import gurobipy

from pangebin.yaml import YAMLInterface


# TODO add to CLI gurobi config
class Gurobi(YAMLInterface):
    """Gurobi config."""

    DEFAULT_MIP_GAP = None
    DEFAULT_TIME_LIMIT = None
    DEFAULT_THREADS = None

    KEY_MIP_GAP = "mip_gap"
    KEY_TIME_LIMIT = "time_limit"
    KEY_THREADS = "threads"

    DEFAULT_YAML_FILE = Path("gurobi_config.yaml")

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Gurobi:
        """Convert dict to object."""
        return cls(
            mip_gap=config_dict.get(cls.KEY_MIP_GAP, cls.DEFAULT_MIP_GAP),
            time_limit=config_dict.get(cls.KEY_TIME_LIMIT, cls.DEFAULT_TIME_LIMIT),
            threads=config_dict.get(cls.KEY_THREADS, cls.DEFAULT_THREADS),
        )

    def __init__(
        self,
        mip_gap: float | None = DEFAULT_MIP_GAP,
        time_limit: float | None = DEFAULT_TIME_LIMIT,
        threads: int | None = DEFAULT_THREADS,
    ) -> None:
        """Initialize object."""
        self.__mip_gap = mip_gap
        self.__time_limit = time_limit
        self.__threads = threads

    def mip_gap(self) -> float | None:
        """Get mip gap."""
        return self.__mip_gap

    def time_limit(self) -> float | None:
        """Get time limit."""
        return self.__time_limit

    def threads(self) -> int | None:
        """Get threads."""
        return self.__threads

    def to_dict(self) -> dict[str, Any]:
        """Convert object to dict."""
        return {
            self.KEY_MIP_GAP: self.__mip_gap,
            self.KEY_TIME_LIMIT: self.__time_limit,
            self.KEY_THREADS: self.__threads,
        }


def configurate_global_gurobi(config: Gurobi) -> None:
    """Configure global Gurobi parameters."""
    mip_gap = config.mip_gap()
    if mip_gap is not None:
        gurobipy.setParam(gurobipy.GRB.Param.MIPGap, mip_gap)
    time_limit = config.time_limit()
    if time_limit is not None:
        gurobipy.setParam(gurobipy.GRB.Param.TimeLimit, time_limit)
    threads = config.threads()
    if threads is not None:
        gurobipy.setParam(gurobipy.GRB.Param.Threads, threads)


if __name__ == "__main__":
    gurobi_config = Gurobi()
    gurobi_config.to_yaml(Gurobi.DEFAULT_YAML_FILE)
