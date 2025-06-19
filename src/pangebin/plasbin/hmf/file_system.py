"""HMF file system logics."""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path

from . import bins


class Root:
    """HMF approach root file."""

    CCOMP_DIR_PREFIX = "ccomp"
    BEST_INSTANCES_YAML = "best_instances.yaml"

    def __init__(self, root_directory: Path) -> None:
        """Initialize object."""
        self.__root_directory = root_directory

    def dir(self) -> Path:
        """Get root directory path."""
        return self.__root_directory

    def best_instances_yaml(self) -> Path:
        """Get best instances YAML file path."""
        return self.dir() / self.BEST_INSTANCES_YAML

    def ccomp_file_system(self, ccomp_number: int) -> ConnectedComponent:
        """Get connected component file system."""
        return ConnectedComponent(
            self.dir() / f"{self.CCOMP_DIR_PREFIX}_{ccomp_number}",
        )


class ConnectedComponent:
    """Connected component file system."""

    class __TopologyPrefix(StrEnum):
        """Topology prefix."""

        CIRCULAR = "c"
        PARTIALLY_CIRCULAR = "pc"

        @classmethod
        def from_topology(
            cls,
            topology: bins.Topology,
        ) -> ConnectedComponent.__TopologyPrefix:
            """Convert topology to prefix."""
            match topology:
                case bins.Topology.CIRCULAR:
                    return cls.CIRCULAR
                case bins.Topology.PARTIALLY_CIRCULAR:
                    return cls.PARTIALLY_CIRCULAR

    class __SeedConstraintPrefix(StrEnum):
        """Seed constraint prefix."""

        SEED = "s"
        FREE_OF_SEED = "fos"

        @classmethod
        def from_seed_constraint(
            cls,
            seed_constraint: bins.SeedConstraint,
        ) -> ConnectedComponent.__SeedConstraintPrefix:
            """Convert seed constraint to prefix."""
            match seed_constraint:
                case bins.SeedConstraint.REQUIRED:
                    return cls.SEED
                case bins.SeedConstraint.NOT_REQUIRED:
                    return cls.FREE_OF_SEED

    def __init__(self, directory: Path) -> None:
        """Initialize object."""
        self.__directory = directory

    def dir(self) -> Path:
        """Get root directory path."""
        return self.__directory

    def bin_class_file_system(
        self,
        topology: bins.Topology,
        seed_constraint: bins.SeedConstraint,
        number_of_bins: int,
    ) -> BinClass:
        """Get bin class file system."""
        return BinClass(
            self.dir()
            / Path(
                f"{self.__TopologyPrefix.from_topology(topology)}"
                f"_{self.__SeedConstraintPrefix.from_seed_constraint(seed_constraint)}"
                f"_{number_of_bins}",
            ),
        )


class BinClass:
    """Multi-flow bin class file system manager."""

    GUROBI_LOG_NAME = Path("gurobi.log")
    __BIN_DIR_PREFIX = "bin"

    def __init__(self, directory: Path) -> None:
        """Initialize object."""
        self.__directory = directory

    def dir(self) -> Path:
        """Get root directory path."""
        return self.__directory

    def gurobi_log(self) -> Path:
        """Get Gurobi log file path."""
        return self.dir() / self.GUROBI_LOG_NAME

    def bin_file_system(self, bin_number: int) -> Bin:
        """Get bin file system manager."""
        return Bin(self.dir() / f"{self.__BIN_DIR_PREFIX}_{bin_number}")


class Bin:
    """File system for one bin."""

    BIN_STATS_FILENAME = Path("bin_stats.yaml")
    MILP_STATS_FILENAME = Path("milp_stats.yaml")
    BIN_SEQ_NORMCOV_FILENAME = Path("bin_seq_normcov.tsv")

    def __init__(self, directory: Path) -> None:
        """Initialize object."""
        self.__directory = directory

    def dir(self) -> Path:
        """Get directory path."""
        return self.__directory

    def bin_stats_path(self) -> Path:
        """Get bin stats YAML file path."""
        return self.dir() / self.BIN_STATS_FILENAME

    def milp_stats_path(self) -> Path:
        """Get MILP stats YAML file path."""
        return self.dir() / self.MILP_STATS_FILENAME

    def bin_seq_normcov_path(self) -> Path:
        """Get bin sequences with normalized coverages file path."""
        return self.dir() / self.BIN_SEQ_NORMCOV_FILENAME
