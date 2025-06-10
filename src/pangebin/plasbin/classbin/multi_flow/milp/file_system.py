"""Multi-flow binning MILP file system logics."""

from __future__ import annotations

from pathlib import Path


class Manager:
    """Multi-flow binning file system manager."""

    __NB_FLOW_DIR_PREFIX = "nb_flow"
    __NB_CIRCULAR_DIR_PREFIX = "c"
    __NB_PARTIALLY_CIRCULAR_DIR_PREFIX = "pc"

    def __init__(self, directory: Path) -> None:
        """Initialize object."""
        self.__directory = directory

    def dir(self) -> Path:
        """Get directory path."""
        return self.__directory

    def get_nb_flow_fs_manager(
        self,
        number_of_circular: int,
        number_of_partially_circular: int,
    ) -> NumberOfFlowManager:
        """Get bin file system manager."""
        number_of_flow = number_of_circular + number_of_partially_circular
        return NumberOfFlowManager(
            self.dir()
            / (
                f"{self.__NB_FLOW_DIR_PREFIX}_{number_of_flow}"
                f"_{self.__NB_CIRCULAR_DIR_PREFIX}_{number_of_circular}"
                f"_{self.__NB_PARTIALLY_CIRCULAR_DIR_PREFIX}_{number_of_partially_circular}"
            ),
        )


class NumberOfFlowManager:
    """Number of flow file system manager."""

    GUROBI_LOG_FILENAME = "gurobi.log"
    __BIN_DIR_PREFIX = "bin"

    def __init__(self, directory: Path) -> None:
        """Initialize object."""
        self.__directory = directory

    def dir(self) -> Path:
        """Get directory path."""
        return self.__directory

    def gurobi_log(self) -> Path:
        """Get Gurobi log file path."""
        return self.dir() / self.GUROBI_LOG_FILENAME

    def bin_fs_manager(self, bin_number: int) -> BinManager:
        """Get bin file system manager."""
        return BinManager(self.dir() / f"{self.__BIN_DIR_PREFIX}_{bin_number}")


class BinManager:
    """MFB file system manager for one bin."""

    BIN_DIR_PREFIX = "bin"
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
        """Get bin stats YAML file path."""
        return self.dir() / self.BIN_SEQ_NORMCOV_FILENAME
