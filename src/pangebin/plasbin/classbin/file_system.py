"""Plasbin classbin MILP file system logics."""

from __future__ import annotations

from pathlib import Path

from .multi_flow.milp import file_system as mfb_fs


class Manager:
    """MILP file system manager."""

    MILP_DIRNAME = Path("milp")

    def __init__(self, root_directory: Path) -> None:
        """Initialize object."""
        self.__root_directory = root_directory
        self.__classify_fs_manager = ClassifyManager(self.milp_dir())
        # REFACTOR may be simplified
        self.__mfb_fs_manager = mfb_fs.Manager(self.milp_dir())

    def root_dir(self) -> Path:
        """Get root directory path."""
        return self.__root_directory

    def milp_dir(self) -> Path:
        """Get MILP directory path."""
        return self.__root_directory / self.MILP_DIRNAME

    def classify_fs_manager(self) -> ClassifyManager:
        """Get classify file system manager."""
        return self.__classify_fs_manager

    def mfb_fs_manager(self) -> mfb_fs.Manager:
        """Get multi-flow binning MILP file system manager."""
        return self.__mfb_fs_manager


# REFACTOR potentially depreciated
class ClassifyManager:
    """Classify MILP file system manager."""

    CLASSIFY_DIRNAME = Path("classify")

    def __init__(self, root_directory: Path) -> None:
        """Initialize object."""
        self.__root_directory = root_directory

    def classify_dir(self) -> Path:
        """Get classify directory path."""
        return self.__root_directory / self.CLASSIFY_DIRNAME

    def gurobi_log_path(self) -> Path:
        """Get Gurobi log file path."""
        return self.classify_dir() / "gurobi.log"
