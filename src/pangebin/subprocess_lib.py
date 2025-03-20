"""Common subprocess module."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pangebin.yaml import YAMLInterface

if TYPE_CHECKING:
    from collections.abc import Sequence

_LOGGER = logging.getLogger(__name__)


def command_path(command_str: str | Path) -> Path:
    """Get command path.

    Raises
    ------
    CommandNotFoundError
        If command not found.

    """
    cmd_path = shutil.which(command_str)
    if cmd_path is None:
        _LOGGER.critical("Command not found: %s", command_str)
        raise CommandNotFoundError(command_str)
    return Path(cmd_path)


def run_cmd(cli_line: Sequence[object], cmd_str: str) -> None:
    """Run external command."""
    try:
        subprocess.run(  # noqa: S603
            [str(x) for x in cli_line],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _cmd_err = CommandFailedError(cmd_str, exc)
        _LOGGER.critical(str(_cmd_err))
        raise _cmd_err from exc


class CommandNotFoundError(Exception):
    """Command not found error."""

    def __init__(self, command: str | Path) -> None:
        """Initialize."""
        super().__init__()
        self.__command = command

    def __str__(self) -> str:
        """Return the error message."""
        return f"Command not found: {self.__command}"


class CommandFailedError(Exception):
    """Command failed error."""

    def __init__(
        self,
        cmd_str: str,
        called_proc_exc: subprocess.CalledProcessError,
    ) -> None:
        """Initialize."""
        super().__init__()
        self.__cmd_str = cmd_str
        self.__called_proc_exc = called_proc_exc

    def cmd_str(self) -> str:
        """Return the command string."""
        return self.__cmd_str

    def called_proc_exc(self) -> subprocess.CalledProcessError:
        """Return the command."""
        return self.__called_proc_exc

    def __str__(self) -> str:
        """Return the error message."""
        return f"{self.__cmd_str} command failed: {self.__called_proc_exc.stderr}"


class RessourcesConfig(YAMLInterface):
    """Ressources config."""

    DEFAULT_MAX_CORES = 8
    DEFAULT_MAX_MEMORY = 8

    KEY_MAX_CORES = "max_number_of_cores"
    KEY_MAX_MEMORY = "max_memory"

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> RessourcesConfig:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_MAX_CORES, cls.DEFAULT_MAX_CORES),
        )

    def __init__(
        self,
        max_cores: int = DEFAULT_MAX_CORES,
        max_memory: int = DEFAULT_MAX_MEMORY,
    ) -> None:
        """Initialize object.

        Parameters
        ----------
        max_cores : int, optional
            Max number of cores, by default DEFAULT_MAX_CORES
        max_memory : int, optional
            Max memory usage (in GB), by default DEFAULT_MAX_MEMORY
        """
        self.__max_cores = max_cores
        self.__max_memory = max_memory

    def max_cores(self) -> int:
        """Get max number of cores option."""
        return self.__max_cores

    def max_memory(self) -> int:
        """Get max memory option (in GB)."""
        return self.__max_memory

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MAX_CORES: self.__max_cores,
            self.KEY_MAX_MEMORY: self.__max_memory,
        }
