"""Common subprocess module."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

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


class CommandNotFoundError(Exception):
    """Command not found error."""

    def __init__(self, command: str | Path) -> None:
        """Initialize."""
        super().__init__()
        self.__command = command

    def __str__(self) -> str:
        """Return the error message."""
        return f"Command not found: {self.__command}"
