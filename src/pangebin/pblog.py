"""Common logging module."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import typer
from rich.console import Console
from rich.logging import RichHandler

if TYPE_CHECKING:
    from pathlib import Path

__ROOT_APP_NAME = "pangebin"

_LOGFORMAT = "%(asctime)s %(name)s [%(levelname)s] %(message)s"
_LOGFORMAT_RICH = "%(message)s"

_LOGGER = logging.getLogger(__ROOT_APP_NAME)
_LOGGER.setLevel(logging.DEBUG)

CONSOLE = Console()

OPT_DEBUG = typer.Option(
    help="Debug mode",
)


def format_logger(
    debug: bool,  # noqa: FBT001
    log_file: Path | None,
) -> None:
    """Format logger."""
    _LOGGER.handlers.clear()
    _LOGGER.filters.clear()

    log_filter = logging.Filter(__ROOT_APP_NAME)
    _LOGGER.addFilter(log_filter)

    rich_handler = RichHandler(console=CONSOLE)
    if debug:
        rich_handler.setLevel(logging.DEBUG)
    else:
        rich_handler.setLevel(logging.INFO)
    rich_formatter = logging.Formatter(_LOGFORMAT_RICH, datefmt="%Y-%m-%dT%H:%M:%S%z")
    rich_handler.setFormatter(rich_formatter)
    _LOGGER.addHandler(rich_handler)

    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        if debug:
            file_handler.setLevel(logging.DEBUG)
        else:
            file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter(_LOGFORMAT, datefmt="%Y-%m-%dT%H:%M:%S%z")
        file_handler.setFormatter(file_formatter)
        _LOGGER.addHandler(file_handler)


def init_logger(
    logger: logging.Logger,
    first_info_message: str,
    debug: bool,  # noqa: FBT001
    log_file: Path | None = None,
) -> None:
    """Initialize logger."""
    format_logger(debug, log_file)
    logger.info(first_info_message)
