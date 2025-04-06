"""Common logging module."""

from __future__ import annotations

import logging

import typer
from rich.console import Console
from rich.logging import RichHandler

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
) -> None:
    """Format logger."""
    _LOGGER.handlers.clear()
    _LOGGER.filters.clear()
    rich_handler = RichHandler(console=CONSOLE)
    if debug:
        rich_handler.setLevel(logging.DEBUG)
    else:
        rich_handler.setLevel(logging.INFO)
    rich_formatter = logging.Formatter(_LOGFORMAT_RICH, datefmt="%Y-%m-%dT%H:%M:%S%z")
    rich_handler.setFormatter(rich_formatter)
    _LOGGER.addHandler(rich_handler)
    log_filter = logging.Filter(__ROOT_APP_NAME)
    _LOGGER.addFilter(log_filter)


def init_logger(
    logger: logging.Logger,
    first_info_message: str,
    debug: bool,  # noqa: FBT001
) -> None:
    """Initialize logger."""
    format_logger(debug)
    logger.info(first_info_message)
