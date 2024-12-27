"""Common logging module."""

from __future__ import annotations

import logging
import sys

__ROOT_APP_NAME = "pangebin"

_LOGGER = logging.getLogger(__ROOT_APP_NAME)
_LOGGER.setLevel(logging.DEBUG)


def format_logger(
    debug: bool,  # noqa: FBT001
) -> None:
    """Format logger."""
    _LOGGER.handlers.clear()
    _LOGGER.filters.clear()
    handler = logging.StreamHandler(sys.stdout)
    if debug:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s %(name)s [%(levelname)s] %(message)s",
    )
    handler.setFormatter(formatter)
    _LOGGER.addHandler(handler)
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
