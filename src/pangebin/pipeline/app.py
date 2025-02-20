"""Pangebin main pipeline application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from typing import Annotated

import typer

import pangebin.logging as common_log

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


@APP.command()
def run(debug: Annotated[bool, common_log.OPT_DEBUG] = False) -> None:
    """Run the main PangeBin pipeline."""
    common_log.init_logger(_LOGGER, "Running pangebin pipeline.", debug)


@APP.command()
def seed_thresholds(debug: Annotated[bool, common_log.OPT_DEBUG] = False) -> None:
    """Obtain the seed threshold pairs from paired Illumina BioSamples."""
    common_log.init_logger(
        _LOGGER,
        "Obtaining the seed threshold pairs from paired Illumina BioSamples.",
        debug,
    )
