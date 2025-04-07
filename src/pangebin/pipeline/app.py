"""Pangebin main pipeline application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from typing import Annotated

import typer

import pangebin.logging as common_log
import pangebin.pipeline.asm_pbf.app as pipe_asm_pbf_app
import pangebin.pipeline.seed_thresholds.app as pipe_seed_thr_app

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


@APP.command()
def run(debug: Annotated[bool, common_log.OPT_DEBUG] = False) -> None:
    """Run the main PangeBin pipeline."""
    common_log.init_logger(_LOGGER, "Running pangebin pipeline.", debug)


CONFIG_APP = typer.Typer(
    name="configs",
    rich_markup_mode="rich",
    help="Write default configuration files for the pipelines",
)


CONFIG_APP.command(name="seed-thresholds")(pipe_seed_thr_app.write_configs)
CONFIG_APP.add_typer(pipe_asm_pbf_app.CONFIGS_APP)
