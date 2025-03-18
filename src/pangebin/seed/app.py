"""Seed contig applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.logging as common_log
import pangebin.seed.create as seed_create
import pangebin.seed.input_output as seed_io
import pangebin.seed.thresholds.app as seed_thr_app

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")

APP.add_typer(seed_thr_app.APP)


class FromGeneDensityArguments:
    """From gene density arguments."""

    GENE_DENSITY_FILE = typer.Argument(
        help="Gene density file",
    )


class IOOptions:
    """Input/output options."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_FILE = typer.Argument(
        help="Seed sequence output file",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command(name="pos-gd")
def from_positive_gene_densities(
    gene_density_file: Annotated[Path, FromGeneDensityArguments.GENE_DENSITY_FILE],
    output_file: Annotated[Path, IOOptions.OUTPUT_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Extract seed sequences with positive gene density."""
    common_log.init_logger(
        _LOGGER,
        "Extracting seed sequences with positive gene density.",
        debug,
    )
    seed_io.to_tsv(seed_create.from_gene_density(gene_density_file), output_file)
    _LOGGER.info("Write seed sequence identifiers in file: %s", output_file)
    return output_file
