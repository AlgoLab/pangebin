"""Seed contig thresholds applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.seed.thresholds.config as seed_thr_config
import pangebin.seed.thresholds.create as seed_thr_create
import pangebin.seed.thresholds.input_output as seed_thr_io
from pangebin import pblog

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
    """Seed threshold arguments."""

    DATATEST = typer.Argument(
        help="Datatest file",
    )


class Options:
    """Seed threshold options."""

    __RICH_HELP_PANEL = "Seed threshold options"

    MIN_LENGTH = typer.Option(
        help="Minimum length threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    MAX_LENGTH = typer.Option(
        help="Maximum length threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    STEP_LENGTH = typer.Option(
        help="Step length threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    MIN_GENE_DENSITY = typer.Option(
        help="Minimum gene density threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    MAX_GENE_DENSITY = typer.Option(
        help="Maximum gene density threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    STEP_GENE_DENSITY = typer.Option(
        help="Step gene density threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    CONFIG_FILE = typer.Option(
        "--config",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )


class IOOptions:
    """Input/Output options."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        "-o",
        "--outdir",
        help="Output folder",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def thresholds(
    datatest: Annotated[
        Path,
        Arguments.DATATEST,
    ],
    min_length: Annotated[
        int,
        Options.MIN_LENGTH,
    ] = seed_thr_config.ThresholdRanges.DEFAULT_MIN_LENGTH,
    max_length: Annotated[
        int,
        Options.MAX_LENGTH,
    ] = seed_thr_config.ThresholdRanges.DEFAULT_MAX_LENGTH,
    step_length: Annotated[
        int,
        Options.STEP_LENGTH,
    ] = seed_thr_config.ThresholdRanges.DEFAULT_STEP_LENGTH,
    min_gene_density: Annotated[
        float,
        Options.MIN_GENE_DENSITY,
    ] = seed_thr_config.ThresholdRanges.DEFAULT_MIN_GENE_DENSITY,
    max_gene_density: Annotated[
        float,
        Options.MAX_GENE_DENSITY,
    ] = seed_thr_config.ThresholdRanges.DEFAULT_MAX_GENE_DENSITY,
    step_gene_density: Annotated[
        float,
        Options.STEP_GENE_DENSITY,
    ] = seed_thr_config.ThresholdRanges.DEFAULT_STEP_GENE_DENSITY,
    config_file: Annotated[
        Path | None,
        Options.CONFIG_FILE,
    ] = None,
    output_dir: Annotated[
        Path,
        IOOptions.OUTPUT_DIR,
    ] = seed_thr_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> seed_thr_io.Manager:
    """Seed thresholds."""
    pblog.init_logger(_LOGGER, "Creating seed thresholds.", debug)

    threshold_ranges = (
        seed_thr_config.ThresholdRanges.from_yaml(config_file)
        if config_file is not None
        else seed_thr_config.ThresholdRanges(
            min_length=min_length,
            max_length=max_length,
            step_length=step_length,
            min_gene_density=min_gene_density,
            max_gene_density=max_gene_density,
            step_gene_density=step_gene_density,
        )
    )

    io_manager = seed_thr_io.Manager(
        seed_thr_io.Config(output_directory=output_dir),
    )
    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)

    seed_contig_thresholds = seed_thr_create.from_datatest(
        datatest,
        threshold_ranges,
    )

    seed_contig_thresholds.to_yaml(io_manager.seed_threshold_yaml())

    return io_manager
