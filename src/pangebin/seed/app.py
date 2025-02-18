"""Seed sequences applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.logging as common_log
import pangebin.seed.config as seed_config
import pangebin.seed.create as seed_create
import pangebin.seed.input_output as seed_io

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class SeedThresholdsArguments:
    """Seed threshold arguments."""

    DATATEST = typer.Argument(
        help="Datatest file",
    )


class SeedThresholdsOptions:
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


class SeedThresholdsIOOptions:
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
        SeedThresholdsArguments.DATATEST,
    ],
    min_length: Annotated[
        int,
        SeedThresholdsOptions.MIN_LENGTH,
    ] = seed_config.ThresholdRanges.DEFAULT_MIN_LENGTH,
    max_length: Annotated[
        int,
        SeedThresholdsOptions.MAX_LENGTH,
    ] = seed_config.ThresholdRanges.DEFAULT_MAX_LENGTH,
    step_length: Annotated[
        int,
        SeedThresholdsOptions.STEP_LENGTH,
    ] = seed_config.ThresholdRanges.DEFAULT_STEP_LENGTH,
    min_gene_density: Annotated[
        float,
        SeedThresholdsOptions.MIN_GENE_DENSITY,
    ] = seed_config.ThresholdRanges.DEFAULT_MIN_GENE_DENSITY,
    max_gene_density: Annotated[
        float,
        SeedThresholdsOptions.MAX_GENE_DENSITY,
    ] = seed_config.ThresholdRanges.DEFAULT_MAX_GENE_DENSITY,
    step_gene_density: Annotated[
        float,
        SeedThresholdsOptions.STEP_GENE_DENSITY,
    ] = seed_config.ThresholdRanges.DEFAULT_STEP_GENE_DENSITY,
    config_file: Annotated[
        Path | None,
        SeedThresholdsOptions.CONFIG_FILE,
    ] = None,
    output_dir: Annotated[
        Path,
        SeedThresholdsIOOptions.OUTPUT_DIR,
    ] = seed_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Seed thresholds."""
    common_log.init_logger(_LOGGER, "Creating seed thresholds.", debug)

    threshold_ranges = (
        seed_config.ThresholdRanges.from_yaml(config_file)
        if config_file is not None
        else seed_config.ThresholdRanges(
            min_length=min_length,
            max_length=max_length,
            step_length=step_length,
            min_gene_density=min_gene_density,
            max_gene_density=max_gene_density,
            step_gene_density=step_gene_density,
        )
    )

    io_manager = seed_io.ThresholdManager(seed_io.Config(output_directory=output_dir))
    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)

    seed_contig_thresholds = seed_create.thresholds_from_datatest(
        datatest,
        threshold_ranges,
    )

    seed_contig_thresholds.to_yaml(io_manager.seed_threshold_yaml())
