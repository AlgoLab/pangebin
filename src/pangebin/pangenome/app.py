"""Pangenome application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Annotated

import typer

import pangebin.input_output as common_io
import pangebin.pangenome.config as pangenome_config
import pangebin.pangenome.create as pangenome_create
import pangebin.pangenome.input_output as pangenome_io
import pangebin.pblog as common_log

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
    """Standardize arguments."""

    INPUT_MIXED_FASTA = typer.Argument(
        help="Mixed FASTA file (can be gzipped or bgzipped or not)",
    )


class ProcessOpts:
    """Process options."""

    __RICH_HELP_PANEL = "Process options"

    RELEASE = typer.Option(
        help="nf-core/pangenome release",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    PROFILE = typer.Option(
        help="Profile environment",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    RESUME = typer.Option(
        help="Resume workflow",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    OVERIDING_CFG_FILE = typer.Option(
        help="Overiding base nf-core/pangenome configuration file",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--config",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )


class IOOpts:
    """Input/Output options."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        "-o",
        "--outdir",
        help="Output folder",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def pangenome(
    mixed_fasta: Annotated[Path, Arguments.INPUT_MIXED_FASTA],
    # Process options
    release: Annotated[
        str,
        ProcessOpts.RELEASE,
    ] = pangenome_config.Pangenome.DEFAULT_RELEASE,
    profile: Annotated[
        pangenome_config.NextflowProfile,
        ProcessOpts.PROFILE,
    ] = pangenome_config.Pangenome.DEFAULT_PROFILE,
    resume: Annotated[
        bool,
        ProcessOpts.RESUME,
    ] = pangenome_config.Pangenome.DEFAULT_RESUME,
    supplementary_nfcore_pangenome_config_path: Annotated[
        Path | None,
        ProcessOpts.OVERIDING_CFG_FILE,
    ] = pangenome_config.Pangenome.DEFAULT_SUPP_NFCORE_PANGENOME_CONFIG_PATH,
    config_file: Annotated[Path | None, ProcessOpts.CONFIG_FILE] = None,
    # Input/Output options
    outdir: Annotated[Path, IOOpts.OUTPUT_DIR] = pangenome_io.Config.DEFAULT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Produce a pangenome using nf-core/pangenome."""
    common_log.init_logger(_LOGGER, "Producing pangenome.", debug)

    config = (
        pangenome_config.Pangenome.from_yaml(config_file)
        if config_file is not None
        else pangenome_config.Pangenome(
            release=release,
            profile=profile,
            resume=resume,
            supp_nfcore_pangenome_config_path=supplementary_nfcore_pangenome_config_path,
        )
    )

    io_manager = pangenome_io.Manager(
        pangenome_io.Config(
            mixed_fasta_path=mixed_fasta,
            output_directory=outdir,
        ),
    )

    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)
    io_manager.nfcore_pangenome_directory().mkdir(parents=True, exist_ok=True)

    temp_mixed_fasta_gz = False
    if common_io.is_gz_file(mixed_fasta):
        mixed_fasta_gz = mixed_fasta
    else:
        mixed_fasta_gz = common_io.bgzip_file(mixed_fasta)
        temp_mixed_fasta_gz = True

    pangenome_create.add_false_sequence(mixed_fasta_gz)

    pangenome_create.nfcore_pangenome(
        mixed_fasta_gz,
        config,
        io_manager.nfcore_pangenome_directory(),
    )

    if temp_mixed_fasta_gz:
        mixed_fasta_gz.unlink()

    shutil.copy(
        io_manager.nfcore_pangenome_gfa_path(),
        io_manager.pangenome_gfa_path(),
    )

    pangenome_create.rename_gfa_paths(io_manager.pangenome_gfa_path())
