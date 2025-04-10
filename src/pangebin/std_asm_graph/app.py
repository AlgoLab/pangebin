"""Standardize application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.assembly.items as asm_items
import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.pblog as common_log
import pangebin.std_asm_graph.fasta as standardize_fasta
import pangebin.std_asm_graph.input_output as standardize_io
from pangebin.std_asm_graph.config import Config
from pangebin.std_asm_graph.ops import (
    standardize_assembly_graph,
)

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
    """Standardize arguments."""

    INPUT_UNICYCLER_GFA = typer.Argument(
        help="Unicycler GFA assembly graph file",
    )

    INPUT_SKESA_GFA = typer.Argument(
        help="SKESA GFA assembly graph file",
    )


class ConfigOpts:
    """Computationnal options."""

    __RICH_HELP_PANEL = "Configuration options"

    MIN_CONTIG_LENGTH = typer.Option(
        help="Minimum contig length threshold",
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
def std_asm_graph(
    skesa_gfa_path: Annotated[Path, Arguments.INPUT_SKESA_GFA],
    unicycler_gfa_path: Annotated[Path, Arguments.INPUT_UNICYCLER_GFA],
    # Config
    min_contig_length: Annotated[
        int,
        ConfigOpts.MIN_CONTIG_LENGTH,
    ] = Config.DEFAULT_MIN_CONTIG_LENGTH,
    config_file: Annotated[Path | None, ConfigOpts.CONFIG_FILE] = None,
    # IO options
    outdir: Annotated[Path, IOOpts.OUTPUT_DIR] = standardize_io.Config.DEFAULT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Standardize GFA assembly graphs."""
    common_log.init_logger(_LOGGER, "Standardizing GFA Assembly graphs.", debug)

    outdir.mkdir(parents=True, exist_ok=True)

    config = (
        Config.from_yaml(config_file)
        if config_file is not None
        else Config(min_contig_length)
    )
    _LOGGER.debug("Config: %s", config)

    io_manager = standardize_io.Manager(
        standardize_io.Config(output_directory=outdir),
    )

    if not gfa_ops.is_skesa_gfa_fixed(skesa_gfa_path):
        skesa_gfa_path = gfa_ops.fix_skesa_gfa(skesa_gfa_path)

    skesa_gfa = gfa_io.from_file(skesa_gfa_path)
    unicyler_gfa = gfa_io.from_file(unicycler_gfa_path)

    standardize_assembly_graph(skesa_gfa, asm_items.Identifier.SKESA, config)
    standardize_assembly_graph(unicyler_gfa, asm_items.Identifier.UNICYCLER, config)

    skesa_gfa.to_file(io_manager.skesa_gfa_path())
    unicyler_gfa.to_file(io_manager.unicycler_gfa_path())

    gfa_io.gfa_file_to_fasta_file(
        io_manager.skesa_gfa_path(),
        io_manager.skesa_fasta_path(),
    )
    gfa_io.gfa_file_to_fasta_file(
        io_manager.unicycler_gfa_path(),
        io_manager.unicycler_fasta_path(),
    )

    standardize_fasta.fastas_to_pansn_mixed_fasta(
        io_manager.skesa_fasta_path(),
        io_manager.unicycler_fasta_path(),
        io_manager.mixed_fasta_path(),
    )
