"""Pan-assembly application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913
from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gfa.header as gfa_header
import pangebin.gfa.input_output as gfa_io
import pangebin.logging as common_log
import pangebin.panassembly.input_output as panassembly_io
from pangebin.panassembly.create import pangenome_graph_into_panassembly_graph

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
    """Pangenome assembly arguments."""

    PANGENOME_GFA = typer.Argument(
        help="Pangenome GFA file",
    )

    STANDARDIZED_UNICYCLER_GFA = typer.Argument(
        help="Unicycler GFA standardized graph",
    )

    STANDARDIZED_SKESA_GFA = typer.Argument(
        help="Skesa GFA standardized graph",
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
def panassembly(
    pangenome_gfa_path: Annotated[Path, Arguments.PANGENOME_GFA],
    standardized_skesa_gfa_path: Annotated[Path, Arguments.STANDARDIZED_SKESA_GFA],
    standardized_unicycler_gfa_path: Annotated[
        Path,
        Arguments.STANDARDIZED_UNICYCLER_GFA,
    ],
    outdir: Annotated[Path, IOOpts.OUTPUT_DIR] = panassembly_io.Config.DEFAULT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Produce a pangenome assembly from the original assemblers and the pangenome."""
    common_log.init_logger(_LOGGER, "Producing pangenome assembly.", debug)

    standardized_skesa_gfa = gfa_io.from_file(standardized_skesa_gfa_path)
    if not gfa_header.is_standardized(standardized_skesa_gfa):
        _LOGGER.error("Skesa GFA is not standardized.")
        raise typer.Exit(1)

    standardized_unicycler_gfa = gfa_io.from_file(standardized_unicycler_gfa_path)
    if not gfa_header.is_standardized(standardized_unicycler_gfa):
        _LOGGER.error("Unicycler GFA is not standardized.")
        raise typer.Exit(1)

    outdir.mkdir(parents=True, exist_ok=True)

    io_manager = panassembly_io.Manager(
        panassembly_io.Config(output_directory=outdir),
    )

    pangenome_gfa = gfa_io.from_file(pangenome_gfa_path)
    pangenome_graph_into_panassembly_graph(
        pangenome_gfa,
        standardized_skesa_gfa,
        standardized_unicycler_gfa,
    )

    pangenome_gfa.to_file(io_manager.panassembly_gfa_path())
