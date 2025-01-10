"""Preprocess application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.logging as common_log
import pangebin.std_asm_graph.input_output as preprocess_io
from pangebin import assembler, fasta
from pangebin.std_asm_graph.config import Config
from pangebin.std_asm_graph.ops import transform_small_contigs_into_links

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


@dataclass
class Arguments:
    """Preprocess arguments."""

    INPUT_UNICYCLER_GFA = typer.Argument(
        help="Unicycler GFA assembly graph file",
    )

    INPUT_SKESA_GFA = typer.Argument(
        help="Skesa GFA assembly graph file",
    )


@dataclass
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
    unicycler_gfa_path: Annotated[Path, Arguments.INPUT_UNICYCLER_GFA],
    skesa_gfa_path: Annotated[Path, Arguments.INPUT_SKESA_GFA],
    min_contig_length: Annotated[
        int,
        ConfigOpts.MIN_CONTIG_LENGTH,
    ] = Config.DEFAULT_MIN_CONTIG_LENGTH,
    config_file: Annotated[Path | None, ConfigOpts.CONFIG_FILE] = None,
    outdir: Annotated[Path, IOOpts.OUTPUT_DIR] = preprocess_io.Config.DEFAULT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Preprocess GFA Assembly files."""
    common_log.init_logger(_LOGGER, "Preprocessing GFA Assembly files.", debug)

    outdir.mkdir(parents=True, exist_ok=True)

    config = (
        Config.from_yaml(config_file)
        if config_file is not None
        else Config(
            min_contig_length=min_contig_length,
        )
    )

    io_manager = preprocess_io.Manager(
        preprocess_io.Config(output_directory=outdir),
    )

    unicyler_gfa = gfa_io.from_file(unicycler_gfa_path)

    if not gfa_ops.is_skesa_gfa_fixed(skesa_gfa_path):
        gfa_ops.fix_skesa_gfa(skesa_gfa_path)

    skesa_gfa = gfa_io.from_file(skesa_gfa_path)

    gfa_ops.set_segment_length_tags(unicyler_gfa)
    gfa_ops.set_segment_length_tags(skesa_gfa)

    gfa_ops.rename_contigs(unicyler_gfa, assembler.Prefix.UNICYCLER)
    gfa_ops.rename_contigs(skesa_gfa, assembler.Prefix.SKESA)

    gfa_ops.convert_kmer_coverage_to_normalized_coverage(skesa_gfa)

    transform_small_contigs_into_links(unicyler_gfa, config.min_contig_length())
    transform_small_contigs_into_links(skesa_gfa, config.min_contig_length())

    gfa_ops.set_preprocessed_header_tag(unicyler_gfa)
    gfa_ops.set_preprocessed_header_tag(skesa_gfa)

    unicyler_gfa.to_file(io_manager.unicycler_gfa_path())
    skesa_gfa.to_file(io_manager.skesa_gfa_path())

    gfa_io.gfa_to_fasta_file(unicyler_gfa, io_manager.unicycler_fasta_path())
    gfa_io.gfa_to_fasta_file(skesa_gfa, io_manager.skesa_fasta_path())

    fasta.merge_two_fastas(
        io_manager.unicycler_fasta_path(),
        io_manager.skesa_fasta_path(),
        io_manager.mixed_fasta_path(),
    )
