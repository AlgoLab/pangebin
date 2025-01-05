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
from pangebin import fasta
from pangebin.assembler import ContigPrefix
from pangebin.preprocess.config import Config
from pangebin.preprocess.create import transform_small_contigs_into_links
from pangebin.preprocess.input_output import PreprocessIOConfig, PreprocessIOManager

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


@dataclass
class PreprocessArgs:
    """GFAUtils arguments."""

    __IO_CAT = "I/O options"
    __CONFIG_CAT = "Config options"

    ARG_INPUT_UNICYCLER_GFA = typer.Argument(
        help="Unicler GFA assembly graph file",
    )

    ARG_INPUT_SKESA_GFA = typer.Argument(
        help="Skesa GFA assembly graph file",
    )

    OPT_MIN_CONTIG_LENGTH = typer.Option(
        help="Minimum contig length threshold",
        rich_help_panel=__CONFIG_CAT,
    )

    OPT_CONFIG_FILE = typer.Option(
        "--config",
        help="The config filepath",
        rich_help_panel=__CONFIG_CAT,
    )

    OPT_OUTPUT_DIR = typer.Option(
        "-o",
        "--outdir",
        # DOCU: here we will place unicycler.gfa skesa.gfa, unicycler.fasta, skesa.fasta, mixed.fasta.gz
        help="Output folder",
        rich_help_panel=__IO_CAT,
    )

    OPT_DEBUG = typer.Option(
        help="Debug mode",
    )


@APP.command()
def preprocess(
    unicycler_gfa_path: Annotated[Path, PreprocessArgs.ARG_INPUT_UNICYCLER_GFA],
    skesa_gfa_path: Annotated[Path, PreprocessArgs.ARG_INPUT_SKESA_GFA],
    min_contig_length: Annotated[
        int,
        PreprocessArgs.OPT_MIN_CONTIG_LENGTH,
    ] = Config.DEFAULT_MIN_CONTIG_LENGTH,
    config_file: Annotated[
        Path | None,
        PreprocessArgs.OPT_CONFIG_FILE,
    ] = None,
    outdir: Annotated[
        Path,
        PreprocessArgs.OPT_OUTPUT_DIR,
    ] = PreprocessIOConfig.DEFAULT_DIR,
    debug: Annotated[bool, PreprocessArgs.OPT_DEBUG] = False,
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

    preprocess_io_manager = PreprocessIOManager(
        PreprocessIOConfig(output_directory=outdir),
    )

    unicyler_gfa = gfa_io.from_file(unicycler_gfa_path)

    if not gfa_ops.is_skesa_gfa_fixed(skesa_gfa_path):
        gfa_ops.fix_skesa_gfa(skesa_gfa_path)

    skesa_gfa = gfa_io.from_file(skesa_gfa_path)

    gfa_ops.rename_contigs(unicyler_gfa, ContigPrefix.UNICYCLER)
    gfa_ops.rename_contigs(skesa_gfa, ContigPrefix.SKESA)
    gfa_ops.convert_kc_to_dp(skesa_gfa)

    transform_small_contigs_into_links(unicyler_gfa, config.min_contig_length())
    transform_small_contigs_into_links(skesa_gfa, config.min_contig_length())

    unicyler_gfa.to_file(preprocess_io_manager.unicycler_gfa_path())
    skesa_gfa.to_file(preprocess_io_manager.skesa_gfa_path())

    gfa_io.gfa_to_fasta_file(unicyler_gfa, preprocess_io_manager.unicycler_fasta_path())
    gfa_io.gfa_to_fasta_file(skesa_gfa, preprocess_io_manager.skesa_fasta_path())

    fasta.merge_two_fastas(
        preprocess_io_manager.unicycler_fasta_path(),
        preprocess_io_manager.skesa_fasta_path(),
        preprocess_io_manager.mixed_fasta_path(),
    )
