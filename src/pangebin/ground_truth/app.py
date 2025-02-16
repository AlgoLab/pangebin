"""Ground truth applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.ground_truth.config as gt_config
import pangebin.ground_truth.create as gt_create
import pangebin.logging as common_log

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class GroundTruthArguments:
    """Ground truth arguments."""

    CONTIGS_FASTA = typer.Argument(
        help="Contigs FASTA file",
    )
    PLASMID_GENBANK_IDS = typer.Argument(
        help="GenBank IDs of plasmid sequences",
    )
    OUTPUT_DIR = typer.Argument(
        help="Output directory",
    )


class GroundTruthOptions:
    """Ground truth options."""

    __RICH_HELP_PANEL = "Ground truth options"

    MIN_PIDENT = typer.Option(
        help="Minimum percent identity threshold (between 0 and 100)",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    MIN_CONTIG_COVERAGE = typer.Option(
        help="Minimum contig coverage threshold (between 0 and 1)",
        rich_help_panel=__RICH_HELP_PANEL,
    )
    CONFIG_FILE = typer.Option(
        "--config",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    EMAIL_ADDRESS = typer.Option(
        help="Email address to fetch NCBI database",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def create(
    contigs_fasta_file: Annotated[Path, GroundTruthArguments.CONTIGS_FASTA],
    plasmid_genbank_ids: Annotated[list[str], GroundTruthArguments.PLASMID_GENBANK_IDS],
    output_dir: Annotated[Path, GroundTruthArguments.OUTPUT_DIR],
    min_pident: Annotated[
        float,
        GroundTruthOptions.MIN_PIDENT,
    ] = gt_config.Config.DEFAULT_MIN_PIDENT,
    min_contig_coverage: Annotated[
        float,
        GroundTruthOptions.MIN_CONTIG_COVERAGE,
    ] = gt_config.Config.DEFAULT_MIN_CONTIG_COVERAGE,
    email_address: Annotated[str | None, GroundTruthOptions.EMAIL_ADDRESS] = None,
    config_file: Annotated[Path | None, GroundTruthOptions.CONFIG_FILE] = None,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Create ground truth."""
    common_log.init_logger(_LOGGER, "Creating ground truth.", debug)
    config = (
        gt_config.Config.from_yaml(config_file)
        if config_file
        else gt_config.Config(
            min_pident=min_pident,
            min_contig_coverage=min_contig_coverage,
        )
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    gt_create.separate_plasmid_and_other_contigs(
        contigs_fasta_file,
        plasmid_genbank_ids,
        output_dir,
        config=config,
        email_address=email_address,
    )
