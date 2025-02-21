"""Ground truth applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.entrez as pg_entrez
import pangebin.ground_truth.config as gt_config
import pangebin.ground_truth.create as gt_create
import pangebin.ground_truth.input_output as gt_io
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


class GroundTruthIOOptions:
    """Input/Output options."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        "-o",
        "--outdir",
        help="Output folder",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def create(
    contigs_fasta_file: Annotated[Path, GroundTruthArguments.CONTIGS_FASTA],
    plasmid_genbank_ids: Annotated[list[str], GroundTruthArguments.PLASMID_GENBANK_IDS],
    min_pident: Annotated[
        float,
        GroundTruthOptions.MIN_PIDENT,
    ] = gt_config.Config.DEFAULT_MIN_PIDENT,
    min_contig_coverage: Annotated[
        float,
        GroundTruthOptions.MIN_CONTIG_COVERAGE,
    ] = gt_config.Config.DEFAULT_MIN_CONTIG_COVERAGE,
    config_file: Annotated[Path | None, GroundTruthOptions.CONFIG_FILE] = None,
    # Entrez configuration
    entrez_email: Annotated[
        str | None,
        pg_entrez.AppOptions.EMAIL,
    ] = pg_entrez.Config.DEFAULT_EMAIL,
    entrez_tool: Annotated[
        str,
        pg_entrez.AppOptions.TOOL,
    ] = pg_entrez.Config.DEFAULT_TOOL,
    entrez_api_key: Annotated[
        str | None,
        pg_entrez.AppOptions.API_KEY,
    ] = pg_entrez.Config.DEFAULT_API_KEY,
    entrez_max_tries: Annotated[
        int,
        pg_entrez.AppOptions.MAX_TRIES,
    ] = pg_entrez.Config.DEFAULT_MAX_TRIES,
    entrez_sleep_between_tries: Annotated[
        int,
        pg_entrez.AppOptions.SLEEP_BETWEEN_TRIES,
    ] = pg_entrez.Config.DEFAULT_SLEEP_BETWEEN_TRIES,
    entrez_config_file: Annotated[
        Path | None,
        pg_entrez.AppOptions.CONFIG_FILE,
    ] = None,
    output_dir: Annotated[
        Path,
        GroundTruthIOOptions.OUTPUT_DIR,
    ] = gt_io.Config.DEFAULT_OUTPUT_DIR,
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
    entrez_config = (
        pg_entrez.Config.from_yaml(entrez_config_file)
        if entrez_config_file
        else pg_entrez.Config(
            email=entrez_email,
            tool=entrez_tool,
            api_key=entrez_api_key,
            max_tries=entrez_max_tries,
            sleep_between_tries=entrez_sleep_between_tries,
        )
    )

    io_manager = gt_io.Manager(gt_io.Config(output_directory=output_dir))
    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)

    plasmid_contigs, non_plasmid_contigs = gt_create.separate_plasmid_and_other_contigs(
        contigs_fasta_file,
        plasmid_genbank_ids,
        config=config,
        entrez_config=entrez_config,
        merged_plasmid_fasta=io_manager.merged_plasmid_fasta(),
        mapping_sam=io_manager.mapping_sam(),
        filtered_mapping_sam=io_manager.filtered_mapping_sam(),
    )

    gt_io.plasmid_contigs_to_file(plasmid_contigs, io_manager.plasmid_contigs_file())
    gt_io.non_plasmid_contigs_to_file(
        non_plasmid_contigs,
        io_manager.non_plasmid_contigs_file(),
    )
