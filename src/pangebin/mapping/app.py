"""Mapping applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.mapping.filter as map_filter
import pangebin.mapping.input_output as map_io
import pangebin.pblog as common_log
from pangebin.mapping import create

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class BlastArguments:
    """Mapping arguments."""

    QUERY_FASTA_FILE = typer.Argument(
        help="Query FASTA file",
    )

    SUBJECT_FASTA_FILE = typer.Argument(
        help="Subject FASTA file",
    )

    OUTPUT_FILE = typer.Argument(
        help="Output file",
    )


@APP.command()
def blast(
    query_fasta_file: Annotated[Path, BlastArguments.QUERY_FASTA_FILE],
    subject_fasta_file: Annotated[Path, BlastArguments.SUBJECT_FASTA_FILE],
    out_mapping_file: Annotated[Path, BlastArguments.OUTPUT_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Map a query FASTA file to a subject FASTA file."""
    common_log.init_logger(_LOGGER, "Mapping query to subject.", debug)
    create.blast_map(query_fasta_file, subject_fasta_file, out_mapping_file)
    return out_mapping_file


class FilterArguments:
    """Filter arguments."""

    INPUT_SAM = typer.Argument(
        help="Input mapping file",
    )

    FILTERED_SAM = typer.Argument(
        help="Output file",
    )


class FilterOptions:
    """Filter options."""

    __RICH_HELP_PANEL = "Filter options"

    QUERY_FASTA = typer.Option(
        help="Query FASTA file",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    SUBJECT_FASTA = typer.Option(
        help="Subject FASTA file",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    MIN_LENGTH = typer.Option(
        help="Minimum length threshold",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    MIN_PIDENT = typer.Option(
        help="Minimum percent identity threshold (between 0 and 100)",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    MIN_Q_COVERAGE = typer.Option(
        help=(
            "Minimum query coverage threshold (between 0 and 1)"
            " (only if query FASTA file is provided)"
        ),
        rich_help_panel=__RICH_HELP_PANEL,
    )

    MIN_S_COVERAGE = typer.Option(
        help=(
            "Minimum subject coverage threshold (between 0 and 1)"
            " (only if subject FASTA file is provided)"
        ),
        rich_help_panel=__RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--config",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def filter(  # noqa: A001
    input_sam: Annotated[Path, FilterArguments.INPUT_SAM],
    filtered_sam: Annotated[Path | None, FilterArguments.FILTERED_SAM] = None,
    query_fasta: Annotated[Path | None, FilterOptions.QUERY_FASTA] = None,
    subject_fasta: Annotated[Path | None, FilterOptions.SUBJECT_FASTA] = None,
    min_length: Annotated[
        int,
        FilterOptions.MIN_LENGTH,
    ] = map_filter.Config.DEFAULT_MIN_LENGTH,
    min_pident: Annotated[
        float,
        FilterOptions.MIN_PIDENT,
    ] = map_filter.Config.DEFAULT_MIN_PIDENT,
    min_q_cov: Annotated[
        float,
        FilterOptions.MIN_Q_COVERAGE,
    ] = map_filter.Config.DEFAULT_MIN_Q_COV,
    min_s_cov: Annotated[
        float,
        FilterOptions.MIN_S_COVERAGE,
    ] = map_filter.Config.DEFAULT_MIN_S_COV,
    config_file: Annotated[Path | None, FilterOptions.CONFIG_FILE] = None,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Filter a SAM file."""
    common_log.init_logger(_LOGGER, "Filtering SAM file.", debug)
    config = (
        map_filter.Config.from_yaml(config_file)
        if config_file
        else map_filter.Config(
            min_length=min_length,
            min_pident=min_pident,
            min_q_cov=min_q_cov,
            min_s_cov=min_s_cov,
        )
    )
    if filtered_sam is None:
        filtered_sam = input_sam.with_suffix(".filtered.sam")

    map_io.to_filtered_sam_file(
        input_sam,
        filtered_sam,
        config,
        query_fasta=query_fasta,
        subject_fasta=subject_fasta,
    )

    return filtered_sam
