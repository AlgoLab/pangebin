"""Mapping applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.logging as common_log
from pangebin.mapping import create

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
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
    query_fasta_file: Annotated[Path, Arguments.QUERY_FASTA_FILE],
    subject_fasta_file: Annotated[Path, Arguments.SUBJECT_FASTA_FILE],
    out_mapping_file: Annotated[Path, Arguments.OUTPUT_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Map a query FASTA file to a subject FASTA file."""
    common_log.init_logger(_LOGGER, "Mapping query to subject.", debug)
    create.blast_map(query_fasta_file, subject_fasta_file, out_mapping_file)
