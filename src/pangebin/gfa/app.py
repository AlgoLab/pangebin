"""GFA application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gfa.header as gfa_header
import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.gfa.segment as gfa_segment
import pangebin.logging as common_log
from pangebin.gfa import iter as gfa_iter

APP = typer.Typer(rich_markup_mode="rich")

_LOGGER = logging.getLogger(__name__)


class CheckSkesaGFAArgs:
    """Check Skeza GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )


@APP.command()
def check_skesa(
    in_gfa: Annotated[Path, CheckSkesaGFAArgs.ARG_IN_GFA],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> bool:
    """Check a Skeza GFA file."""
    common_log.init_logger(_LOGGER, "Checking Skeza GFA file.", debug)

    if not in_gfa.exists():
        _LOGGER.error("Input GFA file does not exist: %s", in_gfa)
        raise typer.Exit(1)

    if gfa_ops.is_skesa_gfa_fixed(in_gfa):
        _LOGGER.info("Skesa GFA file is already fixed.")
        return True
    _LOGGER.info("Skesa GFA file is not fixed.")
    return False


class FixSkesaGFAArgs:
    """Fix Skeza GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_OUT_GFA = typer.Argument(
        help="Output GFA file, must be different from input if provided",
    )


@APP.command()
def fix_skesa(
    in_gfa: Annotated[Path, FixSkesaGFAArgs.ARG_IN_GFA],
    out_gfa: Annotated[Path | None, FixSkesaGFAArgs.ARG_OUT_GFA] = None,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Fix a Skeza GFA file."""
    common_log.init_logger(_LOGGER, "Fixing Skeza GFA file.", debug)

    if not in_gfa.exists():
        _LOGGER.error("Input GFA file does not exist: %s", in_gfa)
        raise typer.Exit(1)

    if gfa_ops.is_skesa_gfa_fixed(in_gfa):
        _LOGGER.info("Skesa GFA file is already fixed.")
        return
    try:
        out_gfa = gfa_ops.fix_skesa_gfa(in_gfa, out_gfa_path=out_gfa)
    except ValueError as e:
        raise typer.Exit(1) from e

    _LOGGER.info("Fixed Skeza GFA file: %s", out_gfa)


class ToFASTAArguments:
    """GFA to FASTA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )


class ToFASTAOptions:
    """GFA to FASTA options."""

    __RICH_HELP_PANEL = "GFA to FASTA options"

    OPT_ATTRIBUTE_STRING_SEPARATOR = typer.Option(
        help="String separator for attributes",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def to_fasta(
    gfa_path: Annotated[Path, ToFASTAArguments.ARG_IN_GFA],
    attribute_string_separator: Annotated[
        str,
        ToFASTAOptions.OPT_ATTRIBUTE_STRING_SEPARATOR,
    ] = gfa_segment.DEFAULT_ATTRIBUTE_STR_SEP,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Convert GFA to FASTA."""
    common_log.init_logger(_LOGGER, "Converting GFA to FASTA.", debug)

    if not gfa_path.exists():
        _LOGGER.error("Input GFA file does not exist: %s", gfa_path)
        raise typer.Exit(1)

    for seq_record in gfa_iter.sequence_records(
        gfa_path,
        sep=attribute_string_separator,
    ):
        sys.stdout.write(seq_record.format("fasta"))


class ISGFAStandardizeArgs:
    """Argument for checking if a GFA is standardized."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )


@APP.command()
def is_standardized(
    gfa_path: Annotated[Path, ISGFAStandardizeArgs.ARG_IN_GFA],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Check if a GFA is standardized."""
    common_log.init_logger(_LOGGER, "Checking if GFA is standardized.", debug)

    if not gfa_path.exists():
        _LOGGER.error("Input GFA file does not exist: %s", gfa_path)
        raise typer.Exit(1)

    gfa = gfa_io.from_file(gfa_path)

    if gfa_header.is_standardized(gfa):
        _LOGGER.info("GFA is standardized.")
    else:
        _LOGGER.info("GFA is not standardized.")
