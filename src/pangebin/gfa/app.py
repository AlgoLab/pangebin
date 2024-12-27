"""GFA application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import gfapy  # type: ignore[import-untyped]
import typer

import pangebin.gfa.ops as gfa_ops
import pangebin.logging as common_log
from pangebin.gfa import iter as gfa_iter

APP = typer.Typer(rich_markup_mode="rich")

_LOGGER = logging.getLogger(__name__)


@dataclass
class CheckSkesaGFAArgs:
    """Check Skeza GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    OPT_DEBUG = typer.Option(
        help="Debug mode",
    )


@dataclass
class FixSkesaGFAArgs:
    """Fix Skeza GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_OUT_GFA = typer.Argument(
        help="Output GFA file, must be different from input if provided",
    )

    OPT_DEBUG = typer.Option(
        help="Debug mode",
    )


@APP.command()
def check_skesa(
    in_gfa: Annotated[Path, CheckSkesaGFAArgs.ARG_IN_GFA],
    debug: Annotated[bool, CheckSkesaGFAArgs.OPT_DEBUG] = False,
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


@APP.command()
def fix_skesa(
    in_gfa: Annotated[Path, FixSkesaGFAArgs.ARG_IN_GFA],
    out_gfa: Annotated[Path | None, FixSkesaGFAArgs.ARG_OUT_GFA] = None,
    debug: Annotated[bool, CheckSkesaGFAArgs.OPT_DEBUG] = False,
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
        gfa_ops.fix_skesa_gfa(in_gfa, out_gfa_path=out_gfa)
    except ValueError as e:
        raise typer.Exit(1) from e

    _LOGGER.info("Fixed Skeza GFA file: %s", out_gfa if out_gfa is not None else in_gfa)


@dataclass
class ToFASTAArgs:
    """GFA to FASTA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    OPT_DEBUG = typer.Option(
        help="Debug mode",
    )


@APP.command()
def to_fasta(
    gfa_path: Annotated[Path, ToFASTAArgs.ARG_IN_GFA],
    debug: Annotated[bool, ToFASTAArgs.OPT_DEBUG] = False,
) -> None:
    """Convert GFA to FASTA."""
    common_log.init_logger(_LOGGER, "Converting GFA to FASTA.", debug)

    if not gfa_path.exists():
        _LOGGER.error("Input GFA file does not exist: %s", gfa_path)
        raise typer.Exit(1)

    gfa = gfapy.Gfa.from_file(gfa_path)

    for seq_record in gfa_iter.iter_gfa_to_fasta(gfa):
        sys.stdout.write(seq_record.format("fasta"))
