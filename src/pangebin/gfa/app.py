"""GFA application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gfa.header as gfa_header
import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.gfa.segment as gfa_segment
import pangebin.gfa.views as gfa_views
import pangebin.input_output as root_io
import pangebin.pblog as common_log
from pangebin.gfa import iter as gfa_iter

APP = typer.Typer(rich_markup_mode="rich")

_LOGGER = logging.getLogger(__name__)


class CheckSKESAGFAArgs:
    """Check Skeza GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )


@APP.command()
def check_skesa(
    in_gfa: Annotated[Path, CheckSKESAGFAArgs.ARG_IN_GFA],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> bool:
    """Check a Skeza GFA file."""
    common_log.init_logger(_LOGGER, "Checking Skeza GFA file.", debug)

    if not in_gfa.exists():
        _LOGGER.error("Input GFA file does not exist: %s", in_gfa)
        raise typer.Exit(1)

    if gfa_ops.is_skesa_gfa_fixed(in_gfa):
        _LOGGER.info("SKESA GFA file is already fixed.")
        return True
    _LOGGER.info("SKESA GFA file is not fixed.")
    return False


class FixSKESAGFAArgs:
    """Fix SKESA GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_OUT_GFA = typer.Argument(
        help="Output GFA file, must be different from input if provided",
    )


@APP.command()
def fix_skesa(
    in_gfa: Annotated[Path, FixSKESAGFAArgs.ARG_IN_GFA],
    out_gfa: Annotated[Path | None, FixSKESAGFAArgs.ARG_OUT_GFA] = None,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Fix a SKESA GFA file."""
    common_log.init_logger(_LOGGER, "Fixing SKESA GFA file.", debug)

    if not in_gfa.exists():
        _LOGGER.error("Input GFA file does not exist: %s", in_gfa)
        raise typer.Exit(1)

    if gfa_ops.is_skesa_gfa_fixed(in_gfa):
        _LOGGER.info("SKESA GFA file is already fixed.")
        return
    try:
        out_gfa = gfa_ops.fix_skesa_gfa(in_gfa, out_gfa_path=out_gfa)
    except ValueError as e:
        raise typer.Exit(1) from e

    _LOGGER.info("Fixed SKESA GFA file: %s", out_gfa)


class FixUnicyclerGFAArgs:
    """Fix Unicycler GFA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_OUT_GFA = typer.Argument(
        help="Output GFA file, must be different from input if provided",
    )


@APP.command()
def fix_unicycler(
    in_gfa: Annotated[Path, FixUnicyclerGFAArgs.ARG_IN_GFA],
    out_gfa: Annotated[Path | None, FixUnicyclerGFAArgs.ARG_OUT_GFA] = None,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Fix a Unicycler GFA file."""
    common_log.init_logger(_LOGGER, "Fixing Unicycler GFA file.", debug)

    if not in_gfa.exists():
        _LOGGER.error("Input GFA file does not exist: %s", in_gfa)
        raise typer.Exit(1)

    if gfa_ops.is_unicycler_gfa_fixed(in_gfa):
        _LOGGER.info("Unicycler GFA file is already fixed.")
        return
    try:
        out_gfa = gfa_ops.fix_unicycler_gfa(in_gfa, out_gfa_path=out_gfa)
    except ValueError as e:
        raise typer.Exit(1) from e

    _LOGGER.info("Fixed Unicycler GFA file: %s", out_gfa)


class ToFASTAArguments:
    """GFA to FASTA arguments."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_OUT_FASTA = typer.Argument(
        help="Output FASTA file, must be different from input if provided",
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
    fasta_path: Annotated[Path, ToFASTAArguments.ARG_OUT_FASTA],
    attribute_string_separator: Annotated[
        str,
        ToFASTAOptions.OPT_ATTRIBUTE_STRING_SEPARATOR,
    ] = gfa_segment.DEFAULT_ATTRIBUTE_STR_SEP,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Convert GFA to FASTA (write to stdout)."""
    common_log.init_logger(_LOGGER, "Converting GFA to FASTA.", debug)

    if not gfa_path.exists():
        _LOGGER.error("Input GFA file does not exist: %s", gfa_path)
        raise typer.Exit(1)

    with root_io.open_file_write(fasta_path) as f_out:
        for seq_record in gfa_iter.sequence_records(
            gfa_path,
            sep=attribute_string_separator,
        ):
            f_out.write(seq_record.format("fasta"))

    _LOGGER.info("Resulting FASTA file: %s", fasta_path)


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


class PrintStatsArgs:
    """Argument for printing GFA stats."""

    GFA_FILE = typer.Argument(
        help="GFA file path",
    )


@APP.command()
def stats(
    gfa_path: Annotated[Path, PrintStatsArgs.GFA_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Print GFA stats."""
    common_log.init_logger(_LOGGER, "Printing GFA stats.", debug)

    if not gfa_path.exists():
        _LOGGER.error("Input GFA file does not exist: %s", gfa_path)
        raise typer.Exit(1)

    gfa_views.print_stats(gfa_path)


class SubRadiusArgs:
    """Argument for extracting subgraph with neighbor radius."""

    GFA_FILE = typer.Argument(
        help="GFA file path",
    )

    SEGMENTS_FILE = typer.Argument(
        help="Segments serving as centers",
    )

    RADIUS = typer.Argument(
        help="Radius",
    )

    SUB_GFA_FILE = typer.Argument(
        help="Output subgraph GFA file, must be different from input if provided",
    )


@APP.command()
def sub_radius(
    gfa_path: Annotated[Path, SubRadiusArgs.GFA_FILE],
    sub_gfa_path: Annotated[Path, SubRadiusArgs.SUB_GFA_FILE],
    radius: Annotated[int, SubRadiusArgs.RADIUS],
    segments_paths: Annotated[list[Path], SubRadiusArgs.SEGMENTS_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Extract subgraph with neighbor radius."""
    common_log.init_logger(_LOGGER, "Extracting subgraph with neighbor radius.", debug)

    if not gfa_path.exists():
        _LOGGER.critical("Input GFA file does not exist: %s", gfa_path)
        raise typer.Exit(1)

    if radius < 0:
        _LOGGER.critical("Radius must be positive.")
        raise typer.Exit(1)

    for p in segments_paths:
        if not p.exists():
            _LOGGER.critical("Input segments file does not exist: %s", p)
            raise typer.Exit(1)

    gfa_graph = gfa_io.from_file(gfa_path)

    centers: list[str] = []
    for p in segments_paths:
        with p.open() as f_in:
            centers.extend([line.strip() for line in f_in])

    sub_graph = gfa_ops.sub_radius_graph(gfa_graph, centers, radius)

    with root_io.open_file_write(sub_gfa_path) as f_out:
        for line in sub_graph.lines:
            f_out.write(f"{line}\n")


class RemoveSmallSequencesArgs:
    """Argument for removing small sequences."""

    ARG_IN_GFA = typer.Argument(
        help="Input GFA file",
    )

    ARG_OUT_GFA = typer.Argument(
        help="Output GFA file, must be different from input if provided",
    )


class RemoveSmallSequencesOpts:
    """Options for removing small sequences."""

    MIN_LENGTH_DEF = 100
    MIN_LENGTH = typer.Option(
        "--min-length",
        "-m",
        help="Minimum length threshold (threshold kept)",
    )


@APP.command("rm-small-seq")
def remove_small_sequences(
    in_gfa: Annotated[Path, RemoveSmallSequencesArgs.ARG_IN_GFA],
    out_gfa: Annotated[Path | None, RemoveSmallSequencesArgs.ARG_OUT_GFA] = None,
    min_length: Annotated[
        int,
        RemoveSmallSequencesOpts.MIN_LENGTH,
    ] = RemoveSmallSequencesOpts.MIN_LENGTH_DEF,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Remove small sequences by preserving the walks."""
    common_log.init_logger(_LOGGER, "Removing small sequences in a GFA.", debug)

    if not in_gfa.exists():
        _LOGGER.error("Input GFA file does not exist: %s", in_gfa)
        raise typer.Exit(1)

    if out_gfa is None:
        out_gfa = in_gfa
    elif out_gfa == in_gfa:
        _LOGGER.error("Output GFA file must be different from input if provided.")
        raise typer.Exit(1)

    graph = gfa_io.from_file(in_gfa)

    gfa_ops.transform_small_contigs_into_links(graph, min_length)

    with root_io.open_file_write(out_gfa) as f_out:
        for line in graph.lines:
            f_out.write(f"{line}\n")

    if out_gfa == in_gfa:
        _LOGGER.info("Inplace filtered GFA file: %s", out_gfa)
    else:
        _LOGGER.info("New filtered GFA file: %s", out_gfa)
