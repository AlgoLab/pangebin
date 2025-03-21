"""GC content appllications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.logging as common_log
from pangebin.gc_content import create, items
from pangebin.gc_content import input_output as io

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


DEFAULT_GC_CONTENT_INTERVAL_FILE = (
    Path(__file__).parent / "default_gc_content_intervals.txt"
)


class GFAScoresArguments:
    """GFA GC content arguments."""

    GFA_FILE = typer.Argument(
        help="GFA assembly graph file",
    )

    GC_SCORES_TSV = typer.Argument(
        help="Output TSV file with the sequences and their GC scores",
    )


class ScoresOptions:
    """GC content options."""

    _RICH_HELP_PANEL = "GFA GC content options"

    GC_CONTENT_INTERVAL_TSV = typer.Option(
        help="GC content intervals TSV file",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    PSEUDO_COUNT = typer.Option(
        help="Pseudo count for GC probabilities",
        rich_help_panel=_RICH_HELP_PANEL,
    )


@APP.command()
def from_gfa(
    gfa_file: Annotated[Path, GFAScoresArguments.GFA_FILE],
    gc_scores_tsv: Annotated[
        Path,
        GFAScoresArguments.GC_SCORES_TSV,
    ],
    gc_content_interval_tsv: Annotated[
        Path,
        ScoresOptions.GC_CONTENT_INTERVAL_TSV,
    ] = DEFAULT_GC_CONTENT_INTERVAL_FILE,
    pseudo_count: Annotated[
        int,
        ScoresOptions.PSEUDO_COUNT,
    ] = create.DEFAULT_PSEUDO_COUNT,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Compute GC scores from GFA assembly graph."""
    common_log.init_logger(
        _LOGGER,
        "Computing GC scores from GFA assembly graph.",
        debug,
    )

    # REFACTOR: check file exist in class methods (?)
    if not gfa_file.exists():
        _LOGGER.error("Input GFA file does not exist: %s", gfa_file)
        raise typer.Exit(1)

    if not gc_content_interval_tsv.exists():
        _LOGGER.error(
            "GC content interval file does not exist: %s",
            gc_content_interval_tsv,
        )
        raise typer.Exit(1)

    with io.IntervalStepsReader.open(gc_content_interval_tsv) as interval_steps_reader:
        gc_content_intervals = items.Intervals(
            iter(interval_steps_reader),
        )

    with io.ScoresWriter.open(gc_content_intervals, gc_scores_tsv) as writer:
        for sequence_gc_scores in create.gfa_file_to_gc_scores(
            gfa_file,
            gc_content_intervals,
            pseudo_count=pseudo_count,
        ):
            writer.write_sequence_gc_scores(sequence_gc_scores)
