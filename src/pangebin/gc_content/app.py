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

    SEQUENCES_PROBA_SCORES_TSV = typer.Argument(
        help="Output TSV file with the sequences and their probability scores",
    )


class ScoresOptions:
    """GC content options."""

    _RICH_HELP_PANEL = "GFA GC content options"

    GC_CONTENT_INTERVAL_FILE = typer.Option(
        help="GC content intervals file",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    PSEUDO_COUNT = typer.Option(
        help="Pseudo count for GC probabilities",
        rich_help_panel=_RICH_HELP_PANEL,
    )


@APP.command()
def from_gfa(
    gfa_file: Annotated[Path, GFAScoresArguments.GFA_FILE],
    sequences_proba_scores_tsv: Annotated[
        Path,
        GFAScoresArguments.SEQUENCES_PROBA_SCORES_TSV,
    ],
    gc_content_interval_file: Annotated[
        Path | None,
        ScoresOptions.GC_CONTENT_INTERVAL_FILE,
    ] = None,
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

    if gc_content_interval_file is None:
        gc_content_interval_file = DEFAULT_GC_CONTENT_INTERVAL_FILE

    if not gc_content_interval_file.exists():
        _LOGGER.error(
            "GC content interval file does not exist: %s",
            gc_content_interval_file,
        )
        raise typer.Exit(1)

    gc_content_intervals = items.Intervals.from_file(gc_content_interval_file)

    with io.Writer.open(gc_content_intervals, sequences_proba_scores_tsv) as writer:
        for sequence_proba_scores in create.gfa_file_to_gc_scores(
            gfa_file,
            gc_content_intervals,
            pseudo_count=pseudo_count,
        ):
            writer.write_sequence_proba_scores(sequence_proba_scores)
