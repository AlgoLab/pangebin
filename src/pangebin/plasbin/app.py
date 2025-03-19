"""Plasbin application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gc_content.input_output as gc_io
import pangebin.gfa.input_output as gfa_io
import pangebin.logging as common_log
import pangebin.plasbin.bins.input_output as bin_io
import pangebin.plasbin.config as pb_cfg
import pangebin.plasbin.create as pb_create
import pangebin.plasbin.input_output as pb_io
import pangebin.plasmidness as plm
import pangebin.seed.input_output as seed_io

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
    """PangeBin-flow arguments."""

    PANASSEMBLY_GFA = typer.Argument(
        help="Pan-assembly GFA file",
    )

    GC_SCORES_TSV = typer.Argument(
        help="TSV file with the sequences and their GC scores",
    )

    PLASMIDNESS_TSV = typer.Argument(
        help="TSV file with the sequences and their plasmidness scores",
    )

    SEEDS_TSV = typer.Argument(
        help="TSV file with the seeds",
    )


class MILPOptions:
    """MILP options."""

    _RICH_HELP_PANEL = "MILP options"

    MCF_COEFFICIENT = typer.Option(
        help="MCF coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    MGC_COEFFICIENT = typer.Option(
        help="MGC coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--milp-cfg",
        help="The configuration file path",
        rich_help_panel=_RICH_HELP_PANEL,
    )


class IOOptions:
    """Input-output options."""

    _RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        help="Output directory",
        rich_help_panel=_RICH_HELP_PANEL,
    )


@APP.command()
def plasbin(
    panassembly_gfa: Annotated[Path, Arguments.PANASSEMBLY_GFA],
    gc_scores_tsv: Annotated[Path, Arguments.GC_SCORES_TSV],
    plasmidness_tsv: Annotated[Path, Arguments.PLASMIDNESS_TSV],
    seeds_tsv: Annotated[Path, Arguments.SEEDS_TSV],
    # MILP options
    mcf_coefficient: Annotated[
        float,
        MILPOptions.MCF_COEFFICIENT,
    ] = pb_cfg.Config.DEFAULT_GAMMA_MCF,
    mgc_coefficient: Annotated[
        float,
        MILPOptions.MGC_COEFFICIENT,
    ] = pb_cfg.Config.DEFAULT_GAMMA_MGC,
    milp_cfg_yaml: Annotated[Path | None, MILPOptions.CONFIG_FILE] = None,
    # IO options
    outdir: Annotated[Path, IOOptions.OUTPUT_DIR] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """PlasBin application."""
    common_log.init_logger(_LOGGER, "Running PlasBin.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    io_manager = pb_io.Manager(pb_io.Config(output_directory=outdir))

    milp_config = (
        pb_cfg.Config.from_yaml(milp_cfg_yaml)
        if milp_cfg_yaml is not None
        else pb_cfg.Config(
            gamma_mcf=mcf_coefficient,
            gamma_mgc=mgc_coefficient,
        )
    )

    with gc_io.Reader.open(gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    with plm.Reader.open(plasmidness_tsv) as plasmidness_fin:
        plasmidness = list(plasmidness_fin)
    with seed_io.Reader.open(seeds_tsv) as seeds_fin:
        seeds = list(seeds_fin)

    for k, (bin_stats, seq_normcovs) in enumerate(
        pb_create.plasbin(
            gfa_io.from_file(panassembly_gfa),
            intervals,
            gc_scores,
            plasmidness,
            seeds,
            milp_config,
            io_manager.config().output_directory(),
        ),
    ):
        io_manager.bin_outdir(k).mkdir(parents=True, exist_ok=True)
        bin_stats.to_yaml(io_manager.bin_stats_path(k))
        with bin_io.Writer.open(io_manager.bin_seq_normcov_path(k)) as fout:
            for seq_normcov in seq_normcovs:
                fout.write_sequence_normcov(
                    seq_normcov.identifier(),
                    seq_normcov.normalized_coverage(),
                )
