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
import pangebin.plasmidness.input_output as plm_io
import pangebin.seed.input_output as seed_io

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class AsmArguments:
    """PangeBin-flow on assembly arguments."""

    ASSEMBLY_GFA = typer.Argument(
        help="Assembly GFA file",
    )

    SEED_CONTIGS_TSV = typer.Argument(
        help="TSV file with the seed contigs",
    )

    CONTIG_GC_SCORES_TSV = typer.Argument(
        help="TSV file with the contigs and their GC scores",
    )

    CONTIG_PLASMIDNESS_TSV = typer.Argument(
        help="TSV file with the contigs and their plasmidness scores",
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


@APP.command("asm")
def plasbin_assembly(
    assembly_gfa: Annotated[Path, AsmArguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, AsmArguments.SEED_CONTIGS_TSV],
    contig_gc_scores_tsv: Annotated[Path, AsmArguments.CONTIG_GC_SCORES_TSV],
    contig_plasmidness_tsv: Annotated[Path, AsmArguments.CONTIG_PLASMIDNESS_TSV],
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
    common_log.init_logger(_LOGGER, "Running PlasBin on assembly.", debug)
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

    with gc_io.ScoresReader.open(contig_gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    with plm_io.Reader.open(contig_plasmidness_tsv) as plasmidness_fin:
        plasmidness = list(plasmidness_fin)
    with seed_io.Reader.open(seed_contigs_tsv) as seeds_fin:
        seeds = list(seeds_fin)

    for k, (bin_stats, seq_normcovs, log_files) in enumerate(
        pb_create.plasbin_assembly(
            gfa_io.from_file(assembly_gfa),
            seeds,
            intervals,
            gc_scores,
            plasmidness,
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
        io_manager.move_gurobi_logs(log_files)


class PanasmArguments:
    """PangeBin-flow on pan-assembly arguments."""

    PANASSEMBLY_GFA = typer.Argument(
        help="Pan-assembly GFA file",
    )

    SEED_FRAGMENTS_TSV = typer.Argument(
        help="TSV file with the seed fragments",
    )

    FRAGMENT_GC_SCORES_TSV = typer.Argument(
        help="TSV file with the fragments and their GC scores",
    )

    FRAGMENT_PLASMIDNESS_TSV = typer.Argument(
        help="TSV file with the fragments and their plasmidness scores",
    )


@APP.command("panasm")
def plasbin_panassembly(
    panassembly_gfa: Annotated[Path, PanasmArguments.PANASSEMBLY_GFA],
    seed_fragments_tsv: Annotated[Path, PanasmArguments.SEED_FRAGMENTS_TSV],
    fragment_gc_scores_tsv: Annotated[Path, PanasmArguments.FRAGMENT_GC_SCORES_TSV],
    fragment_plasmidness_tsv: Annotated[Path, PanasmArguments.FRAGMENT_PLASMIDNESS_TSV],
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
    common_log.init_logger(_LOGGER, "Running PlasBin on pan-assembly.", debug)
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

    with gc_io.ScoresReader.open(fragment_gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    with plm_io.Reader.open(fragment_plasmidness_tsv) as plasmidness_fin:
        plasmidness = list(plasmidness_fin)
    with seed_io.Reader.open(seed_fragments_tsv) as seeds_fin:
        seeds = list(seeds_fin)

    for k, (bin_stats, seq_normcovs, log_files) in enumerate(
        pb_create.plasbin_panassembly(
            gfa_io.from_file(panassembly_gfa),
            seeds,
            intervals,
            gc_scores,
            plasmidness,
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
        io_manager.move_gurobi_logs(log_files)
