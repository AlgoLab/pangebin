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
import pangebin.plasbin.decomp.milp.input_output as lp_io
import pangebin.plasbin.input_output as pb_io
import pangebin.plasbin.milp.config as pb_lp_cfg
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net
import pangebin.plasmidness.input_output as plm_io
import pangebin.seed.input_output as seed_io
from pangebin.plasbin.decomp import create

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(
    name="decomp",
    help="Hiearchical decomposition binning method",
    rich_markup_mode="rich",
)


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


@APP.command("asm")
def plasbin_assembly(
    assembly_gfa: Annotated[Path, AsmArguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, AsmArguments.SEED_CONTIGS_TSV],
    contig_gc_scores_tsv: Annotated[Path, AsmArguments.CONTIG_GC_SCORES_TSV],
    contig_plasmidness_tsv: Annotated[Path, AsmArguments.CONTIG_PLASMIDNESS_TSV],
    # Binning options
    sink_arcs_domain: Annotated[
        net.SinkArcsDomain,
        pb_cfg.BinningOptions.SINK_ARCS_DOMAIN,
    ] = pb_cfg.Binning.DEFAULT_SINK_ARCS_DOMAIN,
    min_flow: Annotated[
        float,
        pb_cfg.BinningOptions.MIN_FLOW,
    ] = pb_cfg.Binning.DEFAULT_MIN_FLOW,
    min_cumulative_len: Annotated[
        int,
        pb_cfg.BinningOptions.MIN_CUMULATIVE_LENGTH,
    ] = pb_cfg.Binning.DEFAULT_MIN_CUMULATIVE_LENGTH,
    circular: Annotated[
        bool,
        pb_cfg.BinningOptions.CIRCULAR,
    ] = pb_cfg.Binning.DEFAULT_CIRCULAR,
    obj_fun_domain: Annotated[
        pb_lp_obj.ObjectiveFunctionDomain,
        pb_cfg.BinningOptions.OBJ_FUN_DOMAIN,
    ] = pb_cfg.Binning.DEFAULT_OBJ_FUN_DOMAIN,
    gamma_mcf: Annotated[
        float,
        pb_cfg.BinningOptions.GAMMA_MCF,
    ] = pb_cfg.Binning.DEFAULT_GAMMA_MCF,
    gamma_mgc: Annotated[
        float,
        pb_cfg.BinningOptions.GAMMA_MGC,
    ] = pb_cfg.Binning.DEFAULT_GAMMA_MGC,
    binning_cfg_yaml: Annotated[Path | None, pb_cfg.BinningOptions.CONFIG_FILE] = None,
    # Gurobi options
    mip_gap: Annotated[
        float | None,
        pb_lp_cfg.GurobiOptions.MIP_GAP,
    ] = pb_lp_cfg.Gurobi.DEFAULT_MIP_GAP,
    time_limit: Annotated[
        float | None,
        pb_lp_cfg.GurobiOptions.TIME_LIMIT,
    ] = pb_lp_cfg.Gurobi.DEFAULT_TIME_LIMIT,
    threads: Annotated[
        int | None,
        pb_lp_cfg.GurobiOptions.THREADS,
    ] = pb_lp_cfg.Gurobi.DEFAULT_THREADS,
    gurobi_cfg_yaml: Annotated[Path | None, pb_lp_cfg.GurobiOptions.CONFIG_FILE] = None,
    # IO options
    outdir: Annotated[
        Path,
        pb_io.IOOptions.OUTPUT_DIR,
    ] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """PlasBin an assembly graph."""
    common_log.init_logger(_LOGGER, "Running PlasBin on assembly.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    io_manager = pb_io.Manager(pb_io.Config(output_directory=outdir))

    binning_config = (
        pb_cfg.Binning.from_yaml(binning_cfg_yaml)
        if binning_cfg_yaml is not None
        else pb_cfg.Binning(
            sink_arcs_domain=sink_arcs_domain,
            min_flow=min_flow,
            min_cumulative_len=min_cumulative_len,
            circular=circular,
            obj_fun_domain=obj_fun_domain,
            gamma_mcf=gamma_mcf,
            gamma_mgc=gamma_mgc,
        )
    )

    gurobi_config = (
        pb_lp_cfg.Gurobi.from_yaml(gurobi_cfg_yaml)
        if gurobi_cfg_yaml is not None
        else pb_lp_cfg.Gurobi(mip_gap=mip_gap, time_limit=time_limit, threads=threads)
    )

    with gc_io.ScoresReader.open(contig_gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    with plm_io.Reader.open(contig_plasmidness_tsv) as plasmidness_fin:
        plasmidness = list(plasmidness_fin)
    with seed_io.Reader.open(seed_contigs_tsv) as seeds_fin:
        seeds = list(seeds_fin)

    for k, (bin_stats, seq_normcovs, all_milp_stats, log_files) in enumerate(
        create.plasbin_assembly(
            gfa_io.from_file(assembly_gfa),
            seeds,
            intervals,
            gc_scores,
            plasmidness,
            binning_config,
            gurobi_config,
            io_manager.config().output_directory(),
        ),
    ):
        io_manager.bin_directory(k).mkdir(parents=True, exist_ok=True)
        bin_stats.to_yaml(io_manager.bin_stats_path(k))
        all_milp_stats.to_yaml(io_manager.milp_stats_path(k))
        with bin_io.Writer.open(io_manager.bin_seq_normcov_path(k)) as fout:
            for seq_normcov in seq_normcovs:
                fout.write_sequence_normcov(
                    seq_normcov.identifier(),
                    seq_normcov.normalized_coverage(),
                )
        io_manager.move_gurobi_logs(
            log_files,
            lp_io.Manager.attributes_from_gurobi_log_path,
        )


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
    # Binning options
    sink_arcs_domain: Annotated[
        net.SinkArcsDomain,
        pb_cfg.BinningOptions.SINK_ARCS_DOMAIN,
    ] = pb_cfg.Binning.DEFAULT_SINK_ARCS_DOMAIN,
    min_flow: Annotated[
        float,
        pb_cfg.BinningOptions.MIN_FLOW,
    ] = pb_cfg.Binning.DEFAULT_MIN_FLOW,
    min_cumulative_len: Annotated[
        int,
        pb_cfg.BinningOptions.MIN_CUMULATIVE_LENGTH,
    ] = pb_cfg.Binning.DEFAULT_MIN_CUMULATIVE_LENGTH,
    circular: Annotated[
        bool,
        pb_cfg.BinningOptions.CIRCULAR,
    ] = pb_cfg.Binning.DEFAULT_CIRCULAR,
    obj_fun_domain: Annotated[
        pb_lp_obj.ObjectiveFunctionDomain,
        pb_cfg.BinningOptions.OBJ_FUN_DOMAIN,
    ] = pb_cfg.Binning.DEFAULT_OBJ_FUN_DOMAIN,
    gamma_mcf: Annotated[
        float,
        pb_cfg.BinningOptions.GAMMA_MCF,
    ] = pb_cfg.Binning.DEFAULT_GAMMA_MCF,
    gamma_mgc: Annotated[
        float,
        pb_cfg.BinningOptions.GAMMA_MGC,
    ] = pb_cfg.Binning.DEFAULT_GAMMA_MGC,
    binning_cfg_yaml: Annotated[Path | None, pb_cfg.BinningOptions.CONFIG_FILE] = None,
    # Gurobi options
    mip_gap: Annotated[
        float | None,
        pb_lp_cfg.GurobiOptions.MIP_GAP,
    ] = pb_lp_cfg.Gurobi.DEFAULT_MIP_GAP,
    time_limit: Annotated[
        float | None,
        pb_lp_cfg.GurobiOptions.TIME_LIMIT,
    ] = pb_lp_cfg.Gurobi.DEFAULT_TIME_LIMIT,
    threads: Annotated[
        int | None,
        pb_lp_cfg.GurobiOptions.THREADS,
    ] = pb_lp_cfg.Gurobi.DEFAULT_THREADS,
    gurobi_cfg_yaml: Annotated[Path | None, pb_lp_cfg.GurobiOptions.CONFIG_FILE] = None,
    # IO options
    outdir: Annotated[
        Path,
        pb_io.IOOptions.OUTPUT_DIR,
    ] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """PlasBin a pan-assembly graph."""
    common_log.init_logger(_LOGGER, "Running PlasBin on pan-assembly.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    io_manager = pb_io.Manager(pb_io.Config(output_directory=outdir))

    binning_config = (
        pb_cfg.Binning.from_yaml(binning_cfg_yaml)
        if binning_cfg_yaml is not None
        else pb_cfg.Binning(
            sink_arcs_domain=sink_arcs_domain,
            min_flow=min_flow,
            min_cumulative_len=min_cumulative_len,
            circular=circular,
            obj_fun_domain=obj_fun_domain,
            gamma_mcf=gamma_mcf,
            gamma_mgc=gamma_mgc,
        )
    )
    _LOGGER.debug("Binning config:\n%s", binning_config.to_dict())

    gurobi_config = (
        pb_lp_cfg.Gurobi.from_yaml(gurobi_cfg_yaml)
        if gurobi_cfg_yaml is not None
        else pb_lp_cfg.Gurobi(mip_gap=mip_gap, time_limit=time_limit, threads=threads)
    )
    _LOGGER.debug("Gurobi config:\n%s", gurobi_config.to_dict())

    with gc_io.ScoresReader.open(fragment_gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    with plm_io.Reader.open(fragment_plasmidness_tsv) as plasmidness_fin:
        plasmidness = list(plasmidness_fin)
    with seed_io.Reader.open(seed_fragments_tsv) as seeds_fin:
        seeds = list(seeds_fin)

    for k, (bin_stats, seq_normcovs, all_milp_stats, log_files) in enumerate(
        create.plasbin_panassembly(
            gfa_io.from_file(panassembly_gfa),
            seeds,
            intervals,
            gc_scores,
            plasmidness,
            binning_config,
            gurobi_config,
            io_manager.config().output_directory(),
        ),
    ):
        io_manager.bin_directory(k).mkdir(parents=True, exist_ok=True)
        bin_stats.to_yaml(io_manager.bin_stats_path(k))
        all_milp_stats.to_yaml(io_manager.milp_stats_path(k))
        with bin_io.Writer.open(io_manager.bin_seq_normcov_path(k)) as fout:
            for seq_normcov in seq_normcovs:
                fout.write_sequence_normcov(
                    seq_normcov.identifier(),
                    seq_normcov.normalized_coverage(),
                )
        io_manager.move_gurobi_logs(
            log_files,
            lp_io.Manager.attributes_from_gurobi_log_path,
        )
