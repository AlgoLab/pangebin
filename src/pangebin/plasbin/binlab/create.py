"""Plasbin binlab creation module."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import gurobipy as gp
import rich.progress as rich_prog

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.binlab.config as binlab_cfg
import pangebin.plasbin.binlab.milp.input_output as lp_io
import pangebin.plasbin.binlab.milp.models as lp_mod
import pangebin.plasbin.binlab.milp.views as lp_views
import pangebin.plasbin.bins.items as bins_items
import pangebin.plasbin.config as pb_cfg
import pangebin.plasbin.milp.config as pb_lp_cfg
import pangebin.plasbin.milp.results as pb_lp_res
import pangebin.plasbin.network as net
from pangebin.pblog import CONSOLE

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    import gfapy  # type: ignore[import-untyped]

_LOGGER = logging.getLogger(__name__)


def plasbin_assembly(  # noqa: PLR0913
    assembly_graph: gfapy.Gfa,
    seed_contigs: Iterable[str],
    gc_intervals: gc_items.Intervals,
    contig_gc_scores: Iterable[gc_items.SequenceGCScores],
    contig_plasmidness: Iterable[tuple[str, float]],
    plasbin_config: pb_cfg.Binning,
    binlab_config: binlab_cfg.Binlab,
    gurobi_config: pb_lp_cfg.Gurobi,
    output_directory: Path = lp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[
        bins_items.Stats,
        Iterable[bins_items.SequenceNormCoverage],
        lp_views.StatsContainer,
        list[Path],
    ]
]:
    """Bin contigs of a standardized assembly graph.

    Yields
    ------
    bins_items.Stats
        Stats of the bin
    iterator on bins_items.FragmentNormCoverage
        Contigs and their normalized coverage in the bin
    lp_views.StatsContainer
        MILP stats
    list of Path
        Paths to the log files

    Warning
    -------
    The GFA assembly graph will mute.
    """
    return plasbin(
        net.Network.from_asm_graph(
            assembly_graph,
            seed_contigs,
            contig_gc_scores,
            contig_plasmidness,
            plasbin_config.sink_arcs_domain(),
        ),
        gc_intervals,
        plasbin_config,
        binlab_config,
        gurobi_config,
        output_directory,
    )


def plasbin_panassembly(  # noqa: PLR0913
    panasm_graph: gfapy.Gfa,
    seed_fragments: Iterable[str],
    gc_intervals: gc_items.Intervals,
    fragment_gc_scores: Iterable[gc_items.SequenceGCScores],
    fragment_plasmidness: Iterable[tuple[str, float]],
    plasbin_config: pb_cfg.Binning,
    binlab_config: binlab_cfg.Binlab,
    gurobi_config: pb_lp_cfg.Gurobi,
    output_directory: Path = lp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[
        bins_items.Stats,
        Iterable[bins_items.SequenceNormCoverage],
        lp_views.StatsContainer,
        list[Path],
    ]
]:
    """Bin fragments of a pan-assembly graph.

    Yields
    ------
    bins_items.Stats
        Stats of the bin
    iterator on bins_items.FragmentNormCoverage
        Fragments and their normalized coverage in the bin
    lp_views.StatsContainer
        MILP stats
    list of Path
        Paths to the log files

    Warning
    -------
    The GFA pan-assembly graph will mute.
    """
    return plasbin(
        net.Network.from_panasm_graph(
            panasm_graph,
            seed_fragments,
            fragment_gc_scores,
            fragment_plasmidness,
            plasbin_config.sink_arcs_domain(),
        ),
        gc_intervals,
        plasbin_config,
        binlab_config,
        gurobi_config,
        output_directory,
    )


def plasbin(  # noqa: PLR0913
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    plasbin_config: pb_cfg.Binning,
    binlab_config: binlab_cfg.Binlab,
    gurobi_config: pb_lp_cfg.Gurobi,
    output_directory: Path,
) -> Iterator[
    tuple[
        bins_items.Stats,
        Iterable[bins_items.SequenceNormCoverage],
        lp_views.StatsContainer,
        list[Path],
    ]
]:
    """Bin DNA sequences of a GFA graph.

    Yields
    ------
    bins_items.Stats
        Stats of the bin
    iterator on bins_items.FragmentNormCoverage
        Fragments and their coverage of the bin
    lp_views.StatsContainer
        MILP stats
    list of Path
        Paths to the log files

    Warning
    -------
    The GFA graph will mute.
    """
    gp.setParam(gp.GRB.Param.LogToConsole, 0)
    pb_lp_cfg.configurate_global_gurobi(gurobi_config)

    io_manager = lp_io.Manager(output_directory)

    with rich_prog.Progress(console=CONSOLE) as progress:
        remaining_seeds = len(network.seeds())
        binning_task = progress.add_task("Binning", total=remaining_seeds)

        bin_number = 0
        full_feasible_approach = True
        while network.seeds() and full_feasible_approach:
            result = _hierarchical_binning(
                network,
                gc_intervals,
                plasbin_config,
                binlab_config,
                io_manager,
                bin_number,
            )
            if result is not None:
                milp_stats, milp_result_values, log_files = result
                fragment_norm_coverages, norm_coverage = (
                    pb_lp_res.fragment_norm_coverages(
                        milp_result_values,
                        plasbin_config.circular(),
                    )
                )
                yield (
                    bins_items.Stats(
                        milp_result_values.cumulative_length(),
                        milp_result_values.total_flow(),
                        norm_coverage,
                        milp_result_values.gc_interval(),
                    ),
                    fragment_norm_coverages,
                    milp_stats,
                    log_files,
                )
                pb_lp_res.update_network(
                    network,
                    milp_result_values,
                    plasbin_config.min_flow(),
                )
                bin_number += 1

                progress.update(
                    binning_task,
                    advance=(remaining_seeds - len(network.seeds())),
                )
                remaining_seeds = len(network.seeds())
            else:
                full_feasible_approach = False
                progress.update(binning_task, advance=remaining_seeds)

    _LOGGER.info("Find %u bins.", bin_number)


# REFACTOR return other than None
def _hierarchical_binning(  # noqa: PLR0913
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    plasbin_config: pb_cfg.Binning,
    binlab_config: binlab_cfg.Binlab,
    io_manager: lp_io.Manager,
    bin_number: int,
) -> tuple[lp_views.StatsContainer, pb_lp_res.Pangebin, list[Path]] | None:
    log_files: list[Path] = []
    #
    # MBS model
    #
    mbs_model, mbs_vars = lp_mod.mbs(
        network,
        plasbin_config.min_flow(),
        plasbin_config.min_cumulative_len(),
        plasbin_config.circular(),
        plasbin_config.obj_fun_domain(),
    )
    log_files.append(_run_model(bin_number, mbs_model, lp_mod.Names.MBS, io_manager))

    if mbs_model.Status != gp.GRB.OPTIMAL:
        return None

    mbs_stats = lp_views.mbs_stats_from_opt_vars(
        network,
        mbs_vars,
        plasbin_config.obj_fun_domain(),
    )
    #
    # MLS model
    #
    mls_model, mls_vars = lp_mod.mls_from_mbs(
        mbs_model,
        mbs_vars,
        plasbin_config.obj_fun_domain(),
        network,
        gc_intervals,
        binlab_config.gamma_mbs(),
    )
    log_files.append(_run_model(bin_number, mls_model, lp_mod.Names.MLS, io_manager))

    if mls_model.Status != gp.GRB.OPTIMAL:
        return None

    mls_stats = lp_views.mls_stats_from_opt_vars(
        network,
        gc_intervals,
        mls_vars,
        plasbin_config.obj_fun_domain(),
    )
    #
    # MRBS model
    #
    mrbs_model, mrbs_vars = lp_mod.mrbs_from_mls(
        mls_model,
        mls_vars,
        network,
        gc_intervals,
        plasbin_config.obj_fun_domain(),
    )
    log_files.append(_run_model(bin_number, mrbs_model, lp_mod.Names.MRBS, io_manager))

    if mrbs_model.Status != gp.GRB.OPTIMAL:
        return None

    mrbs_stats = lp_views.mrbs_stats_from_opt_vars(
        network,
        gc_intervals,
        mrbs_vars,
        plasbin_config.obj_fun_domain(),
    )

    milp_result_values = pb_lp_res.Pangebin.from_optimal_variables(
        network,
        gc_intervals,
        mls_vars.flow(),
        mls_vars.frag(),
        mls_vars.gc(),
    )
    return (
        lp_views.StatsContainer(mbs_stats, mls_stats, mrbs_stats),
        milp_result_values,
        log_files.copy(),
    )


def _run_model(
    bin_number: int,
    model: gp.Model,
    model_name: lp_mod.Names,
    io_manager: lp_io.Manager,
) -> Path:
    """Run model."""
    log_file = io_manager.gurobi_log_path(bin_number, model_name)
    model.Params.LogFile = str(log_file)
    model.optimize()
    return log_file
