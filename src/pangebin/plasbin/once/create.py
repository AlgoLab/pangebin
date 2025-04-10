"""Plasbin once creation module."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import gurobipy as gp
import rich.progress as rich_prog

import pangebin.gc_content.items as gc_items
import pangebin.plasbin.bins.items as bins_items
import pangebin.plasbin.config as pb_cfg
import pangebin.plasbin.milp.config as pb_lp_cfg
import pangebin.plasbin.milp.results as pb_lp_res
import pangebin.plasbin.network as net
import pangebin.plasbin.once.milp.input_output as once_lp_io
import pangebin.plasbin.once.milp.models as lp_mod
import pangebin.plasbin.once.milp.views as once_lp_views
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
    gurobi_config: pb_lp_cfg.Gurobi,
    output_directory: Path = once_lp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[
        bins_items.Stats,
        Iterable[bins_items.SequenceNormCoverage],
        once_lp_views.MGCLBStats,
        Path,
    ]
]:
    """Bin contigs of a standardized assembly graph.

    Yields
    ------
    bins_items.Stats
        Stats of the bin
    iterator on bins_items.FragmentNormCoverage
        Contigs and their normalized coverage in the bin
    lp_views.MGCLBStats
        MILP stats
    Path
        Log file path

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
    gurobi_config: pb_lp_cfg.Gurobi,
    output_directory: Path = once_lp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[
        bins_items.Stats,
        Iterable[bins_items.SequenceNormCoverage],
        once_lp_views.MGCLBStats,
        Path,
    ]
]:
    """Bin fragments of a pan-assembly graph.

    Yields
    ------
    bins_items.Stats
        Stats of the bin
    iterator on bins_items.FragmentNormCoverage
        Fragments and their normalized coverage in the bin
    lp_views.MGCLBStats
        MILP stats
    Path
        Log file path

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
        gurobi_config,
        output_directory,
    )


def plasbin(
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    plasbin_config: pb_cfg.Binning,
    gurobi_config: pb_lp_cfg.Gurobi,
    output_directory: Path,
) -> Iterator[
    tuple[
        bins_items.Stats,
        Iterable[bins_items.SequenceNormCoverage],
        once_lp_views.MGCLBStats,
        Path,
    ]
]:
    """Bin DNA sequences of a GFA graph.

    Yields
    ------
    bins_items.Stats
        Stats of the bin
    iterator on bins_items.FragmentNormCoverage
        Fragments and their coverage of the bin
    lp_views.MGCLBStats
        MILP stats
    Path
        Log file path

    Warning
    -------
    The GFA graph will mute.
    """
    gp.setParam(gp.GRB.Param.LogToConsole, 0)
    pb_lp_cfg.configurate_global_gurobi(gurobi_config)

    io_manager = once_lp_io.Manager(output_directory)

    with rich_prog.Progress(console=CONSOLE) as progress:
        remaining_seeds = len(network.seeds())
        binning_task = progress.add_task("Binning", total=remaining_seeds)

        bin_number = 0
        full_feasible_approach = True
        while network.seeds() and full_feasible_approach:
            result = _milp_binning(
                network,
                gc_intervals,
                plasbin_config,
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
def _milp_binning(
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    plasbin_config: pb_cfg.Binning,
    io_manager: once_lp_io.Manager,
    bin_number: int,
) -> tuple[once_lp_views.MGCLBStats, pb_lp_res.Pangebin, Path] | None:
    #
    # MGCLB model
    #
    mgclb_model, mgclb_vars = lp_mod.mgclb(
        network,
        gc_intervals,
        plasbin_config.min_flow(),
        plasbin_config.min_cumulative_len(),
        plasbin_config.circular(),
        plasbin_config.obj_fun_domain(),
    )
    log_file = _run_model(bin_number, mgclb_model, lp_mod.Names.MGCLB, io_manager)

    if mgclb_model.Status != gp.GRB.OPTIMAL:
        return None

    mgclb_stats = once_lp_views.mgclb_stats_from_opt_vars(
        network,
        gc_intervals,
        mgclb_vars,
        plasbin_config.obj_fun_domain(),
    )

    milp_result_values = pb_lp_res.Pangebin.from_optimal_variables(
        network,
        gc_intervals,
        mgclb_vars.flow(),
        mgclb_vars.frag(),
        mgclb_vars.gc(),
    )
    return (mgclb_stats, milp_result_values, log_file)


def _run_model(
    bin_number: int,
    model: gp.Model,
    model_name: lp_mod.Names,
    io_manager: once_lp_io.Manager,
) -> Path:
    """Run model."""
    log_file = io_manager.gurobi_log_path(bin_number, model_name)
    model.Params.LogFile = str(log_file)
    model.optimize()
    return log_file
