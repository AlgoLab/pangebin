"""Plasbin creation module."""

import logging
from collections.abc import Iterable, Iterator
from pathlib import Path

import gfapy  # type: ignore[import-untyped]
import gurobipy
import rich.progress as rich_prog

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.bins.items as bins_item
import pangebin.plasbin.milp.input_output as milp_io
import pangebin.plasbin.milp.models as milp_models
import pangebin.plasbin.milp.results as milp_results
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.milp.views as milp_views
import pangebin.plasbin.network as pb_network
from pangebin.logging import CONSOLE
from pangebin.plasbin.config import Config

_LOGGER = logging.getLogger(__name__)


def plasbin_assembly(  # noqa: PLR0913
    assembly_graph: gfapy.Gfa,
    seed_contigs: Iterable[str],
    gc_intervals: gc_items.Intervals,
    contig_gc_scores: Iterable[gc_items.SequenceGCScores],
    contig_plasmidness: Iterable[tuple[str, float]],
    config: Config,
    output_directory: Path = milp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[bins_item.Stats, Iterable[bins_item.FragmentNormCoverage], list[Path]]
]:
    """Bin contigs of a standardized assembly graph.

    Yields
    ------
    bins_item.Stats
        Stats of the bin
    iterator on bins_item.FragmentNormCoverage
        Contigs and their normalized coverage in the bin
    list of Path
        Paths to the log files

    Warning
    -------
    The GFA assembly graph will mute.
    """
    return plasbin(
        pb_network.Network.from_asm_graph(
            assembly_graph,
            seed_contigs,
            contig_gc_scores,
            contig_plasmidness,
        ),
        gc_intervals,
        config,
        output_directory,
    )


def plasbin_panassembly(  # noqa: PLR0913
    panasm_graph: gfapy.Gfa,
    seed_fragments: Iterable[str],
    gc_intervals: gc_items.Intervals,
    fragment_gc_scores: Iterable[gc_items.SequenceGCScores],
    fragment_plasmidness: Iterable[tuple[str, float]],
    config: Config,
    output_directory: Path = milp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[bins_item.Stats, Iterable[bins_item.FragmentNormCoverage], list[Path]]
]:
    """Bin fragments of a pan-assembly graph.

    Yields
    ------
    bins_item.Stats
        Stats of the bin
    iterator on bins_item.FragmentNormCoverage
        Fragments and their normalized coverage in the bin
    list of Path
        Paths to the log files

    Warning
    -------
    The GFA pan-assembly graph will mute.
    """
    return plasbin(
        pb_network.Network.from_panasm_graph(
            panasm_graph,
            seed_fragments,
            fragment_gc_scores,
            fragment_plasmidness,
        ),
        gc_intervals,
        config,
        output_directory,
    )


def plasbin(
    network: pb_network.Network,
    gc_intervals: gc_items.Intervals,
    config: Config,
    output_directory: Path,
) -> Iterator[
    tuple[bins_item.Stats, Iterable[bins_item.FragmentNormCoverage], list[Path]]
]:
    """Bin DNA sequences of a GFA graph.

    Yields
    ------
    bins_item.Stats
        Stats of the bin
    iterator on bins_item.FragmentNormCoverage
        Fragments and their coverage of the bin
    list of Path
        Paths to the log files

    Warning
    -------
    The GFA graph will mute.
    """
    gurobipy.setParam(gurobipy.GRB.Param.LogToConsole, 0)
    # HACK Gurobi: params
    # REFACTOR use Gurobi Config and rename other config
    gurobipy.setParam(gurobipy.GRB.Param.MIPGap, 0.0)

    io_manager = milp_io.Manager(output_directory)

    with rich_prog.Progress(console=CONSOLE) as progress:
        remaining_seeds = len(network.seeds())
        binning_task = progress.add_task("Binning", total=remaining_seeds)

        bin_number = 0
        found_positive_flow = True
        while network.seeds() and found_positive_flow:
            result = _hierarchical_binning(
                network,
                gc_intervals,
                config,
                io_manager,
                bin_number,
            )
            if result is not None:
                bin_stats, milp_result_values, log_files = result
                yield (
                    bin_stats,
                    _fragment_norm_coverages(milp_result_values),
                    log_files,
                )
                _update_network(network, milp_result_values)
                bin_number += 1

                progress.update(
                    binning_task,
                    advance=(remaining_seeds - len(network.seeds())),
                )
                remaining_seeds = len(network.seeds())
            else:
                # REFACTOR rename this boolean while condition
                found_positive_flow = False
                progress.update(binning_task, advance=remaining_seeds)

    _LOGGER.info("Find %u bins.", bin_number)


def plasbin_multiobj(  # noqa: PLR0913
    asm_graph: gfapy.Gfa,
    seed_fragments: Iterable[str],
    gc_intervals: gc_items.Intervals,
    fragment_gc_scores: Iterable[gc_items.SequenceGCScores],
    fragment_plasmidness: Iterable[tuple[str, float]],
    config: Config,
    output_directory: Path = milp_io.Manager.DEFAULT_OUTPUT_DIR,
) -> Iterator[
    tuple[bins_item.Stats, Iterable[bins_item.FragmentNormCoverage], list[Path]]
]:
    """Bin fragments of a GFA assembly graph using multi objective."""
    network = pb_network.Network.from_asm_graph(
        asm_graph,
        seed_fragments,
        fragment_gc_scores,
        fragment_plasmidness,
    )
    log_files: list[Path] = []

    gurobipy.setParam(gurobipy.GRB.Param.LogToConsole, 0)
    # TODO use config MILP
    gurobipy.setParam(gurobipy.GRB.Param.MIPGap, 0.0)

    io_manager = milp_io.Manager(output_directory)

    with rich_prog.Progress(console=CONSOLE) as progress:
        remaining_seeds = len(network.seeds())
        binning_task = progress.add_task("Binning", total=remaining_seeds)

        bin_number = 0
        is_feasible = True
        while network.seeds() and is_feasible:
            result = None

            m, var = milp_models.multiobjective(network, gc_intervals)
            log_files.append(
                io_manager.gurobi_log_path(bin_number, milp_models.Names.MCF),
            )
            m.Params.LogFile = str(log_files[-1])
            m.optimize()
            if m.Status == gurobipy.GRB.OPTIMAL:
                milp_result_values = milp_results.Pangebin.from_optimal_variables(
                    network,
                    gc_intervals,
                    var.mcf_vars(),
                    var.mgc_vars(),
                )
                result = (
                    bins_item.Stats(
                        _cumulative_fragment_id_length(
                            milp_result_values.active_fragments(),
                            network,
                        ),
                        milp_result_values.total_flow(),
                        milp_result_values.gc_interval(),
                        milp_views.MCFStats(0),
                        milp_views.MGCStats(0, 0),
                        milp_views.MPSStats(0, 0, 0),
                    ),
                    milp_result_values,
                    log_files.copy(),
                )

            if result is not None:
                bin_stats, milp_result_values, log_files = result
                yield (
                    bin_stats,
                    _fragment_norm_coverages(milp_result_values),
                    log_files,
                )
                _update_network(network, milp_result_values)
                bin_number += 1

                progress.update(
                    binning_task,
                    advance=(remaining_seeds - len(network.seeds())),
                )
                remaining_seeds = len(network.seeds())
            else:
                # REFACTOR rename this boolean while condition
                is_feasible = False
                progress.update(binning_task, advance=remaining_seeds)

    _LOGGER.info("Find %u bins.", bin_number)


# REFACTOR return other than None
def _hierarchical_binning(
    network: pb_network.Network,
    gc_intervals: gc_items.Intervals,
    config: Config,
    io_manager: milp_io.Manager,
    bin_number: int,
) -> tuple[bins_item.Stats, milp_results.Pangebin, list[Path]] | None:
    log_files: list[Path] = []
    #
    # MCF model
    #
    mcf_model, mcf_vars = milp_models.mcf(network)
    log_files.append(
        io_manager.gurobi_log_path(bin_number, milp_models.Names.MCF),
    )
    mcf_model.Params.LogFile = str(log_files[-1])
    mcf_model.optimize()

    if mcf_model.Status != gurobipy.GRB.OPTIMAL:
        return None

    mcf_stats = milp_views.MCFStats(
        milp_vars.coverage_score(network, mcf_vars).getValue(),
    )
    #
    # MGC model
    #
    # BUG User mip start violoates constraint while should not
    # mcf_vars.start_with_previous_values(mcf_vars)
    mgc_model, mgc_vars = milp_models.mgc_from_mcf(
        mcf_model,
        mcf_vars,
        network,
        gc_intervals,
        config.gamma_mgc(),
    )
    log_files.append(
        io_manager.gurobi_log_path(bin_number, milp_models.Names.MGC),
    )
    mgc_model.Params.LogFile = str(log_files[-1])
    mgc_model.optimize()

    if mgc_model.Status != gurobipy.GRB.OPTIMAL:
        return None

    mgc_stats = milp_views.MGCStats(
        milp_vars.coverage_score(network, mgc_vars.mcf_vars()).getValue(),
        milp_vars.gc_probability_score(network, gc_intervals, mgc_vars).getValue(),
    )
    #
    # MPS model
    #
    # BUG User mip start violoates constraint while should not
    # mgc_vars.start_with_previous_values(mgc_vars)
    mps_model, mps_vars = milp_models.mps_from_mgc(
        mgc_model,
        mgc_vars,
        network,
        gc_intervals,
        config.gamma_mgc(),
    )
    log_files.append(
        io_manager.gurobi_log_path(bin_number, milp_models.Names.MPS),
    )
    mps_model.Params.LogFile = str(log_files[-1])
    mps_model.optimize()

    if mps_model.Status != gurobipy.GRB.OPTIMAL:
        return None

    mps_obj_val = milp_vars.plasmidness_score(
        network,
        gc_intervals,
        mps_vars,
    ).getValue()
    #
    # MPS' model
    #
    # BUG User mip start violoates constraint while should not
    # mps_vars.start_with_previous_values(mps_vars)
    mps_prime_model, _ = milp_models.mps_prime_from_mps(
        mps_model,
        mps_vars,
        network,
        gc_intervals,
    )
    log_files.append(
        io_manager.gurobi_log_path(bin_number, milp_models.Names.MPS_PRIME),
    )
    mps_prime_model.Params.LogFile = str(log_files[-1])
    mps_prime_model.optimize()

    if mps_prime_model.Status != gurobipy.GRB.OPTIMAL:
        return None

    mps_stats = milp_views.MPSStats(
        milp_vars.coverage_score(network, mps_vars.mcf_vars()).getValue(),
        milp_vars.gc_probability_score(
            network,
            gc_intervals,
            mps_vars.mgc_vars(),
        ).getValue(),
        mps_obj_val,
    )

    milp_result_values = milp_results.Pangebin.from_optimal_variables(
        network,
        gc_intervals,
        mps_vars.mcf_vars(),
        mps_vars.mgc_vars(),
    )
    return (
        bins_item.Stats(
            _cumulative_fragment_id_length(
                milp_result_values.active_fragments(),
                network,
            ),
            milp_result_values.total_flow(),
            milp_result_values.gc_interval(),
            mcf_stats,
            mgc_stats,
            mps_stats,
        ),
        milp_result_values,
        log_files.copy(),
    )


def _cumulative_fragment_id_length(
    active_fragments: Iterable[str],
    network: pb_network.Network,
) -> int:
    """Get cumulative fragment length."""
    return sum(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in active_fragments
    )


def _fragment_norm_coverages(
    milp_result_values: milp_results.Pangebin,
) -> Iterable[bins_item.FragmentNormCoverage]:
    """Get fragment normalized coverages."""
    return (
        bins_item.FragmentNormCoverage(
            frag_id,
            incoming_flow / milp_result_values.total_flow(),
        )
        for frag_id, incoming_flow in milp_result_values.fragments_incoming_flow()
    )


def _update_network(
    network: pb_network.Network,
    milp_result_values: milp_results.Pangebin,
) -> None:
    """Update network."""
    for frag_id, incoming_flow in milp_result_values.fragments_incoming_flow():
        network.reduce_coverage(frag_id, incoming_flow)
