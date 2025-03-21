"""Plasbin creation module."""

import logging
from collections.abc import Iterable, Iterator
from pathlib import Path

import gfapy  # type: ignore[import-untyped]
import gurobipy
import rich.progress as rich_prog

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.bins.item as bins_item
import pangebin.plasbin.milp.input_output as milp_io
import pangebin.plasbin.milp.models as milp_models
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

    io_manager = milp_io.Manager(output_directory)

    with rich_prog.Progress(console=CONSOLE) as progress:
        remaining_seeds = len(network.seeds())
        binning_task = progress.add_task("Binning", total=remaining_seeds)

        bin_number = 0
        while network.seeds():
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
            mcf_stats = milp_views.MCFStats(mcf_model.ObjVal)
            #
            # MGC model
            #
            mcf_vars.start_with_previous_values(mcf_vars)
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
            mgc_stats = milp_views.MGCStats(
                milp_vars.coverage_score(network, mcf_vars).getValue(),
                mgc_model.ObjVal,
            )
            #
            # MPS model
            #
            mgc_vars.start_with_previous_values(mgc_vars)
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
            mps_obj_val = mps_model.ObjVal
            #
            # MPS' model
            #
            mps_vars.start_with_previous_values(mps_vars)
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

            mps_stats = milp_views.MPSStats(
                milp_vars.coverage_score(network, mcf_vars).getValue(),
                milp_vars.gc_probability_score(
                    network,
                    gc_intervals,
                    mgc_vars,
                ).getValue(),
                mps_obj_val,
            )
            yield (
                bins_item.Stats(
                    _cumulative_fragment_id_length(mps_vars.mcf_vars(), network),
                    milp_vars.total_flow_value(mps_vars.mcf_vars()),
                    milp_vars.active_gc_content_interval(
                        gc_intervals,
                        mps_vars.mgc_vars(),
                    ),
                    mcf_stats,
                    mgc_stats,
                    mps_stats,
                ),
                _fragment_norm_coverages(mps_vars.mcf_vars(), network),
                log_files.copy(),
            )
            _update_network(network, mps_vars.mcf_vars())
            bin_number += 1

            progress.update(
                binning_task,
                advance=(remaining_seeds - len(network.seeds())),
            )
            remaining_seeds = len(network.seeds())
    _LOGGER.info("Find %u bins.", bin_number)


def _cumulative_fragment_id_length(
    mcf_vars: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> int:
    """Get cumulative fragment length."""
    return sum(
        gfa_segment.length(network.gfa_graph().segment(frag_id))
        for frag_id in milp_vars.active_fragments(network, mcf_vars)
    )


def _fragment_norm_coverages(
    mcf_vars: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> Iterator[bins_item.FragmentNormCoverage]:
    """Get fragment normalized coverages."""
    total_flow = milp_vars.total_flow_value(mcf_vars)
    return (
        bins_item.FragmentNormCoverage(
            frag_id,
            milp_vars.incoming_flow_forward_reverse(
                frag_id,
                network,
                mcf_vars,
            ).getValue()
            / total_flow,
        )
        for frag_id in milp_vars.active_fragments(network, mcf_vars)
    )


def _update_network(
    network: pb_network.Network,
    mcf_vars: milp_vars.MaxCovFlow,
) -> None:
    """Update network."""
    for frag_id in milp_vars.active_fragments(network, mcf_vars):
        network.reduce_coverage(
            frag_id,
            milp_vars.incoming_flow_forward_reverse(
                frag_id,
                network,
                mcf_vars,
            ).getValue(),
        )
