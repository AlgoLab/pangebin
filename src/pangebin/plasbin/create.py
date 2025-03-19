"""Plasbin creation module."""

from collections.abc import Iterable, Iterator
from pathlib import Path

import gfapy  # type: ignore[import-untyped]

import pangebin.gc_content.items as gc_items
import pangebin.gfa.segment as gfa_segment
import pangebin.plasbin.bins.item as bins_item
import pangebin.plasbin.input_output as pb_io
import pangebin.plasbin.milp.models as milp_models
import pangebin.plasbin.milp.variables as milp_vars
import pangebin.plasbin.milp.views as milp_views
import pangebin.plasbin.network as pb_network
from pangebin.plasbin.config import Config


def plasbin(  # noqa: PLR0913
    panasm_graph: gfapy.Gfa,
    gc_intervals: gc_items.Intervals,
    gc_scores: Iterable[gc_items.SequenceGCScores],
    plasmidness: Iterable[tuple[str, float]],
    seeds: Iterable[str],
    config: Config,
    log_directory: Path = pb_io.Config.DEFAULT_OUTPUT_DIR,
) -> Iterator[tuple[bins_item.Stats, Iterable[bins_item.FragmentNormCoverage]]]:
    """Bin fragments from a pan-assembly graph.

    Warning
    -------
    The pan-assembly graph can mute.
    """
    network = pb_network.Network(panasm_graph, seeds, gc_scores, plasmidness)
    bin_number = 0
    while network.seeds():
        #
        # MCF model
        #
        mcf_model, mcf_vars = milp_models.mcf(network)
        mcf_model.Params.LogToConsole = 0
        mcf_model.Params.LogFile = str(
            pb_io.gurobi_log_path(log_directory, bin_number, milp_models.Names.MCF),
        )
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
        mgc_model.Params.LogFile = str(
            pb_io.gurobi_log_path(log_directory, bin_number, milp_models.Names.MGC),
        )
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
        mps_model.Params.LogFile = str(
            pb_io.gurobi_log_path(log_directory, bin_number, milp_models.Names.MPS),
        )
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
        mps_prime_model.Params.LogFile = str(
            pb_io.gurobi_log_path(
                log_directory,
                bin_number,
                milp_models.Names.MPS_PRIME,
            ),
        )
        mps_prime_model.optimize()

        mps_stats = milp_views.MPSStats(
            milp_vars.coverage_score(network, mcf_vars).getValue(),
            milp_vars.gc_probability_score(network, gc_intervals, mgc_vars).getValue(),
            mps_obj_val,
        )
        yield (
            bins_item.Stats(
                cumulative_fragment_id_length(mps_vars.mcf_vars(), network),
                milp_vars.total_flow_value(mps_vars.mcf_vars()),
                milp_vars.active_gc_content_interval(gc_intervals, mps_vars.mgc_vars()),
                mcf_stats,
                mgc_stats,
                mps_stats,
            ),
            fragment_norm_coverages(mps_vars.mcf_vars(), network),
        )
        update_network(network, mps_vars.mcf_vars())
        bin_number += 1


def cumulative_fragment_id_length(
    mcf_vars: milp_vars.MaxCovFlow,
    network: pb_network.Network,
) -> int:
    """Get cumulative fragment length."""
    return sum(
        gfa_segment.length(network.panasm_graph().segment(frag_id))
        for frag_id in milp_vars.active_fragments(network, mcf_vars)
    )


def fragment_norm_coverages(
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


def update_network(network: pb_network.Network, mcf_vars: milp_vars.MaxCovFlow) -> None:
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
