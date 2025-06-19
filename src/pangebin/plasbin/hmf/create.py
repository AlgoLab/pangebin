"""HMF creation module."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import gurobipy as gp
import rich.progress as rich_prog

import pangebin.gc_content.items as gc_items
import pangebin.gfa.ops as gfa_ops
import pangebin.plasbin.milp.config as cmn_lp_cfg
import pangebin.plasbin.milp.results as cmn_lp_res
import pangebin.plasbin.network as net
from pangebin.pblog import CONSOLE

from . import bins, managers, results
from . import config as cfg
from . import file_system as fs
from .milp import config as lp_cfg
from .milp import constraints as lp_cst
from .milp import views as lp_views

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    import gfapy  # type: ignore[import-untyped]

_LOGGER = logging.getLogger(__name__)


def plasbin_assembly(  # noqa: PLR0913
    assembly_graph: gfapy.Gfa,
    seed_contigs: Iterable[str],
    contig_gc_scores: Iterable[gc_items.SequenceGCScores],
    contig_plasmidness: Iterable[tuple[str, float]],
    config: cfg.Config,
    gurobi_config: cmn_lp_cfg.Gurobi,
    output_directory: Path,
) -> results.RootReader:
    """Bin contigs of a standardized assembly graph.

    Warning
    -------
    The GFA assembly graph will mute.

    Returns
    -------
    results.RootReader
        Reader for the results
    """
    cmn_lp_cfg.configurate_global_gurobi(gurobi_config)

    best_instances = results.Root.new()

    hmf_root_file_system = fs.Root(output_directory)
    hmf_root_file_system.dir().mkdir(parents=True)

    for ccomp_idx, (subgraph, seeds, gc_scores, plasmidness) in enumerate(
        components(
            assembly_graph,
            seed_contigs,
            contig_gc_scores,
            contig_plasmidness,
        ),
    ):
        _LOGGER.debug("Processing connected component %d", ccomp_idx)
        best_ccomp = bin_ccomp(
            net.Network.from_asm_graph(
                subgraph,
                seeds,
                gc_scores,
                plasmidness,
                config.network().sink_arcs_domain(),
            ),
            config.bin_properties(),
            config.model(),
            hmf_root_file_system.ccomp_file_system(ccomp_idx),
        )
        best_instances.add_connected_component(best_ccomp)

    best_instances.to_yaml(hmf_root_file_system.best_instances_yaml())

    return results.RootReader(best_instances, hmf_root_file_system)


def components(
    graph: gfapy.Gfa,
    seeds: Iterable[str],
    gc_scores: Iterable[gc_items.SequenceGCScores],
    all_plasmidness: Iterable[tuple[str, float]],
) -> list[
    tuple[
        gfapy.Gfa,
        list[str],  # seeds
        list[gc_items.SequenceGCScores],  # gc_scores
        list[tuple[str, float]],  # plasmidness
    ]
]:
    """Get GFA connected components (and attributes).

    Returns
    -------
    list[
        tuple[
            Gfa,
            list[str],
            list[SequenceGCScores],
            list[tuple[str, float]],
        ]
    ]
        Sub gfa, and associated seeds, gc_scores, and plasmidness
    """
    # OPTIMIZE can do better
    sub_graphs = gfa_ops.connected_components(graph)
    seg_cc: dict[str, int] = {
        seg_name: k
        for k, sub_graph in enumerate(sub_graphs)
        for seg_name in sub_graph.segment_names
    }

    sub_seeds: list[list[str]] = [[] for _ in range(len(sub_graphs))]
    for seg_name in seeds:
        sub_seeds[seg_cc[seg_name]].append(seg_name)

    sub_gc_scores: list[list[gc_items.SequenceGCScores]] = [
        [] for _ in range(len(sub_graphs))
    ]
    for gc_score in gc_scores:
        sub_gc_scores[seg_cc[gc_score.sequence_id()]].append(gc_score)

    sub_plasmidness: list[list[tuple[str, float]]] = [
        [] for _ in range(len(sub_graphs))
    ]
    for plasmidness in all_plasmidness:
        sub_plasmidness[seg_cc[plasmidness[0]]].append(plasmidness)

    return [
        (sub_graphs[k], sub_seeds[k], sub_gc_scores[k], sub_plasmidness[k])
        for k in range(len(sub_graphs))
    ]


def bin_ccomp(
    network: net.Network,
    bin_properties: cfg.BinProperties,
    model_config: lp_cfg.Config,
    ccomp_fs_manager: fs.ConnectedComponent,
) -> results.ConnectedComponent:
    """Bin DNA sequences of a GFA graph.

    Warning
    -------
    The network and the associated GFA will mute.
    """
    ccomp_fs_manager.dir().mkdir()

    # Test circular bins
    # * With seeds
    # * Without seeds

    best_instances = results.ConnectedComponent.new()

    for bin_class_manager in managers.iter_bin_class_manager(
        bin_properties,
        network,
        model_config,
    ):
        best_bin_class_instance = search_best_bin_class_instance(
            network,
            bin_class_manager,
            ccomp_fs_manager,
        )
        if best_bin_class_instance is not None:
            best_instances.best_for_topology(
                bin_class_manager.topology(),
            ).set_best_for_seed_constraint(
                bin_class_manager.seed_constraint(),
                best_bin_class_instance,
            )
            _remove_solution_from_network(
                network,
                results.BinClassReader(
                    best_bin_class_instance,
                    ccomp_fs_manager.bin_class_file_system(
                        bin_class_manager.topology(),
                        bin_class_manager.seed_constraint(),
                        best_bin_class_instance.number_of_bins(),
                    ),
                ),
            )

    return best_instances


def search_best_bin_class_instance(
    network: net.Network,
    bin_class_manager: managers.BinClass,
    ccomp_fs_manager: fs.ConnectedComponent,
) -> results.FeasibleInstance | None:
    """Search best bin class instance."""
    best_feasible_instance: results.FeasibleInstance | None = None

    with rich_prog.Progress(console=CONSOLE) as circ_progress:
        circ_task = circ_progress.add_task(
            _progress_message(bin_class_manager),
            total=bin_class_manager.stats().max_number_of_bins(),
        )

        test_new_bin = bool(bin_class_manager.stats().max_number_of_bins())
        best_multi_flow_objective: float | None = None
        while test_new_bin:
            bin_class_manager.new_bin()

            bin_class_fs_manager = ccomp_fs_manager.bin_class_file_system(
                bin_class_manager.topology(),
                bin_class_manager.seed_constraint(),
                bin_class_manager.stats().number_of_active_bins(),
            )

            if bin_class_manager.stats().number_of_active_bins() > 1:
                lp_cst.plasmidness_score_order(
                    bin_class_manager.model().gurobi_model(),
                    bin_class_manager.model().bins_variables(),
                    network,
                    bin_class_manager.model().config().obj_fun_domain(),
                    (
                        bin_class_manager.stats().number_of_active_bins() - 2,
                        bin_class_manager.stats().number_of_active_bins(),
                    ),
                )
                # FIXME embete le gap

            feasible_instance = _solve_instance(
                network,
                bin_class_manager,
                bin_class_fs_manager,
            )
            # REFACTOR avoid if-else here, flatten with function
            # Infeasible because: no circular or not better
            if feasible_instance is None or (
                best_multi_flow_objective is not None
                # BUG TMP ROUND
                and round(bin_class_manager.model().gurobi_model().ObjVal, 4)
                <= round(best_multi_flow_objective, 4)
            ):
                test_new_bin = False
            else:
                best_feasible_instance = feasible_instance
                # FIXME avoid same obj value
                # Because of gap, it can be better while it is the same splitted bin
                best_multi_flow_objective = (
                    bin_class_manager.model().gurobi_model().ObjVal
                )

                bin_class_manager.set_objective_lb(best_multi_flow_objective)

                circ_progress.update(
                    circ_task,
                    advance=1,
                )

                test_new_bin = (
                    bin_class_manager.stats().number_of_active_bins()
                    < bin_class_manager.stats().max_number_of_bins()
                )

        circ_progress.update(
            circ_task,
            advance=(
                bin_class_manager.stats().max_number_of_bins()
                - bin_class_manager.stats().number_of_active_bins()
            ),
        )
    return best_feasible_instance


def _solve_instance(
    network: net.Network,
    bin_class_manager: managers.BinClass,
    bin_class_file_system: fs.BinClass,
) -> results.FeasibleInstance | None:
    bin_class_file_system.dir().mkdir()
    _solve_mfb_model(bin_class_manager.model().gurobi_model(), bin_class_file_system)

    if bin_class_manager.model().gurobi_model().Status != gp.GRB.OPTIMAL:
        return None

    for k in range(bin_class_manager.stats().number_of_active_bins()):
        bin_file_system = bin_class_file_system.bin_file_system(k)
        bin_file_system.dir().mkdir()

        milp_result = cmn_lp_res.Pangebin.from_optimal_vars_without_gc_intervals(
            network,
            bin_class_manager.model().bins_variables()[k].flows(),
            bin_class_manager.model().bins_variables()[k].sub_frag(),
        )
        milp_stats = lp_views.stats_from_opt_bin_vars(
            network,
            bin_class_manager.model().bins_variables()[k],
            bin_class_manager.model().config().obj_fun_domain(),
        )

        fragment_norm_coverages, norm_coverage = cmn_lp_res.fragment_norm_coverages(
            milp_result,
        )

        bin_stats = bins.Stats(
            bin_class_manager.topology(),
            bin_class_manager.seed_constraint(),
            milp_result.cumulative_length(),
            milp_result.total_flow(),
            norm_coverage,
        )

        results.BinWriter(bin_file_system).write(
            bin_stats,
            milp_stats,
            fragment_norm_coverages,
        )

    return results.FeasibleInstance(bin_class_manager.stats().number_of_active_bins())


def _solve_mfb_model(
    model: gp.Model,
    bin_class_fs: fs.BinClass,
) -> Path:
    """Run model."""
    log_file = bin_class_fs.gurobi_log()
    model.Params.LogFile = str(log_file)
    model.optimize()
    return log_file


def _progress_message(bin_class_manager: managers.BinClass) -> str:
    """Return progress message."""
    match bin_class_manager.topology():
        case bins.Topology.CIRCULAR:
            _firs_part = "Binning circular bins"
        case bins.Topology.PARTIALLY_CIRCULAR:
            _firs_part = "Binning partially circular bins"

    match bin_class_manager.seed_constraint():
        case bins.SeedConstraint.REQUIRED:
            _second_part = "with seeds"
        case bins.SeedConstraint.NOT_REQUIRED:
            _second_part = "free of seeds"

    return f"{_firs_part} {_second_part}"


def _remove_solution_from_network(
    network: net.Network,
    best_bin_class_reader: results.BinClassReader,
) -> None:
    """Remove solution instance from network."""
    for bin_result_reader in best_bin_class_reader.bin_readers():
        normalization_factor = bin_result_reader.bin_stats().normalizing_coverage()
        for seq_normcov in bin_result_reader.iter_seq_normcov():
            network.reduce_coverage(
                seq_normcov.identifier(),
                seq_normcov.normalized_coverage() * normalization_factor,
            )
