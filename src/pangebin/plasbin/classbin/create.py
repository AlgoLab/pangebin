"""Plasbin classbin creation module."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import gurobipy as gp
import rich.progress as rich_prog

import pangebin.gc_content.items as gc_items
import pangebin.gfa.ops as gfa_ops
import pangebin.plasbin.bins.items as bins_items
import pangebin.plasbin.classbin.file_system as classbin_lp_io
import pangebin.plasbin.config as cmn_cfg
import pangebin.plasbin.milp.config as cmn_lp_cfg
import pangebin.plasbin.milp.results as cmn_lp_res
import pangebin.plasbin.network as net
from pangebin.pblog import CONSOLE

from .multi_flow.milp import constraints as mfb_cst
from .multi_flow.milp import file_system as mfb_fs
from .multi_flow.milp import models as mfb_model
from .multi_flow.milp import objectives as mfb_obj
from .multi_flow.milp import variables as mfb_var
from .multi_flow.milp import views as mfb_views

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
    plasbin_config: cmn_cfg.Binning,
    gurobi_config: cmn_lp_cfg.Gurobi,
    output_directory: Path,
) -> Iterator[
    tuple[
        list[
            tuple[
                bins_items.Stats,
                Iterable[bins_items.SequenceNormCoverage],
                mfb_views.Stats,
                Path,
            ]
        ],
        mfb_fs.NumberOfFlowManager,
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


# def plasbin_panassembly(
#     panasm_graph: gfapy.Gfa,
#     seed_fragments: Iterable[str],
#     gc_intervals: gc_items.Intervals,
#     fragment_gc_scores: Iterable[gc_items.SequenceGCScores],
#     fragment_plasmidness: Iterable[tuple[str, float]],
#     plasbin_config: pb_cfg.Binning,
#     gurobi_config: pb_lp_cfg.Gurobi,
#     output_directory: Path,
# ) -> Iterator[
#     tuple[
#         bins_items.Stats,
#         Iterable[bins_items.SequenceNormCoverage],
#         classbin_lp_views.ClassifyStats,
#         Path,
#     ]
# ]:
#     """Bin fragments of a pan-assembly graph.

#     Yields
#     ------
#     bins_items.Stats
#         Stats of the bin
#     iterator on bins_items.FragmentNormCoverage
#         Fragments and their normalized coverage in the bin
#     lp_views.MGCLBStats
#         MILP stats
#     Path
#         Log file path

#     Warning
#     -------
#     The GFA pan-assembly graph will mute.
#     """
#     return plasbin(
#         net.Network.from_panasm_graph(
#             panasm_graph,
#             seed_fragments,
#             fragment_gc_scores,
#             fragment_plasmidness,
#             plasbin_config.sink_arcs_domain(),
#         ),
#         gc_intervals,
#         plasbin_config,
#         gurobi_config,
#         output_directory,
#     )


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


def plasbin(
    network: net.Network,
    gc_intervals: gc_items.Intervals,
    plasbin_config: cmn_cfg.Binning,
    gurobi_config: cmn_lp_cfg.Gurobi,
    output_directory: Path,
) -> Iterator[
    tuple[
        list[
            tuple[
                bins_items.Stats,
                Iterable[bins_items.SequenceNormCoverage],
                mfb_views.Stats,
                Path,
            ]
        ],
        mfb_fs.NumberOfFlowManager,
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
    cmn_lp_cfg.configurate_global_gurobi(gurobi_config)

    lp_io_manager = classbin_lp_io.Manager(output_directory)
    # REFACTOR create lp io manager dir here?

    model, bins_vars, flow_union_frag_vars, state_csts = mfb_model.new(
        network,
        gc_intervals,
        plasbin_config.min_flow(),
        0,  # BUG TMP PLASMIDNESS COEFFICIENT
        plasbin_config.min_cumulative_len(),
        plasbin_config.obj_fun_domain(),
    )

    obj_lb_cst: gp.Constr | None = None

    with rich_prog.Progress(console=CONSOLE) as progress:
        remaining_bins = len(network.seeds())
        binning_task = progress.add_task("Binning", total=remaining_bins)

        number_of_circular_flows = 0
        number_of_partially_circular_flows = 0

        bin_number = 0
        test_new_circular_bin = True
        best_multi_flow_objective: float | None = None
        while test_new_circular_bin:
            state_csts[bin_number].activate_bin()
            state_csts[bin_number].define_circular()

            number_of_circular_flows += 1
            number_of_flows = (
                number_of_circular_flows + number_of_partially_circular_flows
            )
            nb_flow_fs_manager = lp_io_manager.mfb_fs_manager().get_nb_flow_fs_manager(
                number_of_circular_flows,
                number_of_partially_circular_flows,
            )
            result = _solve_multi_flows(
                network,
                plasbin_config,
                model,
                bins_vars,
                number_of_flows,
                nb_flow_fs_manager,
                state_csts[bin_number].is_partially_circular(),
            )
            # Infeasible because: no circular or not better
            if result is None or (
                best_multi_flow_objective is not None
                # BUG TMP ROUND
                and round(model.ObjVal, 4) <= round(best_multi_flow_objective, 4)
            ):
                test_new_circular_bin = False
                number_of_circular_flows -= 1
            else:
                # FIXME avoid same obj value
                # Because of gap, it can be better while it is the same splitted bin
                model.write(str(lp_io_manager.milp_dir() / "model.sol"))
                milp_results_values, log_file = result
                best_multi_flow_objective = model.ObjVal
                last_circular_plasmidness_score = (
                    mfb_obj.circular_flows_plasmidness_score(
                        bins_vars,
                        network,
                        plasbin_config.obj_fun_domain(),
                        number_of_circular_flows,
                    ).getValue()
                )
                if obj_lb_cst is None:
                    obj_lb_cst = mfb_cst.objective_lower_bound(
                        model,
                        bins_vars,
                        flow_union_frag_vars,
                        network,
                        gc_intervals,
                        plasbin_config.obj_fun_domain(),
                        best_multi_flow_objective,
                    )
                else:
                    mfb_cst.update_binning_objective_lower_bound(
                        model,
                        obj_lb_cst,
                        best_multi_flow_objective,
                    )
                bins_results = []
                for k in range(number_of_flows):
                    milp_stats, milp_result_values = milp_results_values[k]
                    fragment_norm_coverages, norm_coverage = (
                        cmn_lp_res.fragment_norm_coverages(milp_result_values)
                    )
                    bins_results.append(
                        (
                            bins_items.Stats(
                                milp_result_values.cumulative_length(),
                                milp_result_values.total_flow(),
                                norm_coverage,
                                milp_result_values.gc_interval(),
                            ),
                            fragment_norm_coverages,
                            milp_stats,
                            log_file,
                        ),
                    )
                yield bins_results, nb_flow_fs_manager
                progress.update(
                    binning_task,
                    advance=1,
                )
                remaining_bins -= 1
                bin_number += 1
                # REFACTOR use better condition with object wrapper
                if bin_number == state_csts.max_number_of_flows():
                    test_new_circular_bin = False
                else:
                    # DOCU plasmidness score order for C
                    mfb_cst.plasmidness_score_order(
                        model,
                        bins_vars,
                        network,
                        plasbin_config.obj_fun_domain(),
                        (bin_number - 1, bin_number + 1),
                    )
                    # FIXME embete le gap

        if number_of_circular_flows:
            _LOGGER.info(
                "Number of circular flows: %d",
                number_of_circular_flows,
            )
            # DOCU plasmidness score for C lb
            mfb_cst.circular_flows_plasmidness_score_lb(
                model,
                bins_vars,
                network,
                plasbin_config.obj_fun_domain(),
                number_of_circular_flows,
                last_circular_plasmidness_score,
                # REFACTOR unbound, use state programming
            )
            # DOCU plasmidness score order for C
            # TODO add it in C loop?

        test_new_partially_circular_bin = bin_number < state_csts.max_number_of_flows()
        # if test_new_partially_circular_bin:
        #     mfb_cst.each_seed_must_be_in_at_least_one_bin(model, bins_vars, network)
        # BUG TMP
        if test_new_partially_circular_bin and number_of_circular_flows:
            model.read(str(lp_io_manager.milp_dir() / "model.sol"))
        while test_new_partially_circular_bin:
            state_csts[bin_number].activate_bin()
            state_csts[bin_number].define_partially_circular()
            number_of_partially_circular_flows += 1
            number_of_flows = (
                number_of_circular_flows + number_of_partially_circular_flows
            )
            nb_flow_fs_manager = lp_io_manager.mfb_fs_manager().get_nb_flow_fs_manager(
                number_of_circular_flows,
                number_of_partially_circular_flows,
            )
            result = _solve_multi_flows(
                network,
                plasbin_config,
                model,
                bins_vars,
                number_of_flows,
                nb_flow_fs_manager,
                state_csts[bin_number].is_partially_circular(),
            )
            # Infeasible because: no other solution or not better
            if result is None or (
                best_multi_flow_objective is not None
                # BUG TMP ROUND
                and round(model.ObjVal, 4) <= round(best_multi_flow_objective, 4)
            ):
                progress.update(binning_task, advance=remaining_bins)
                remaining_bins = 0
                number_of_partially_circular_flows -= 1
                test_new_partially_circular_bin = False
            else:
                milp_results_values, log_file = result
                best_multi_flow_objective = model.ObjVal
                if obj_lb_cst is None:
                    obj_lb_cst = mfb_cst.objective_lower_bound(
                        model,
                        bins_vars,
                        flow_union_frag_vars,
                        network,
                        gc_intervals,
                        plasbin_config.obj_fun_domain(),
                        best_multi_flow_objective,
                    )
                else:
                    mfb_cst.update_binning_objective_lower_bound(
                        model,
                        obj_lb_cst,
                        best_multi_flow_objective,
                    )
                bins_results = []
                for k in range(number_of_flows):
                    milp_stats, milp_result_values = milp_results_values[k]
                    fragment_norm_coverages, norm_coverage = (
                        cmn_lp_res.fragment_norm_coverages(milp_result_values)
                    )
                    bins_results.append(
                        (
                            bins_items.Stats(
                                milp_result_values.cumulative_length(),
                                milp_result_values.total_flow(),
                                norm_coverage,
                                milp_result_values.gc_interval(),
                            ),
                            fragment_norm_coverages,
                            milp_stats,
                            log_file,
                        ),
                    )
                yield bins_results, nb_flow_fs_manager
                progress.update(binning_task, advance=1)
                remaining_bins -= 1
                bin_number += 1
                if bin_number == state_csts.max_number_of_flows():
                    test_new_partially_circular_bin = False
                else:
                    # DOCU plasmidness score order for PC
                    mfb_cst.plasmidness_score_order(
                        model,
                        bins_vars,
                        network,
                        plasbin_config.obj_fun_domain(),
                        (bin_number - 1, bin_number + 1),
                    )

        _LOGGER.info(
            "Number of partially circular flows: %d",
            number_of_partially_circular_flows,
        )
    # FIXME possibly change
    _LOGGER.info("Find %u bins.", bin_number)


# REFACTOR return other than None (class containning log file e.g.)
def _solve_multi_flows(
    network: net.Network,
    plasbin_config: cmn_cfg.Binning,
    lp_model: gp.Model,
    bins_vars: list[mfb_var.BinVariables],
    number_of_flows: int,
    io_manager: mfb_fs.NumberOfFlowManager,
    is_partially_circular: bool,  # noqa: FBT001
) -> tuple[list[tuple[mfb_views.Stats, cmn_lp_res.Pangebin]], Path] | None:
    log_file = _solve_mfb_model(lp_model, io_manager)

    if lp_model.Status != gp.GRB.OPTIMAL:
        return None

    bin_results: list[tuple[mfb_views.Stats, cmn_lp_res.Pangebin]] = []
    for k in range(number_of_flows):
        result_stats = mfb_views.stats_from_opt_bin_vars(
            network,
            bins_vars[k],
            plasbin_config.obj_fun_domain(),
            is_partially_circular,
        )
        milp_result_values = cmn_lp_res.Pangebin.from_optimal_vars_without_gc_intervals(
            network,
            bins_vars[k].flows(),
            bins_vars[k].sub_frag(),
        )

        bin_results.append((result_stats, milp_result_values))
    return (bin_results, log_file)


def _solve_mfb_model(
    model: gp.Model,
    io_manager: mfb_fs.NumberOfFlowManager,
) -> Path:
    """Run model."""
    io_manager.dir().mkdir(parents=True)
    log_file = io_manager.gurobi_log()
    model.Params.LogFile = str(log_file)
    model.optimize()
    return log_file
