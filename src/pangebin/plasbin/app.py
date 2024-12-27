"""PlasBin-Flow application module."""


# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import networkx as nx  # type: ignore[import-untyped]
import pandas as pd
import typer
from gurobipy import GRB, LinExpr, Model

from pangebin.plasbin import model_setup
from pangebin.plasbin.data_utils import (
    ASS_PENALTY_KEY,  # REFACTOR config in utils
    COV_KEY,
    DEFAULT_HEAD_STR,
    DEFAULT_MIN_CTG_LEN,
    DEFAULT_PLS_SCORE,
    DEFAULT_SEED_LEN_THRESHOLD,
    DEFAULT_SEED_SCORE_THRESHOLD,
    DEFAULT_SINK,
    DEFAULT_SOURCE,
    DEFAULT_TAIL_STR,
    LEN_KEY,
    SCORE_KEY,
    SEED_KEY,
    Assembler,  # REFACTOR config in gfa_fasta_utils
    BinMode,
    get_capacities,  # REFACTOR graph ops in utils
    get_seeds,
    log_data,  # REFACTOR log ops in utils
    read_ctgs_data,  # REFACTOR io ops in utils
    read_gc_data,
    read_links_data,
)
from pangebin.plasbin.log_errors_utils import (
    check_file,
    check_number_range,  # REFACTOR io ops in log_errors_utils
)

APP = typer.Typer(rich_markup_mode="rich")

SOURCE = DEFAULT_SOURCE
SINK = DEFAULT_SINK

MAX_THREADS = 32

DEFAULT_SCORE_OFFSET = 0.5
DEFAULT_ALPHA1 = 1
DEFAULT_ALPHA2 = 1
DEFAULT_ALPHA3 = 1
DEFAULT_ALPHA4 = 1
DEFAULT_RMITER_MAX = 50
DEFAULT_MIN_PLS_LEN = 1500
DEFAULT_GUROBI_MIP_GAP = 0.05
DEFAULT_GUROBI_TIME_LIMIT = 120


@dataclass
class PlasBinFlowArgs:
    """PlasBin-Flow args and opts."""

    __INPUT_CAT = "Input options"

    ARG_ASM_GRAPH_FILE = typer.Argument(
        help="Path to assembly graph file",
        rich_help_panel=__INPUT_CAT,
    )
    ARG_GC_PROB_FILE = typer.Argument(
        help="Path to GC probabilities file",
        rich_help_panel=__INPUT_CAT,
    )
    ARG_PLASMID_SCORE_FILE = typer.Argument(
        help="Path to plasmid score file",
        rich_help_panel=__INPUT_CAT,
    )
    ARG_ASSEMBLER = typer.Argument(
        help="Name of the assembler",
        rich_help_panel=__INPUT_CAT,
    )

    __OUTPUT_CAT = "Output options"

    ARG_OUT_DIR = typer.Option(
        help="Path to output dir",
        rich_help_panel=__OUTPUT_CAT,
    )  # REFACTOR should be an option
    ARG_OUT_FILE = typer.Option(
        help="Name of output file",
        rich_help_panel=__OUTPUT_CAT,
    )  # REFACTOR should not be an option and depend on ARG_OUT_DIR
    ARG_LOG_FILE = typer.Option(
        help="Path to log file",
        rich_help_panel=__OUTPUT_CAT,
    )  # REFACTOR should not be an option and depend on ARG_OUT_DIR

    __OBJECTIVE_FUNCTION_CAT = "Objective function parameters"

    OPT_PLASMID_SCORE_OFFSET = typer.Option(
        help="Offset for the plasmid score term",
        rich_help_panel=__OBJECTIVE_FUNCTION_CAT,
    )
    OPT_ALPHA_I = typer.Option(
        help="Weight of flow term",
        rich_help_panel=__OBJECTIVE_FUNCTION_CAT,
    )
    OPT_ALPHA_II = typer.Option(
        help="Weight of GC content term",
        rich_help_panel=__OBJECTIVE_FUNCTION_CAT,
    )
    OPT_ALPHA_III = typer.Option(
        help="Weight of plasmid score term",
        rich_help_panel=__OBJECTIVE_FUNCTION_CAT,
    )
    OPT_ALPHA_IV = typer.Option(
        help="Weight of assembly PENALTY term",
        rich_help_panel=__OBJECTIVE_FUNCTION_CAT,
    )

    __SEED_CAT = "Seed contig threshold parameters"

    OPT_SEED_LEN = typer.Option(
        help="Seed length threshold",
        rich_help_panel=__SEED_CAT,
    )
    OPT_SEED_SCORE = typer.Option(
        help="Seed plasmid score threshold",
        rich_help_panel=__SEED_CAT,
    )

    __GUROBI_CAT = "Gurobi stopping parameters"

    OPT_GUROBI_MIP_GAP = typer.Option(
        help="MIPGap parameter for Gurobi",
        rich_help_panel=__GUROBI_CAT,
    )
    OPT_GUROBI_TIME_LIMIT = typer.Option(
        help="Time limit for Gurobi (in seconds)",
        rich_help_panel=__GUROBI_CAT,
    )

    __OTHER_CAT = "Other parameters for miscellaneous purposes"

    # REFACTOR should rename this "Max plasmid search iterations"?
    OPT_RM_ITER_MAX = typer.Option(
        help="Number of iterations to remove circular components",
        rich_help_panel=__OTHER_CAT,
    )
    OPT_GC_INTERVALS = typer.Option(
        help="GC intervals file",
        rich_help_panel=__OTHER_CAT,
    )
    OPT_MIN_PLASMID_LENGTH = typer.Option(
        help="Minimum plasmid length to be reported",
        rich_help_panel=__OTHER_CAT,
    )
    OPT_MIN_CONTIG_LENGTH = typer.Option(
        help="Minimum contig length to account for plasmid score",
        rich_help_panel=__OTHER_CAT,
    )
    OPT_DEFAULT_PLASMID_SCORE = typer.Option(
        help="Default plasmid score",
        rich_help_panel=__OTHER_CAT,
    )
    OPT_BINNING_MODE = typer.Option(
        help="Binning mode (standard[std] / overlap[ovl] / naive[nve])",
        rich_help_panel=__OTHER_CAT,
    )
    OPT_CONTIG_PENALTY = typer.Option(
        help="Consider penalty for contigs",
        rich_help_panel=__OTHER_CAT,
    )


@APP.command()
def plasbin(
    # Input
    asm_graph: Annotated[Path, PlasBinFlowArgs.ARG_ASM_GRAPH_FILE],
    gc_prob: Annotated[Path, PlasBinFlowArgs.ARG_GC_PROB_FILE],
    plasmid_scores: Annotated[Path, PlasBinFlowArgs.ARG_PLASMID_SCORE_FILE],
    assembler: Annotated[
        Assembler,
        PlasBinFlowArgs.ARG_ASSEMBLER,
    ],
    # Output
    out_dir: Annotated[Path, PlasBinFlowArgs.ARG_OUT_DIR],
    out_file: Annotated[Path, PlasBinFlowArgs.ARG_OUT_FILE],
    log_file: Annotated[Path, PlasBinFlowArgs.ARG_LOG_FILE],
    # Objective function
    plasmid_score_offset: Annotated[
        float,
        PlasBinFlowArgs.OPT_PLASMID_SCORE_OFFSET,
    ] = DEFAULT_SCORE_OFFSET,
    alpha_i: Annotated[int, PlasBinFlowArgs.OPT_ALPHA_I] = DEFAULT_ALPHA1,
    alpha_ii: Annotated[int, PlasBinFlowArgs.OPT_ALPHA_II] = DEFAULT_ALPHA2,
    alpha_iii: Annotated[int, PlasBinFlowArgs.OPT_ALPHA_III] = DEFAULT_ALPHA3,
    alpha_iv: Annotated[int, PlasBinFlowArgs.OPT_ALPHA_IV] = DEFAULT_ALPHA4,
    seed_len_threshold: Annotated[
        int,
        PlasBinFlowArgs.OPT_SEED_LEN,
    ] = DEFAULT_SEED_LEN_THRESHOLD,
    seed_score_threshold: Annotated[
        float,
        PlasBinFlowArgs.OPT_SEED_SCORE,
    ] = DEFAULT_SEED_SCORE_THRESHOLD,
    # Gurobi
    gurobi_mip_gap: Annotated[
        float,
        PlasBinFlowArgs.OPT_GUROBI_MIP_GAP,
    ] = DEFAULT_GUROBI_MIP_GAP,
    gurobi_time_limit: Annotated[
        float,
        PlasBinFlowArgs.OPT_GUROBI_TIME_LIMIT,
    ] = DEFAULT_GUROBI_TIME_LIMIT,
    # Other
    rm_iter_max: Annotated[  # REFACTOR rename it to `max_plasmid_search_iterations`?
        int,
        PlasBinFlowArgs.OPT_RM_ITER_MAX,
    ] = DEFAULT_RMITER_MAX,
    gc_intervals: Annotated[Path | None, PlasBinFlowArgs.OPT_GC_INTERVALS] = None,
    min_plasmid_length: Annotated[
        int,
        PlasBinFlowArgs.OPT_MIN_PLASMID_LENGTH,
    ] = DEFAULT_MIN_PLS_LEN,
    min_contig_length: Annotated[
        int,
        PlasBinFlowArgs.OPT_MIN_CONTIG_LENGTH,
    ] = DEFAULT_MIN_CTG_LEN,
    default_plasmid_score: Annotated[
        float,
        PlasBinFlowArgs.OPT_DEFAULT_PLASMID_SCORE,
    ] = DEFAULT_PLS_SCORE,
    binning_mode: Annotated[
        BinMode,
        PlasBinFlowArgs.OPT_BINNING_MODE,
    ] = BinMode.STANDARD,
    contig_penalty: Annotated[
        bool,
        PlasBinFlowArgs.OPT_CONTIG_PENALTY,
    ] = True,  # REFACTOR reverse the bool
):
    """PlasBin-Flow."""
    # Naming and creating output files
    out_dir.mkdir(parents=True, exist_ok=True)
    # REFACTOR do good logging
    # Initialize logging
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(name)s - %(levelname)s - %(message)s",
    )
    # REFACTOR change checking method
    # Checking the values of parameters
    check_number_range(
        plasmid_score_offset,
        (0.0, 1.0),
        msg=f'INPUT\tParameter "-p", {plasmid_score_offset}: ',
    )
    check_number_range(
        alpha_i,
        (0.0, None),
        msg=f'INPUT\tParameter "-alpha_i", {alpha_i}: ',
    )
    check_number_range(
        alpha_ii,
        (0.0, None),
        msg=f'INPUT\tParameter "-alpha_ii", {alpha_ii}: ',
    )
    check_number_range(
        alpha_iii,
        (0.0, None),
        msg=f'INPUT\tParameter "-alpha_iii", {alpha_iii}: ',
    )
    check_number_range(
        alpha_iv,
        (0.0, None),
        msg=f'INPUT\tParameter "-alpha_iv", {alpha_iv}: ',
    )

    check_number_range(
        seed_len_threshold,
        (0, None),
        msg=f'INPUT\tParameter "-seed_len", {seed_len_threshold}: ',
    )
    check_number_range(
        seed_score_threshold,
        (0.0, 1.0),
        msg=f'INPUT\tParameter "-seed_score", {seed_score_threshold}: ',
    )
    check_number_range(
        min_plasmid_length,
        (0.0, None),
        msg=f'INPUT\tParameter "-min_pls_len", {min_plasmid_length}: ',
    )
    check_number_range(
        min_contig_length,
        (0.0, None),
        msg=f'INPUT\tParameter "-min_ctg_len", {min_contig_length}: ',
    )
    check_number_range(
        default_plasmid_score,
        (0.0, 1.0),
        msg=f'INPUT\tParameter "-default_pls_score", {default_plasmid_score}: ',
    )
    check_number_range(
        gurobi_mip_gap,
        (0.0, 1.0),
        msg=f'INPUT\tParameter "-gurobi_mip_gap", {gurobi_mip_gap}: ',
    )
    check_number_range(
        gurobi_time_limit,
        (0.0, None),
        msg=f'INPUT\tParameter "-gurobi_time_limit", {gurobi_time_limit}: ',
    )

    # Checking that input files exist and are not empty (warning if empty)
    input_files = [asm_graph, plasmid_scores, gc_prob]
    if gc_intervals is not None:
        input_files.append(gc_intervals)
    for in_file in input_files:
        check_file(in_file)

    # Reading and checking data
    contigs_dict = read_ctgs_data(
        asm_graph,
        plasmid_scores,
        assembler,
        default_pls_score=default_plasmid_score,
        min_ctg_len=min_contig_length,
    )
    seeds_set = get_seeds(
        contigs_dict,
        seed_len=seed_len_threshold,
        seed_score=seed_score_threshold,
    )
    gc_probs, gc_pens = read_gc_data(
        gc_prob,
        gc_intervals,
        contigs_dict.keys(),
    )
    links_list = read_links_data(asm_graph)
    log_data(
        contigs_dict,
        links_list,
        asm_graph,
        plasmid_scores,
    )

    output_bins = out_file.open("w")

    # -----------------------------------------------
    # Main program
    """
    Dictionary for storing bins
    Key = Plasmid id, Value = Dictionary of attributes
    Flow (float), GC bin id, List of contigs with multiplicities
    """
    pbf_bins = {}
    n_bins = 0

    n_comp = 0
    while len(seeds_set) > 0:
        """
        For consistency, both extremities for a contig should be part of exactly the same number of links in a plasmid.
        For each extremity, we make a list of links that involve the extremity.
        """
        extr_dict = {}  # Key: Extremity (e.g. ('1',DEFAULT_HEAD_STR)), Value: List of unordered links incident on extremity

        """
        List of incoming and outgoing edges for each extremity
        Including links from S and links to T
        """
        incoming, outgoing = (
            {},
            {},
        )  # Key: Extremity, Value: List of ordered links in/out of extremity

        contigs_df = pd.DataFrame.from_dict(contigs_dict).T

        max_cov = max(contigs_df[COV_KEY])
        min_cov = min(contigs_df[COV_KEY])

        for c in contigs_dict:
            ext1, ext2 = (c, DEFAULT_HEAD_STR), (c, DEFAULT_TAIL_STR)
            extr_dict[ext1], extr_dict[ext2] = [], []
            incoming[ext1], incoming[ext2] = [], []
            outgoing[ext1], outgoing[ext2] = [], []
            if c in seeds_set:
                extr_dict[ext1], extr_dict[ext2] = [(SOURCE, ext1)], [(SOURCE, ext2)]
                incoming[ext1], incoming[ext2] = [(SOURCE, ext1)], [(SOURCE, ext2)]
                extr_dict[ext1], extr_dict[ext2] = [(ext1, SINK)], [(ext2, SINK)]
            outgoing[ext1], outgoing[ext2] = [(ext1, SINK)], [(ext2, SINK)]

            for link in links_list:
                if link[0] == ext1 or link[1] == ext1:
                    extr_dict[ext1].append(link)
                    extr_dict[ext1].append(link[::-1])
                ## ELSEIF?
                if link[0] == ext2 or link[1] == ext2:
                    extr_dict[ext2].append(link)
                    extr_dict[ext2].append(link[::-1])

        for link in links_list:
            ext1, ext2 = link[0], link[1]
            incoming[ext1].append(link[::-1])
            outgoing[ext1].append(link)
            incoming[ext2].append(link)
            outgoing[ext2].append(link[::-1])

        capacities = get_capacities(links_list, contigs_dict)

        logging.info(f"Number of vertices: {len(contigs_dict)}")
        logging.info(f"Number of edges: {len(links_list)}")

        # -----------------------------------------------
        # Initializing the ILP
        m = Model("Plasmids")
        m.params.Threads = MAX_THREADS
        m.params.LogFile = str(out_dir / "m.log")  # REFACTOR gurobi log
        m.setParam(GRB.Param.TimeLimit, gurobi_time_limit)
        m.setParam(GRB.Param.MIPGap, gurobi_mip_gap)

        # Initializing variables
        contigs = {}  # Key: Contig (e.g. '1'), Value: Gurobi binary variable
        contigs = model_setup.contig_vars(m, contigs_dict, contigs)

        links = {}  # Key: Directed link from one extremity to another
        # (e.g. (('1',DEFAULT_HEAD_STR),('2',DEFAULT_TAIL_STR)) ), Value: Gurobi binary variable
        links = model_setup.link_vars(m, links_list, links, contigs)

        plas_GC = {}  # Key: GC bin, Value: Gurobi binary variable
        contig_GC = {}  # Nested Dict: Key: Contig, Value: Dict of GC bins (Key: GC bin, Value: Gurobi binary variable)
        plas_GC, contig_GC = model_setup.GC_vars(m, gc_probs, plas_GC, contig_GC)

        flows = {}  # Key: Directed link, Value: Gurobi continuous variable
        counted_F = {}  # Key: Directed link, Value: Gurobi continuous variable
        flows, counted_F = model_setup.flow_vars(m, links, flows, counted_F)
        flow_var = m.addVar(vtype=GRB.CONTINUOUS, name="overall-flow")

        # -----------------------------------------------
        # Setting up the expression for the objective function
        expr = LinExpr()

        expr.addTerms(alpha_i, flow_var)
        for c in contigs:
            for b in plas_GC:
                expr.addTerms(alpha_ii * (gc_pens[c][b]), contig_GC[c][b])
            expr.addTerms(
                alpha_iii * (contigs_dict[c][SCORE_KEY] - plasmid_score_offset),
                contigs[c],
            )
            if contig_penalty:
                assembly_penalty = alpha_iv * contigs_dict[c][ASS_PENALTY_KEY]
                assembly_penalty = min(assembly_penalty, 1.0)
                expr.addTerms(-assembly_penalty, contigs[c])

        m.setObjective(expr, GRB.MAXIMIZE)

        # -----------------------------------------------
        # Setting up constraints

        constraint_count = 0

        # Constraint type 1: A link 'e' is in the plasmid only if both its endpoints are in the plasmid.
        m, constraint_count = model_setup.link_inclusion_constr(
            m,
            links,
            contigs,
            constraint_count,
        )

        # Constraint type 2: An extremity is in the plasmid only if at least one link is incident on it.
        m, constraint_count = model_setup.extr_inclusion_constr(
            m,
            links,
            contigs,
            extr_dict,
            constraint_count,
        )

        # Constraint type 3: A contig is considered to be a ”counted” seed if it is eligible to be a seed contig
        #                   and is considered to be part of the solution
        m, constraint_count = model_setup.seed_inclusion_constr(
            m,
            contigs,
            contigs_dict,
            constraint_count,
        )

        # Constraint type 4: 'F' should equal the flow out of SOURCE and into SINK.
        #                   Exactly one edge exits SOURCE and exactly one enters SINK.
        m, constraint_count = model_setup.min_flow_constr(
            m,
            links,
            flows,
            flow_var,
            min_cov,
            constraint_count,
        )

        # Constraint types 5 and 6
        # Conservation constraints: Flow into ('u',DEFAULT_HEAD_STR) (resp. ('u',DEFAULT_TAIL_STR) ) should be equal to flow out of ('u',DEFAULT_TAIL_STR) (resp. ('u',DEFAULT_HEAD_STR) ).
        # Capacity constraints    : The maximum flow into a vertex should be at most the capacity (read depth) of the vertex itself.
        # The maximum flow through an edge has to be at most the capacity (capacities[e]) of the edge.
        m, constraint_count = model_setup.flow_conservation_constraints(
            m,
            links,
            contigs,
            flows,
            incoming,
            outgoing,
            capacities,
            contigs_dict,
            constraint_count,
        )

        # Constraint types 7 and 8
        # 7. The overall flow 'F' through link 'e' is ”counted” only if 'e' is part of the solution.
        # 8. The overall flow 'F' cannot exceed the flow through any active link 'e'.
        m, constraint_count = model_setup.counted_flow_constr(
            m,
            links,
            flows,
            counted_F,
            flow_var,
            max_cov,
            constraint_count,
        )

        # Constraint type 9: Handling the GC-content term in the objective function
        m, constraint_count = model_setup.GC_constr(
            m,
            contig_GC,
            plas_GC,
            contigs,
            constraint_count,
        )

        extra_comps = 1  # default
        rmiter_count = 0
        dc_count = 0
        dc_dict = {}
        while extra_comps >= 1 and rmiter_count <= rm_iter_max:
            # Running the MILP
            start = time.time()
            m.optimize()
            stop = time.time()
            duration = stop - start
            logging.info(f"MILP\tIteration {n_bins}: {duration}")

            rmiter_count += 1

            # Message if solution not obtained
            if m.Status == GRB.Status.INFEASIBLE:
                logging.warning(
                    "MILP\tThe model cannot be solved because it is infeasible",
                )
            elif m.Status == GRB.Status.UNBOUNDED:
                logging.warning(
                    "MILP\tThe model cannot be solved because it is unbounded",
                )
            elif m.Status == GRB.Status.INF_OR_UNBD:
                logging.warning(
                    "MILP\tThe model cannot be solved because it is infeasible or unbounded ",
                )

            # Storing Irreducible Inconsistent Subsystem in case solution is not obtained
            if m.Status in (GRB.Status.INF_OR_UNBD, GRB.Status.INFEASIBLE):
                logging.warning(
                    "MILP\tStoring Irreducible Inconsistent Subsystem in m.ilp",
                )
                m.computeIIS()
                m.write("m.ilp")
                for con in m.getConstrs():
                    if con.IISConstr:
                        logging.warning(f"{con.ConstrName}")
                sys.exit(1)

            logging.info("Solution:")
            logging.info(m.printAttr("x"))  # REFACTOR log debug

            # Flow zero condition
            if flow_var.X == 0:
                sys.exit(1)

            # Finding components in the solution
            flow_graph = nx.DiGraph()
            for c in contigs:
                if contigs[c].X > 0:
                    flow_graph.add_node(c)
                flow_graph.add_node(SOURCE)
                flow_graph.add_node(SINK)

            for e in links:
                if links[e].X > 0:
                    end1, end2 = e[0], e[1]
                    if end1 == SOURCE:
                        c1 = SOURCE
                        ext1 = None
                    else:
                        c1 = end1[0]
                        ext1 = end1[1]
                    if end2 == SINK:
                        c2 = SINK
                        ext2 = None
                    else:
                        c2 = end2[0]
                        ext2 = end2[1]
                    flow_graph.add_edge(c1, c2)
                    nx.set_edge_attributes(
                        flow_graph,
                        {(c1, c2): {"extremities": (ext1, ext2)}},
                    )
            conn_comps = nx.weakly_connected_components(flow_graph)

            comp_count = 0
            if nx.number_weakly_connected_components(flow_graph) > 1:
                logging.info(
                    f"Disconnected component edges removal {rmiter_count} for plasmid bin {n_bins}",
                )
            for comp in conn_comps:
                comp_count += 1
                comp_len = 0
                for node in comp:
                    if node not in (SOURCE, SINK):
                        comp_len += contigs_dict[node][LEN_KEY]
                if SOURCE in comp:
                    ST_comp = flow_graph.subgraph(comp)
                else:
                    disconn_comp = flow_graph.subgraph(comp)
                    # Muting individual edges
                    for edge in disconn_comp.edges:
                        exts = nx.get_edge_attributes(flow_graph, "extremities")[edge]
                        e = ((edge[0], exts[0]), (edge[1], exts[1]))
                        m.addConstr(links[e] == 0, "muted_edge-" + str(e))
                        logging.info(f"{e!s}")

            # Condition to stop iterating: If number of connected components is 1, there are no extra components, thus breaking the while loop.
            extra_comps = comp_count - 1

        # Recording the plasmid bin if the plasmid is long enough
        plasmid_length = 0
        for c in contigs:
            if contigs[c].X > 0:
                plasmid_length += contigs_dict[c][LEN_KEY]
        logging.info(
            f"Plasmid length: {plasmid_length} and min_pls_len: {min_plasmid_length}",
        )
        if plasmid_length >= min_plasmid_length:
            # Recording objective function scores
            GC_sum = 0
            gd_sum = 0
            as_sum = 0
            for c in contigs:
                gc_c_sum = 0
                for b in contig_GC[c]:
                    GC_sum += alpha_ii * (gc_pens[c][b]) * contig_GC[c][b].X
                    gc_c_sum += alpha_ii * (gc_pens[c][b]) * contig_GC[c][b].X
                gd_sum += alpha_iii * (contigs_dict[c][SCORE_KEY] - 0.5) * contigs[c].X
                if contig_penalty:
                    as_sum -= (contigs_dict[c][ASS_PENALTY_KEY]) * contigs[c].X
                lr_gd = (contigs_dict[c][SCORE_KEY] - 0.5) * contigs[c].X

            # Recording components in the solution
            flow_graph = nx.DiGraph()
            for c in contigs:
                if contigs[c].X > 0:
                    flow_graph.add_node(c)
            flow_graph.add_node(SOURCE)
            flow_graph.add_node(SINK)

            for e in links:
                if links[e].X > 0:
                    end1, end2 = e[0], e[1]
                    if end1 == SOURCE:
                        c1 = SOURCE
                        ext1 = None
                    else:
                        c1 = end1[0]
                        ext1 = end1[1]
                    if end2 == SINK:
                        c2 = SINK
                        ext2 = None
                    else:
                        c2 = end2[0]
                        ext2 = end2[1]

                    flow_graph.add_edge(c1, c2)
                    nx.set_edge_attributes(
                        flow_graph,
                        {(c1, c2): {"extremities": (ext1, ext2)}},
                    )
            conn_comps = nx.weakly_connected_components(flow_graph)

            GC_bin = 0
            for b in plas_GC:
                if plas_GC[b].X == 1:
                    GC_bin = b

            comp_count = 0
            for comp in conn_comps:
                n_bins += 1
                pbf_bins[n_bins] = {
                    "Flow": flow_var.X,
                    "GC_bin": GC_bin,
                    "Length": plasmid_length,
                    "Contigs": {},
                }
                comp_count += 1
                comp_len = 0
                for node in comp:
                    if node not in (SOURCE, SINK):
                        if node not in pbf_bins[n_bins]["Contigs"]:
                            pbf_bins[n_bins]["Contigs"][node] = 0

        # ### BUG: the following "for" loop(s) fail(s) if the minimum plasmid length requirment is not met
        #     else:
        #         continue

        # Updating assembly graph and formulation
        for e in flows:
            if e[1] != SINK:
                c = e[1][0]
                contigs_dict[c][COV_KEY] = max(0, contigs_dict[c][COV_KEY] - flows[e].X)
                logging.info("debug!!")
                logging.info(
                    f"{pbf_bins}, {n_bins}, {c}, {contigs_dict[c][COV_KEY]}, {flows[e].X}, {flow_var.X}",
                )
                logging.info("!!debug")

                # if the minimum plasmid length requirment is not met, pbf_bins is empty
                if pbf_bins != {} and c in pbf_bins[n_bins]["Contigs"]:
                    pbf_bins[n_bins]["Contigs"][c] += round(flows[e].X / flow_var.X, 2)

        for c in contigs:
            if contigs[c].X > 0:
                if c in seeds_set and contigs_dict[c][COV_KEY] <= 0.5:
                    contigs_dict[c][SEED_KEY] = 0
                    seeds_set.remove(c)

                if contigs_dict[c][COV_KEY] <= 0.05:
                    del contigs_dict[c]
                    del gc_probs[c]
                    del gc_pens[c]

                    links_list = [
                        e
                        for e in links_list
                        if (c, DEFAULT_HEAD_STR) not in e
                        and (c, DEFAULT_TAIL_STR) not in e
                    ]

    output_bins.write("#Pls_ID\tFlow\tGC_bin\tLength\tContigs\n")
    if binning_mode == BinMode.STANDARD:
        for plastid in pbf_bins:
            fval = "{:.2f}".format(pbf_bins[plastid]["Flow"])
            gcb = pbf_bins[plastid]["GC_bin"]
            pl_length = pbf_bins[plastid]["Length"]

            output_bins.write(
                "P"
                + str(plastid)
                + "\t"
                + str(fval)
                + "\t"
                + str(gcb)
                + "\t"
                + str(pl_length)
                + "\t",
            )
            for nctg, c in enumerate(pbf_bins[plastid]["Contigs"]):
                ctg_mul = pbf_bins[plastid]["Contigs"][c]
                if nctg == 0:
                    output_bins.write(c + ":" + str(ctg_mul))
                else:
                    output_bins.write("," + c + ":" + str(ctg_mul))
            output_bins.write("\n")

    elif binning_mode == BinMode.OVERLAP:
        # extend_bins(ovl, pbf_bins)["Flow"]
        # extend_bins(pbf_bins)["GC_bin"]
        # extend_bins(pbf_bins)["Length"]
        # extend_bins(pbf_bins)["Contigs"][c] for loop
        pass
    elif binning_mode == BinMode.NAIVE:
        # extend_bins
        # extend_bins(nve, pbf_bins)["Flow"]
        # extend_bins(pbf_bins)["GC_bin"]
        # extend_bins(pbf_bins)["Length"]
        # extend_bins(pbf_bins)["Contigs"][c] for loop
        pass
