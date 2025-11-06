"""Common application items for the asm-pbf pipeline."""

import logging
from pathlib import Path

import gfapy  # type: ignore[import-untyped]
import typer

import pangebin.gc_content.app as gc_app
import pangebin.gc_content.create as gc_create
import pangebin.gc_content.input_output as gc_io
import pangebin.gc_content.items as gc_items
import pangebin.gfa.app as gfa_app
import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.pbf_comp.input_output as pbf_comp_io
import pangebin.pbf_comp.ops as pbf_comp_ops
import pangebin.plasbin.milp.config as cmn_lp_cfg

_LOGGER = logging.getLogger(__name__)


class Arguments:
    """Pangebin asm-pbf arguments."""

    ASSEMBLY_GFA = typer.Argument(
        help="Assembly GFA file",
    )

    SEED_CONTIGS_TSV = typer.Argument(
        help="TSV file with the seed contigs",
    )

    CONTIG_PLASMIDNESS_TSV = typer.Argument(
        help="TSV file with the contigs and their plasmidness scores",
    )


class InputOptions:
    """Pangebin asm-pbf input options."""

    _RICH_HELP_PANEL = "Input options"

    SKESA_GFA = typer.Option(
        help="SKESA GFA file flag",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    GC_SCORES_TSV = typer.Option(
        "--gc-scores",
        help="TSV file containing the GC scores",
        rich_help_panel=_RICH_HELP_PANEL,
    )


class GurobiOptions:
    """Gurobi options."""

    _RICH_HELP_PANEL = "Gurobi options"

    GUROBI_CFG_YAML = typer.Option(
        "--gurobi-cfg",
        help="The Gurobi configuration YAML path",
        rich_help_panel=_RICH_HELP_PANEL,
    )


class IOOptions:
    """Input-output options."""

    _RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        help="Output directory",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    BIN_RESULT_TSV = typer.Argument(
        help="Bin result TSV file",
        rich_help_panel=_RICH_HELP_PANEL,
    )


def prepare_graph(
    gfa_path: Path,
    is_skesa: bool,  # noqa: FBT001
    debug: bool,  # noqa: FBT001
) -> gfapy.Gfa:
    """Prepare the input GFA graph.

    Fix SKESA GFA if necessary.
    """
    if is_skesa:
        gfa_app.fix_skesa(gfa_path, gfa_path, debug=debug)
    gfa = gfa_io.from_file(gfa_path)
    gfa_ops.set_segment_length_tags(gfa)
    if is_skesa:
        gfa_ops.convert_kmer_coverage_to_normalized_coverage(gfa)
    return gfa


def prepare_gc_scores_tsv(
    assembly_gfa: Path,
    gc_scores_tsv: Path | None,
) -> tuple[gc_items.Intervals, list[gc_items.SequenceGCScores]]:
    """Prepare the GC scores TSV file."""
    if gc_scores_tsv is None:
        with gc_io.IntervalStepsReader.open(
            gc_app.DEFAULT_GC_CONTENT_INTERVAL_FILE,
        ) as interval_steps_reader:
            intervals = gc_items.Intervals(iter(interval_steps_reader))
        return intervals, list(
            gc_create.gfa_file_to_gc_scores(
                assembly_gfa,
                intervals,
            ),
        )
    with gc_io.ScoresReader.open(gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    return intervals, gc_scores


def init_seeds(
    seed_sequences_tsv: Path,
) -> list[str]:
    """Initialize the seed contigs."""
    with pbf_comp_io.SeedReader.open(seed_sequences_tsv) as seeds_fin:
        return list(seeds_fin)


def init_plasmidness(
    sequence_plasmidness_tsv: Path,
) -> list[tuple[str, float]]:
    """Initialize the plasmidness for the contigs."""
    with pbf_comp_io.PlmReader.open(sequence_plasmidness_tsv) as plasmidness_fin:
        # DOCU consider the plasmidness file gives originally plm between 0 and 1
        # DOCU PBF plasmidness is converted runtime to pangebin plasmidness
        plasmidness = [
            (seq_id, pbf_comp_ops.pbf_to_pg_plasmidness(plm))
            for seq_id, plm in plasmidness_fin
        ]
    # BUG TMP TEST PLASMIDNESS STEP
    _step = 0.001
    _corrected_plasmidness: list[tuple[str, float]] = []
    for ctg_id, plm in plasmidness:
        q, r = divmod(plm, _step)
        if r <= 0:
            new_pl = q * _step if r > -_step / 2 else (q - 1) * _step
        elif r < _step / 2:
            new_pl = q * _step
        else:
            new_pl = (q + 1) * _step
        _corrected_plasmidness.append((ctg_id, new_pl))

    return plasmidness  # _corrected_plasmidness  # BUG BEFORE: plasmidness


def init_gurobi_cfg(gurobi_cfg_yaml: Path | None) -> cmn_lp_cfg.Gurobi:
    """Initialize the Gurobi configuration."""
    gurobi_cfg = (
        cmn_lp_cfg.Gurobi.from_yaml(gurobi_cfg_yaml)
        if gurobi_cfg_yaml is not None
        else cmn_lp_cfg.Gurobi.default()
    )
    _LOGGER.debug("Gurobi config:\n%s", gurobi_cfg.to_dict())
    return gurobi_cfg
