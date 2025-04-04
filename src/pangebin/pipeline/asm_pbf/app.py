"""Pangebin asm-pbf application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Annotated

import gfapy  # type: ignore[import-untyped]
import typer

import pangebin.gc_content.app as gc_app
import pangebin.gc_content.input_output as gc_io
import pangebin.gc_content.items as gc_items
import pangebin.gfa.app as gfa_app
import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.logging as common_log
import pangebin.pbf_comp.app as pbf_comp_app
import pangebin.pbf_comp.input_output as pbf_comp_io
import pangebin.plasbin.bins.input_output as bins_io
import pangebin.plasbin.bins.items as bins_items
import pangebin.plasbin.config as pb_cfg
import pangebin.plasbin.input_output as pb_io
import pangebin.plasbin.milp.config as pb_lp_cfg
import pangebin.plasbin.once.create as once_create
import pangebin.plasbin.once.milp.input_output as once_lp_io
import pangebin.plasbin.once.milp.views as once_lp_views

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(
    name="asm-pbf",
    help="PangeBin asm-pbf pipeline",
    rich_markup_mode="rich",
)


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


class Options:
    """Pangebin asm-pbf options."""

    __RICH_HELP_PANEL = "Configurations"

    SKESA_GFA = typer.Option(
        help="SKESA GFA file flag",
        rich_help_panel=__RICH_HELP_PANEL,
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


@APP.command()
def once(
    assembly_gfa: Annotated[Path, Arguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, Arguments.SEED_CONTIGS_TSV],
    contig_plasmidness_tsv: Annotated[Path, Arguments.CONTIG_PLASMIDNESS_TSV],
    # pbf options
    is_skesa: Annotated[bool, Options.SKESA_GFA] = False,
    # Binning options
    binning_cfg_yaml: Annotated[
        Path | None,
        pb_cfg.BinningOptions.CONFIG_FILE,
    ] = None,
    # MILP options
    gurobi_cfg_yaml: Annotated[Path | None, pb_lp_cfg.GurobiOptions.CONFIG_FILE] = None,
    # IO options
    outdir: Annotated[
        Path,
        IOOptions.OUTPUT_DIR,
    ] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Run Pangebin once approach."""
    common_log.init_logger(_LOGGER, "Running pangebin once approach.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    std_gfa = _prepare_graph(assembly_gfa, is_skesa, debug)

    gc_scores_tsv = outdir / "gc_scores.tsv"
    gc_app.from_gfa(assembly_gfa, gc_scores_tsv, debug=debug)

    io_manager = _init_io_manager(outdir)

    binning_config = _init_binning_cfg(binning_cfg_yaml)
    gurobi_config = _init_gurobi_cfg(gurobi_cfg_yaml)

    seeds, intervals, gc_scores, plasmidness = _init_plasbin_input(
        seed_contigs_tsv,
        gc_scores_tsv,
        contig_plasmidness_tsv,
    )

    for k, (bin_stats, seq_normcovs, all_milp_stats, log_file) in enumerate(
        once_create.plasbin_assembly(
            std_gfa,
            seeds,
            intervals,
            gc_scores,
            plasmidness,
            binning_config,
            gurobi_config,
            io_manager.config().output_directory(),
        ),
    ):
        _write_outputs(
            io_manager,
            k,
            bin_stats,
            seq_normcovs,
            all_milp_stats,
            log_file,
        )

    pbf_comp_app.bins(
        io_manager.config().output_directory(),
        io_manager.config().output_directory() / "bins.tsv",
    )


ARG_CONFIGS_OUTDIR = typer.Argument(help="Output directory")


def write_once_configs(output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR]) -> None:
    """Write the configuration files for the once approach."""
    output_directory.mkdir(parents=True, exist_ok=True)

    binning_config = pb_cfg.Binning.default()
    binning_config.to_yaml(output_directory / pb_cfg.Binning.DEFAULT_YAML_FILE)

    gurobi_config = pb_lp_cfg.Gurobi.default()
    gurobi_config.to_yaml(output_directory / pb_lp_cfg.Gurobi.DEFAULT_YAML_FILE)


def _init_io_manager(outdir: Path) -> pb_io.Manager:
    outdir.mkdir(parents=True, exist_ok=True)
    return pb_io.Manager(pb_io.Config(output_directory=outdir))


def _prepare_graph(gfa_path: Path, is_skesa: bool, debug: bool) -> gfapy.Gfa:
    if is_skesa:
        gfa_app.fix_skesa(gfa_path, gfa_path, debug=debug)
    gfa = gfa_io.from_file(gfa_path)
    gfa_ops.set_segment_length_tags(gfa)
    if is_skesa:
        gfa_ops.convert_kmer_coverage_to_normalized_coverage(gfa)
    return gfa


def _init_binning_cfg(
    binning_cfg_yaml: Path | None,
) -> pb_cfg.Binning:
    binning_config = (
        pb_cfg.Binning.from_yaml(binning_cfg_yaml)
        if binning_cfg_yaml is not None
        else pb_cfg.Binning.default()
    )
    _LOGGER.debug("Binning config:\n%s", binning_config.to_dict())
    return binning_config


def _init_gurobi_cfg(
    gurobi_cfg_yaml: Path | None,
) -> pb_lp_cfg.Gurobi:
    gurobi_config = (
        pb_lp_cfg.Gurobi.from_yaml(gurobi_cfg_yaml)
        if gurobi_cfg_yaml is not None
        else pb_lp_cfg.Gurobi.default()
    )
    _LOGGER.debug("Gurobi config:\n%s", gurobi_config.to_dict())
    return gurobi_config


def _init_plasbin_input(
    seed_sequences_tsv: Path,
    sequence_gc_scores_tsv: Path,
    sequence_plasmidness_tsv: Path,
) -> tuple[
    list[str],
    gc_items.Intervals,
    list[gc_items.SequenceGCScores],
    list[tuple[str, float]],
]:
    with pbf_comp_io.SeedReader.open(seed_sequences_tsv) as seeds_fin:
        seeds = list(seeds_fin)
    with gc_io.ScoresReader.open(sequence_gc_scores_tsv) as gc_scores_fin:
        intervals = gc_scores_fin.intervals()
        gc_scores = list(gc_scores_fin)
    with pbf_comp_io.PlmReader.open(sequence_plasmidness_tsv) as plasmidness_fin:
        plasmidness = list(plasmidness_fin)
    return seeds, intervals, gc_scores, plasmidness


def _write_outputs(
    io_manager: pb_io.Manager,
    k: int,
    bin_stats: bins_items.Stats,
    seq_normcovs: Iterable[bins_items.SequenceNormCoverage],
    milp_stats: once_lp_views.MGCLBStats,
    log_file: Path,
) -> None:
    io_manager.bin_directory(k).mkdir(parents=True, exist_ok=True)
    bin_stats.to_yaml(io_manager.bin_stats_path(k))
    with bins_io.Writer.open(io_manager.bin_seq_normcov_path(k)) as fout:
        for seq_normcov in seq_normcovs:
            fout.write_sequence_normcov(
                seq_normcov.identifier(),
                seq_normcov.normalized_coverage(),
            )
    milp_stats.to_yaml(io_manager.milp_stats_path(k))
    io_manager.move_gurobi_logs(
        [log_file],
        once_lp_io.Manager.attributes_from_gurobi_log_path,
    )
