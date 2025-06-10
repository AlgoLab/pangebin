"""Pangebin asm-pbf application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

import logging
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Annotated

import gfapy  # type: ignore[import-untyped]
import typer

import pangebin.gc_content.app as gc_app
import pangebin.gc_content.create as gc_create
import pangebin.gc_content.input_output as gc_io
import pangebin.gc_content.items as gc_items
import pangebin.gfa.app as gfa_app
import pangebin.gfa.input_output as gfa_io
import pangebin.gfa.ops as gfa_ops
import pangebin.pbf_comp.app as pbf_comp_app
import pangebin.pbf_comp.input_output as pbf_comp_io
import pangebin.pbf_comp.ops as pbf_comp_ops
import pangebin.pblog as common_log
import pangebin.plasbin.binlab.config as binlab_cfg
import pangebin.plasbin.binlab.create as binlab_create
import pangebin.plasbin.binlab.milp.input_output as binlab_lp_io
import pangebin.plasbin.binlab.milp.views as binlab_lp_views
import pangebin.plasbin.bins.input_output as bins_io
import pangebin.plasbin.bins.items as bins_items
import pangebin.plasbin.classbin.create as classbin_create
import pangebin.plasbin.classbin.milp.views as classbin_lp_views
import pangebin.plasbin.classbin.multi_flow.milp.file_system as mfb_lp_fs
import pangebin.plasbin.classbin.multi_flow.milp.views as mfb_lp_views
import pangebin.plasbin.config as pb_cfg
import pangebin.plasbin.decomp.config as decomp_cfg
import pangebin.plasbin.decomp.create as decomp_create
import pangebin.plasbin.decomp.milp.input_output as decomp_lp_io
import pangebin.plasbin.decomp.milp.views as decomp_lp_views
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

    _RICH_HELP_PANEL = "Configurations"

    SKESA_GFA = typer.Option(
        help="SKESA GFA file flag",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    GC_SCORES_TSV = typer.Option(
        "--gc-scores",
        help="TSV file containing the GC scores",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    BINNING_CFG_YAML = typer.Option(
        "--bin-cfg",
        help="The general binning configuration YAML path",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    GUROBI_CFG_YAML = typer.Option(
        "--gurobi-cfg",
        help="The Gurobi configuration YAML path",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    DECOMP_CFG_YAML = typer.Option(
        "--decomp-cfg",
        help="The decomp approach configuration YAML path",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    BINLAB_CFG_YAML = typer.Option(
        "--binlab-cfg",
        help="The binlab approach configuration YAML path",
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


@APP.command(help="Run Pangebin decomp approach.")
def decomp(
    assembly_gfa: Annotated[Path, Arguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, Arguments.SEED_CONTIGS_TSV],
    contig_plasmidness_tsv: Annotated[Path, Arguments.CONTIG_PLASMIDNESS_TSV],
    # pbf options
    is_skesa: Annotated[bool, Options.SKESA_GFA] = False,
    gc_scores_tsv: Annotated[Path | None, Options.GC_SCORES_TSV] = None,
    # Binning options
    binning_cfg_yaml: Annotated[Path | None, Options.BINNING_CFG_YAML] = None,
    # Decomposition options
    decomp_cfg_yaml: Annotated[Path | None, Options.DECOMP_CFG_YAML] = None,
    # MILP options
    gurobi_cfg_yaml: Annotated[Path | None, Options.GUROBI_CFG_YAML] = None,
    # IO options
    outdir: Annotated[Path, IOOptions.OUTPUT_DIR] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Run Pangebin decomp approach.

    Returns
    -------
    Path
        Path to the PlasBin-flow bin file.
    """
    common_log.init_logger(_LOGGER, "Running pangebin decomp approach.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    std_gfa = _prepare_graph(assembly_gfa, is_skesa, debug)

    intervals, gc_scores = _prepare_gc_scores_tsv(assembly_gfa, gc_scores_tsv)

    io_manager = _init_io_manager(outdir)

    binning_config = _init_binning_cfg(binning_cfg_yaml)
    decomp_config = _init_decomp_cfg(decomp_cfg_yaml)
    gurobi_config = _init_gurobi_cfg(gurobi_cfg_yaml)

    seeds, plasmidness = _init_plasbin_input(seed_contigs_tsv, contig_plasmidness_tsv)

    for k, (bin_stats, seq_normcovs, all_milp_stats, log_file) in enumerate(
        decomp_create.plasbin_assembly(
            std_gfa,
            seeds,
            intervals,
            gc_scores,
            plasmidness,
            binning_config,
            decomp_config,
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
            decomp_lp_io.Manager.attributes_from_gurobi_log_path,
            log_file,
        )

    return _convert_pangebin_output_to_pbf_output(
        io_manager.config().output_directory(),
    )


@APP.command(help="Run Pangebin binlab approach.")
def binlab(
    assembly_gfa: Annotated[Path, Arguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, Arguments.SEED_CONTIGS_TSV],
    contig_plasmidness_tsv: Annotated[Path, Arguments.CONTIG_PLASMIDNESS_TSV],
    # pbf options
    is_skesa: Annotated[bool, Options.SKESA_GFA] = False,
    gc_scores_tsv: Annotated[Path | None, Options.GC_SCORES_TSV] = None,
    # Binning options
    binning_cfg_yaml: Annotated[Path | None, Options.BINNING_CFG_YAML] = None,
    # Decomposition options
    binlab_cfg_yaml: Annotated[Path | None, Options.DECOMP_CFG_YAML] = None,
    # MILP options
    gurobi_cfg_yaml: Annotated[Path | None, Options.GUROBI_CFG_YAML] = None,
    # IO options
    outdir: Annotated[
        Path,
        IOOptions.OUTPUT_DIR,
    ] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Run Pangebin binlab approach.

    Returns
    -------
    Path
        Path to the PlasBin-flow bin file.
    """
    common_log.init_logger(_LOGGER, "Running pangebin binlab approach.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    std_gfa = _prepare_graph(assembly_gfa, is_skesa, debug)

    intervals, gc_scores = _prepare_gc_scores_tsv(assembly_gfa, gc_scores_tsv)

    io_manager = _init_io_manager(outdir)

    binning_config = _init_binning_cfg(binning_cfg_yaml)
    binlab_config = _init_binlab_cfg(binlab_cfg_yaml)
    gurobi_config = _init_gurobi_cfg(gurobi_cfg_yaml)

    seeds, plasmidness = _init_plasbin_input(seed_contigs_tsv, contig_plasmidness_tsv)

    for k, (bin_stats, seq_normcovs, all_milp_stats, log_files) in enumerate(
        binlab_create.plasbin_assembly(
            std_gfa,
            seeds,
            intervals,
            gc_scores,
            plasmidness,
            binning_config,
            binlab_config,
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
            binlab_lp_io.Manager.attributes_from_gurobi_log_path,
            log_files,
        )

    pbf_comp_app.bins(
        io_manager.config().output_directory(),
        io_manager.config().output_directory() / "bins.tsv",
    )

    return _convert_pangebin_output_to_pbf_output(
        io_manager.config().output_directory(),
    )


@APP.command(help="Run Pangebin once approach.")
def once(
    assembly_gfa: Annotated[Path, Arguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, Arguments.SEED_CONTIGS_TSV],
    contig_plasmidness_tsv: Annotated[Path, Arguments.CONTIG_PLASMIDNESS_TSV],
    # pbf options
    is_skesa: Annotated[bool, Options.SKESA_GFA] = False,
    gc_scores_tsv: Annotated[Path | None, Options.GC_SCORES_TSV] = None,
    # Binning options
    binning_cfg_yaml: Annotated[Path | None, Options.BINNING_CFG_YAML] = None,
    # Decomposition options
    # MILP options
    gurobi_cfg_yaml: Annotated[Path | None, Options.GUROBI_CFG_YAML] = None,
    # IO options
    outdir: Annotated[
        Path,
        IOOptions.OUTPUT_DIR,
    ] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Run Pangebin once approach.

    Returns
    -------
    Path
        Path to the PlasBin-flow bin file.
    """
    common_log.init_logger(_LOGGER, "Running pangebin once approach.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    std_gfa = _prepare_graph(assembly_gfa, is_skesa, debug)

    intervals, gc_scores = _prepare_gc_scores_tsv(assembly_gfa, gc_scores_tsv)

    io_manager = _init_io_manager(outdir)

    binning_config = _init_binning_cfg(binning_cfg_yaml)
    gurobi_config = _init_gurobi_cfg(gurobi_cfg_yaml)

    seeds, plasmidness = _init_plasbin_input(seed_contigs_tsv, contig_plasmidness_tsv)

    # XXX PAER TESTS
    seeds_set = set(seeds)
    for k, (frag_id, plm) in enumerate(plasmidness):
        if frag_id in seeds_set:
            plasmidness[k] = (frag_id, (plm + 1) / 2)

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
            once_lp_io.Manager.attributes_from_gurobi_log_path,
            [log_file],
        )

    return _convert_pangebin_output_to_pbf_output(
        io_manager.config().output_directory(),
    )


@APP.command(help="Run Pangebin classbin approach.")
def classbin(
    assembly_gfa: Annotated[Path, Arguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, Arguments.SEED_CONTIGS_TSV],
    contig_plasmidness_tsv: Annotated[Path, Arguments.CONTIG_PLASMIDNESS_TSV],
    # pbf options
    is_skesa: Annotated[bool, Options.SKESA_GFA] = False,
    gc_scores_tsv: Annotated[Path | None, Options.GC_SCORES_TSV] = None,
    # Binning options
    binning_cfg_yaml: Annotated[Path | None, Options.BINNING_CFG_YAML] = None,
    # Decomposition options
    # MILP options
    gurobi_cfg_yaml: Annotated[Path | None, Options.GUROBI_CFG_YAML] = None,
    # IO options
    outdir: Annotated[
        Path,
        IOOptions.OUTPUT_DIR,
    ] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Run Pangebin classbin approach.

    Returns
    -------
    Path
        Path to the PlasBin-flow bin file.
    """
    # TODO WORK IN PROGRESS
    common_log.init_logger(_LOGGER, "Running pangebin classbin approach.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    std_gfa = _prepare_graph(assembly_gfa, is_skesa, debug)

    intervals, gc_scores = _prepare_gc_scores_tsv(assembly_gfa, gc_scores_tsv)

    io_manager = _init_io_manager(outdir)

    binning_config = _init_binning_cfg(binning_cfg_yaml)
    gurobi_config = _init_gurobi_cfg(gurobi_cfg_yaml)

    seeds, plasmidness = _init_plasbin_input(seed_contigs_tsv, contig_plasmidness_tsv)

    # XXX PAER TESTS
    # REFACTOR clarify this logic
    seeds_set = set(seeds)
    for k, (frag_id, plm) in enumerate(plasmidness):
        if frag_id in seeds_set:
            plasmidness[k] = (frag_id, (plm + 1) / 2)

    # XXX tmp fix
    _all_nb_flow_fs_managers: list[mfb_lp_fs.NumberOfFlowManager] = []
    for bin_results, nb_flow_fs_manager in classbin_create.plasbin_assembly(
        std_gfa,
        seeds,
        intervals,
        gc_scores,
        plasmidness,
        binning_config,
        gurobi_config,
        io_manager.config().output_directory(),
    ):
        _all_nb_flow_fs_managers.append(nb_flow_fs_manager)
        _write_mfb_outputs(
            nb_flow_fs_manager,
            bin_results,
        )

    if _all_nb_flow_fs_managers:
        _last_nb_flow_fs_manager = _all_nb_flow_fs_managers[-1]
        import shutil

        _nb_bin = 0
        while _last_nb_flow_fs_manager.bin_fs_manager(_nb_bin).dir().exists():
            shutil.copytree(
                _last_nb_flow_fs_manager.bin_fs_manager(_nb_bin).dir(),
                io_manager.bin_directory(_nb_bin),
            )
            _nb_bin += 1

    return _convert_pangebin_output_to_pbf_output(
        io_manager.config().output_directory(),
    )


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


def _prepare_gc_scores_tsv(
    assembly_gfa: Path,
    gc_scores_tsv: Path | None,
) -> tuple[gc_items.Intervals, list[gc_items.SequenceGCScores]]:
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


def _init_binning_cfg(binning_cfg_yaml: Path | None) -> pb_cfg.Binning:
    binning_config = (
        pb_cfg.Binning.from_yaml(binning_cfg_yaml)
        if binning_cfg_yaml is not None
        else pb_cfg.Binning.default()
    )
    _LOGGER.debug("Binning config:\n%s", binning_config.to_dict())
    return binning_config


def _init_decomp_cfg(decomp_cfg_yaml: Path | None) -> decomp_cfg.Decomp:
    decomp_config = (
        decomp_cfg.Decomp.from_yaml(decomp_cfg_yaml)
        if decomp_cfg_yaml is not None
        else decomp_cfg.Decomp.default()
    )
    _LOGGER.debug("Decomp config:\n%s", decomp_config.to_dict())
    return decomp_config


def _init_binlab_cfg(binlab_cfg_yaml: Path | None) -> binlab_cfg.Binlab:
    binlab_config = (
        binlab_cfg.Binlab.from_yaml(binlab_cfg_yaml)
        if binlab_cfg_yaml is not None
        else binlab_cfg.Binlab.default()
    )
    _LOGGER.debug("Binlab config:\n%s", binlab_config.to_dict())
    return binlab_config


def _init_gurobi_cfg(gurobi_cfg_yaml: Path | None) -> pb_lp_cfg.Gurobi:
    gurobi_config = (
        pb_lp_cfg.Gurobi.from_yaml(gurobi_cfg_yaml)
        if gurobi_cfg_yaml is not None
        else pb_lp_cfg.Gurobi.default()
    )
    _LOGGER.debug("Gurobi config:\n%s", gurobi_config.to_dict())
    return gurobi_config


def _init_plasbin_input(
    seed_sequences_tsv: Path,
    sequence_plasmidness_tsv: Path,
) -> tuple[list[str], list[tuple[str, float]]]:
    with pbf_comp_io.SeedReader.open(seed_sequences_tsv) as seeds_fin:
        seeds = list(seeds_fin)
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

    return seeds, _corrected_plasmidness  # BUG BEFORE: plasmidness


def _write_outputs(
    io_manager: pb_io.Manager,
    k: int,
    bin_stats: bins_items.Stats,
    seq_normcovs: Iterable[bins_items.SequenceNormCoverage],
    milp_stats: decomp_lp_views.StatsContainer
    | binlab_lp_views.StatsContainer
    | once_lp_views.MGCLBStats
    | classbin_lp_views.ClassifyStats,
    fn_attributes_from_gurobi_log_path: Callable[[Path], tuple[int, str]]
    | None,  # FIXME tmp fix
    log_files: list[Path],
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

    # FIXME tmp fix
    if fn_attributes_from_gurobi_log_path is not None:
        io_manager.move_gurobi_logs(
            log_files,
            fn_attributes_from_gurobi_log_path,
        )


def _write_mfb_outputs(
    io_manager: mfb_lp_fs.NumberOfFlowManager,
    bin_results: list[
        tuple[
            bins_items.Stats,
            Iterable[bins_items.SequenceNormCoverage],
            mfb_lp_views.Stats,
            Path,
        ]
    ],
) -> None:
    # REFACTOR already done
    io_manager.dir().mkdir(parents=True, exist_ok=True)
    for bin_number, (bin_stats, seq_normcovs, milp_stats, _) in enumerate(
        bin_results,
    ):
        bin_fs_manager = io_manager.bin_fs_manager(bin_number)
        bin_fs_manager.dir().mkdir(parents=True, exist_ok=True)
        bin_stats.to_yaml(bin_fs_manager.bin_stats_path())
        with bins_io.Writer.open(bin_fs_manager.bin_seq_normcov_path()) as fout:
            for seq_normcov in seq_normcovs:
                fout.write_sequence_normcov(
                    seq_normcov.identifier(),
                    seq_normcov.normalized_coverage(),
                )
        milp_stats.to_yaml(bin_fs_manager.milp_stats_path())


def _convert_pangebin_output_to_pbf_output(pangebin_outdir: Path) -> Path:
    pbf_bin_file = pangebin_outdir / "bins.tsv"
    pbf_comp_app.bins(pangebin_outdir, pbf_bin_file)
    return pbf_bin_file


# ------------------------------------------------------------------------------------ #
#                                     Write Configs                                    #
# ------------------------------------------------------------------------------------ #
CONFIGS_APP = typer.Typer(
    name="asm-pbf",
    rich_markup_mode="rich",
    help="Write default configuration files for the asm-pbf pipelines",
)

ARG_CONFIGS_OUTDIR = typer.Argument(help="Output directory")


def _write_binning_config(output_dir: Annotated[Path, ARG_CONFIGS_OUTDIR]) -> None:
    binning_config = pb_cfg.Binning.default()
    binning_cfg_yaml = output_dir / pb_cfg.Binning.DEFAULT_YAML_FILE
    binning_config.to_yaml(binning_cfg_yaml)
    _LOGGER.info("Write general binning config: %s", binning_cfg_yaml)


def _write_decomp_config(output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR]) -> None:
    decomp_config = decomp_cfg.Decomp.default()
    decomp_cfg_yaml = output_directory / decomp_cfg.Decomp.DEFAULT_YAML_FILE
    decomp_config.to_yaml(decomp_cfg_yaml)
    _LOGGER.info("Write general decomp config: %s", decomp_cfg_yaml)


def _write_binlab_config(output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR]) -> None:
    binlab_config = binlab_cfg.Binlab.default()
    binlab_cfg_yaml = output_directory / binlab_cfg.Binlab.DEFAULT_YAML_FILE
    binlab_config.to_yaml(binlab_cfg_yaml)
    _LOGGER.info("Write binlab config: %s", binlab_cfg_yaml)


def _write_gurobi_config(output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR]) -> None:
    gurobi_config = pb_lp_cfg.Gurobi.default()
    gurobi_cfg_yaml = output_directory / pb_lp_cfg.Gurobi.DEFAULT_YAML_FILE
    gurobi_config.to_yaml(gurobi_cfg_yaml)
    _LOGGER.info("Write gurobi config: %s", gurobi_cfg_yaml)


@CONFIGS_APP.command(name="decomp")
def write_decomp_configs(
    output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Write the configuration files for the decomp approach."""
    common_log.init_logger(
        _LOGGER,
        "Write the configuration files for the decomp approach.",
        debug,
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    for write_fn in (_write_binning_config, _write_decomp_config, _write_gurobi_config):
        write_fn(output_directory)


@CONFIGS_APP.command(name="binlab")
def write_binlab_configs(
    output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Write the configuration files for the binlab approach."""
    common_log.init_logger(
        _LOGGER,
        "Write the configuration files for the binlab approach.",
        debug,
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    for write_fn in (_write_binning_config, _write_binlab_config, _write_gurobi_config):
        write_fn(output_directory)


@CONFIGS_APP.command(name="once")
def write_once_configs(
    output_directory: Annotated[Path, ARG_CONFIGS_OUTDIR],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Write the configuration files for the once approach."""
    common_log.init_logger(
        _LOGGER,
        "Write the configuration files for the once approach.",
        debug,
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    for write_fn in (_write_binning_config, _write_gurobi_config):
        write_fn(output_directory)
