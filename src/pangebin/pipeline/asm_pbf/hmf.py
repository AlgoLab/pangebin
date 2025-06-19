"""HMF PlasBin-flow compatibility application wrapper module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.pbf_comp.input_output as pbf_io
import pangebin.pbf_comp.items as pbf_items
import pangebin.plasbin.hmf.config as hmf_cfg
import pangebin.plasbin.hmf.create as hmf_create
import pangebin.plasbin.hmf.results as hmf_res
import pangebin.plasbin.input_output as pb_io
import pangebin.plasbin.milp.config as milp_cfg
from pangebin import pblog

from . import common as cmn

APP = typer.Typer(rich_markup_mode="rich")

_LOGGER = logging.getLogger(__name__)


class HMFOptions:
    """HMF options."""

    _RICH_HELP_PANEL = "Pangebin HMF options"

    BINNING_CFG_YAML = typer.Option(
        "--bin-cfg",
        help="The general binning configuration YAML path",
        rich_help_panel=_RICH_HELP_PANEL,
    )


@APP.command(help="Run Pangebin HMF approach.")
def hmf(
    assembly_gfa: Annotated[Path, cmn.Arguments.ASSEMBLY_GFA],
    seed_contigs_tsv: Annotated[Path, cmn.Arguments.SEED_CONTIGS_TSV],
    contig_plasmidness_tsv: Annotated[Path, cmn.Arguments.CONTIG_PLASMIDNESS_TSV],
    # pbf options
    is_skesa: Annotated[bool, cmn.InputOptions.SKESA_GFA] = False,
    gc_scores_tsv: Annotated[Path | None, cmn.InputOptions.GC_SCORES_TSV] = None,
    # Binning options
    binning_cfg_yaml: Annotated[Path | None, HMFOptions.BINNING_CFG_YAML] = None,
    # Gurobi options
    gurobi_cfg_yaml: Annotated[Path | None, cmn.GurobiOptions.GUROBI_CFG_YAML] = None,
    # IO options
    outdir: Annotated[Path, cmn.IOOptions.OUTPUT_DIR] = pb_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> Path:
    """Run Pangebin HMF approach.

    Returns
    -------
    Path
        Path to the PlasBin-flow bin file.
    """
    pblog.init_logger(_LOGGER, "Running pangebin classbin approach.", debug)
    outdir.mkdir(parents=True, exist_ok=True)

    std_gfa = cmn.prepare_graph(assembly_gfa, is_skesa, debug)

    _, gc_scores = cmn.prepare_gc_scores_tsv(assembly_gfa, gc_scores_tsv)

    hmf_config = _init_binning_cfg(binning_cfg_yaml)
    gurobi_config = cmn.init_gurobi_cfg(gurobi_cfg_yaml)

    seeds = cmn.init_seeds(seed_contigs_tsv)
    plasmidness = cmn.init_plasmidness(contig_plasmidness_tsv)

    # XXX PAER TESTS
    # REFACTOR clarify this logic
    _mean_plasmidness_with_seeds(seeds, plasmidness)

    best_instances_reader = hmf_create.plasbin_assembly(
        std_gfa,
        seeds,
        gc_scores,
        plasmidness,
        hmf_config,
        gurobi_config,
        outdir,
    )

    return _convert_output_to_pbf_output(best_instances_reader)


def _mean_plasmidness_with_seeds(
    seeds: list[str],
    plasmidness: list[tuple[str, float]],
) -> None:
    seeds_set = set(seeds)
    for k, (frag_id, plm) in enumerate(plasmidness):
        if frag_id in seeds_set:
            plasmidness[k] = (frag_id, (plm + 1) / 2)


def _init_binning_cfg(config_yaml: Path | None) -> hmf_cfg.Config:
    config = (
        hmf_cfg.Config.from_yaml(config_yaml)
        if config_yaml is not None
        else hmf_cfg.Config.default()
    )
    _LOGGER.debug("HMF configurations:\n%s", config.to_dict())
    return config


def _convert_output_to_pbf_output(best_instances_reader: hmf_res.RootReader) -> Path:
    _LOGGER.info("Converting PangeBin-flow bins into PlasBin-flow bin info TSV file.")
    pbf_bin_file = best_instances_reader.file_system().dir() / "bins.tsv"
    with pbf_io.BinsWriter.open(pbf_bin_file) as pbf_bin_info_fout:
        for bin_number, bin_reader_with_info in enumerate(
            best_instances_reader.all_bins(),
        ):
            seq_mults = [
                pbf_items.ContigMult(
                    snc.identifier(),
                    snc.normalized_coverage(),
                )
                for snc in bin_reader_with_info.results_reader().iter_seq_normcov()
            ]
            bin_stats = bin_reader_with_info.results_reader().bin_stats()
            pbf_bin_info_fout.write_bin_line(
                pbf_items.PBFBinInfo(
                    f"P{bin_number + 1}",
                    bin_stats.total_flow(),
                    (0.0, 1.0),
                    seq_mults,
                ),
            )

    _LOGGER.info("PlasBin-flow bin info TSV file: %s", pbf_bin_file)
    return pbf_bin_file


def write_default_configs(
    binning_cfg_yaml: Annotated[Path | None, HMFOptions.BINNING_CFG_YAML] = None,
    # Gurobi options
    gurobi_cfg_yaml: Annotated[Path | None, cmn.GurobiOptions.GUROBI_CFG_YAML] = None,
) -> None:
    """Write default HMF configurations."""
    _LOGGER.info("Write default HMF configurations.")
    if binning_cfg_yaml is None:
        binning_cfg_yaml = Path("./hmf_config.yaml")

    _LOGGER.info("Write HMF configurations: %s", binning_cfg_yaml)
    hmf_cfg.Config.default().to_yaml(binning_cfg_yaml)

    if gurobi_cfg_yaml is None:
        gurobi_cfg_yaml = Path("./gurobi_config.yaml")

    _LOGGER.info("Write Gurobi configurations: %s", gurobi_cfg_yaml)
    milp_cfg.Gurobi.default().to_yaml(gurobi_cfg_yaml)
