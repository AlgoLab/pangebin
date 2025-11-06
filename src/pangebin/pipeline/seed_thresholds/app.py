"""Pangebin main pipeline application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
import shutil
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path
from typing import Annotated

import typer

import pangebin.assembly.create as asm_create
import pangebin.database.input_output as db_io
import pangebin.gene_density.app as gd_app
import pangebin.ground_truth.app as gt_app
import pangebin.mapping.app as mapping_app
import pangebin.pipeline.seed_thresholds.config as pipe_seed_thr_cfg
import pangebin.seed.thresholds.app as seed_thr_app
import pangebin.seed.thresholds.input_output as seed_thr_io
import pangebin.seed.thresholds.items as seed_thr_items
import pangebin.sra_tools as sra
from pangebin import pblog

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class Arguments:
    """Seed thresholds arguments."""

    DATATEST_YAML = typer.Argument(
        help="Paired Illumina datatest YAML file",
    )

    GENE_FASTA = typer.Argument(
        help="Genes FASTA file",
    )


class Options:
    """Seed thresholds options."""

    __RICH_HELP_PANEL = "Configurations"

    DEFAULT_CONFIG_FILE = (
        Path(__file__).parent / pipe_seed_thr_cfg.Config.DEFAULT_FILENAME
    )

    CONFIG_FILE = typer.Option(
        "--config",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )


class IOOptions:
    """Seed thresholds I/O options."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTDIR = typer.Option(
        help="Output directory",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def seed_thresholds(
    datatest_yaml: Annotated[Path, Arguments.DATATEST_YAML],
    gene_fasta: Annotated[
        Path,
        Arguments.GENE_FASTA,
    ],
    # Configurations
    config_file: Annotated[
        Path,
        Options.CONFIG_FILE,
    ] = Options.DEFAULT_CONFIG_FILE,
    # IO options
    outdir: Annotated[
        Path,
        IOOptions.OUTDIR,
    ] = seed_thr_io.Config.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> seed_thr_io.Manager:
    """Obtain the seed threshold pairs from paired Illumina BioSamples."""
    pblog.init_logger(
        _LOGGER,
        "Obtaining the seed threshold pairs from paired Illumina BioSamples.",
        debug,
    )
    io_manager = seed_thr_io.Manager(seed_thr_io.Config(outdir))
    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)

    config = pipe_seed_thr_cfg.Config.from_yaml(config_file)

    seed_ctg_thr_dataset_tsv = (
        io_manager.config().output_directory() / "seed_ctg_thr_dataset.tsv"
    )
    #
    # Temporary configuration files
    #
    with __tmp_seed_thresholds_config_files(config, io_manager) as (
        ground_truth_cfg_file,
        entrez_cfg_file,
        gene_on_contigs_sam_filter_cfg_file,
        thr_ranges_cfg_file,
    ):
        for datatest in db_io.illumina_biosamples_from_file(datatest_yaml):
            first_srr_id = datatest.run_sra_ids()[0]
            datatest_outdir = io_manager.config().output_directory() / first_srr_id
            datatest_outdir.mkdir(parents=True, exist_ok=True)
            fastq_1, fastq_2 = sra.download_sra_fastq_paired(
                first_srr_id,
                datatest_outdir,
                remove_dir=True,
            )
            assembly_outdir = datatest_outdir / "assembly"
            assembly_outdir.mkdir(parents=True, exist_ok=True)
            asm_create.unicycler_paired_end_reads(
                fastq_1,
                fastq_2,
                assembly_outdir,
                config.ressources_config(),
            )
            # remove short reads
            fastq_1.unlink()
            fastq_2.unlink()

            asm_fasta = assembly_outdir / "assembly.fasta"
            gt_io_manager = gt_app.create(
                asm_fasta,
                datatest.plasmids(),
                config_file=ground_truth_cfg_file,
                entrez_config_file=entrez_cfg_file,
                output_dir=datatest_outdir,
                debug=debug,
            )

            gene_mapping_on_contigs_sam = mapping_app.blast(
                gene_fasta,
                asm_fasta,
                datatest_outdir / "gene_mapping_on_contigs.sam",
                debug=debug,
            )
            filtered_gene_mapping_on_contigs_sam = mapping_app.filter(
                gene_mapping_on_contigs_sam,
                query_fasta=gene_fasta,
                config_file=gene_on_contigs_sam_filter_cfg_file,
                debug=debug,
            )
            gene_mapping_on_contigs_sam.unlink()

            ctg_gd_tsv = gd_app.fasta(
                asm_fasta,
                filtered_gene_mapping_on_contigs_sam,
                datatest_outdir / "ctg_gene_density.tsv",
                debug=debug,
            )

            seed_ctg_thr_test_item = seed_thr_items.TestItem(
                gt_io_manager.plasmid_contigs_file(),
                gt_io_manager.non_plasmid_contigs_file(),
                ctg_gd_tsv,
            )
            with seed_ctg_thr_dataset_tsv.open("a") as tsv_out:
                tsv_out.write(f"{seed_ctg_thr_test_item.to_dataset_line()}\n")

    return seed_thr_app.thresholds(
        seed_ctg_thr_dataset_tsv,
        config_file=thr_ranges_cfg_file,
        output_dir=io_manager.config().output_directory(),
        debug=debug,
    )


@contextmanager
def __tmp_seed_thresholds_config_files(
    config: pipe_seed_thr_cfg.Config,
    io_manager: seed_thr_io.Manager,
) -> Generator[
    tuple[Path, Path, Path, Path],
    None,
    None,
]:
    ground_truth_cfg_file = config.ground_truth_config().to_yaml(
        io_manager.config().output_directory() / "ground_truth_cfg.yaml",
    )
    entrez_cfg_file = config.entrez_config().to_yaml(
        io_manager.config().output_directory() / "entrez_cfg.yaml",
    )
    gene_on_contigs_sam_filter_cfg_file = config.sam_filter_config().to_yaml(
        io_manager.config().output_directory() / "gene_on_ctgs_sam_filter_cfg.yaml",
    )
    thr_ranges_cfg_file = config.threshold_ranges().to_yaml(
        io_manager.config().output_directory() / "threshold_ranges_cfg.yaml",
    )

    yield (
        ground_truth_cfg_file,
        entrez_cfg_file,
        gene_on_contigs_sam_filter_cfg_file,
        thr_ranges_cfg_file,
    )

    ground_truth_cfg_file.unlink()
    entrez_cfg_file.unlink()
    gene_on_contigs_sam_filter_cfg_file.unlink()
    thr_ranges_cfg_file.unlink()


ARG_CONFIG_PATH = typer.Argument(help="Output directory")


def write_configs(
    config_path: Annotated[
        Path,
        ARG_CONFIG_PATH,
    ] = pipe_seed_thr_cfg.Config.DEFAULT_FILENAME,
) -> None:
    """Write the configuration file for the seed thresholds pipeline."""
    shutil.copyfile(
        Options.DEFAULT_CONFIG_FILE,
        config_path,
    )
