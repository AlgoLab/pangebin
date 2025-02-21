"""Pangebin main pipeline application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Annotated

import typer

import pangebin.assembly.create as asm_create
import pangebin.database.input_output as db_io
import pangebin.entrez as pg_entrez
import pangebin.gene_density.app as gd_app
import pangebin.ground_truth.app as gt_app
import pangebin.logging as common_log
import pangebin.mapping.app as mapping_app
import pangebin.seed.app as seed_app
import pangebin.seed.input_output as seed_io
import pangebin.seed.items as seed_items
import pangebin.sra_tools as sra

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class SeedThresholdsArguments:
    """Seed threshold arguments."""

    DATATEST_YAML = typer.Argument(
        help="Paired Illumina datatest YAML file",
    )

    GENE_FASTA = typer.Argument(
        help="Genes FASTA file",
    )


class SeedThresholdsIOOptions:
    """Seed threshold I/O options."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTDIR = typer.Option(
        help="Output directory",
        rich_help_panel=__RICH_HELP_PANEL,
    )


class SeedThresholdsOptions:
    """Seed threshold options."""

    __RICH_HELP_PANEL = "Configurations"

    GROUND_TRUTH_CFG_FILE = typer.Option(
        "--gd-cfg",
        help="Ground truth config file",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    FILTER_GENE_ON_CONTIGS_SAM_CFG_FILE = typer.Option(
        "--sam-filter-cfg",
        help="Filter gene on contigs SAM config file",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    THR_RANGES_CFG_FILE = typer.Option(
        "--thr-ranges-cfg",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )


DEFAULT_GD_CFG_FILE = Path(__file__).parent / "ground_truth_config.yaml"
DEFAULT_SAM_FILTER_CFG_FILE = Path(__file__).parent / "sam_filter_config.yaml"
DEFAULT_THR_RANGES_CFG_FILE = Path(__file__).parent / "thr_ranges_config.yaml"


@APP.command()
def seed_thresholds(
    datatest_yaml: Annotated[Path, SeedThresholdsArguments.DATATEST_YAML],
    gene_fasta: Annotated[
        Path,
        SeedThresholdsArguments.GENE_FASTA,
    ],
    # Configurations
    entrez_cfg_file: Annotated[
        Path | None,
        pg_entrez.AppOptions.CONFIG_FILE,
    ] = None,
    ground_truth_cfg_file: Annotated[
        Path,
        SeedThresholdsOptions.GROUND_TRUTH_CFG_FILE,
    ] = DEFAULT_GD_CFG_FILE,
    filter_gene_on_contigs_sam_cfg_file: Annotated[
        Path,
        SeedThresholdsOptions.FILTER_GENE_ON_CONTIGS_SAM_CFG_FILE,
    ] = DEFAULT_SAM_FILTER_CFG_FILE,
    thr_ranges_cfg_file: Annotated[
        Path,
        SeedThresholdsOptions.THR_RANGES_CFG_FILE,
    ] = DEFAULT_THR_RANGES_CFG_FILE,
    # IO options
    outdir: Annotated[
        Path,
        SeedThresholdsIOOptions.OUTDIR,
    ] = seed_io.ThresholdsConfig.DEFAULT_OUTPUT_DIR,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Obtain the seed threshold pairs from paired Illumina BioSamples."""
    common_log.init_logger(
        _LOGGER,
        "Obtaining the seed threshold pairs from paired Illumina BioSamples.",
        debug,
    )
    io_manager = seed_io.ThresholdsManager(seed_io.ThresholdsConfig(outdir))
    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)

    seed_ctg_thr_dataset_tsv = (
        io_manager.config().output_directory() / "seed_ctg_thr_dataset.tsv"
    )

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
        )
        # remove short reads
        fastq_1.unlink()
        fastq_2.unlink()

        asm_fasta = assembly_outdir / "assembly.fasta"
        gt_app.create(
            asm_fasta,
            datatest.plasmids(),
            config_file=ground_truth_cfg_file,
            entrez_config_file=entrez_cfg_file,
            output_dir=datatest_outdir,
            debug=debug,
        )

        gene_mapping_on_contigs_sam = datatest_outdir / "gene_mapping_on_contigs.sam"
        mapping_app.blast(
            gene_fasta,
            asm_fasta,
            gene_mapping_on_contigs_sam,
            debug=debug,
        )
        filtered_gene_mapping_on_contigs_sam = (
            datatest_outdir / "gene_mapping_on_contigs.filtered.sam"
        )
        mapping_app.filter(
            gene_mapping_on_contigs_sam,
            filtered_gene_mapping_on_contigs_sam,
            query_fasta=gene_fasta,
            config_file=filter_gene_on_contigs_sam_cfg_file,
            debug=debug,
        )
        gene_mapping_on_contigs_sam.unlink()

        ctg_gd_tsv = datatest_outdir / "ctg_gene_density.tsv"
        gd_app.fasta(
            asm_fasta,
            filtered_gene_mapping_on_contigs_sam,
            ctg_gd_tsv,
            debug=debug,
        )

        seed_ctg_thr_test_item = seed_items.SeedContigThresholdTestItem(
            datatest_outdir / "plasmid_contigs.tsv",
            datatest_outdir / "non_plasmid_contigs.tsv",
            ctg_gd_tsv,
        )
        with seed_ctg_thr_dataset_tsv.open("a") as tsv_out:
            tsv_out.write(f"{seed_ctg_thr_test_item.to_dataset_line()}\n")

    seed_app.thresholds(
        seed_ctg_thr_dataset_tsv,
        config_file=thr_ranges_cfg_file,
        output_dir=io_manager.config().output_directory(),
        debug=debug,
    )


PRINT_CFG_OUTDIR_OPT = typer.Argument(help="Output directory")
DEFAULT_PRINT_CFG_OUTDIR = Path()


def write_configs(
    outdir: Annotated[Path, PRINT_CFG_OUTDIR_OPT] = DEFAULT_GD_CFG_FILE,
) -> None:
    """Write the configuration files for seed thresholds pipeline."""
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(DEFAULT_GD_CFG_FILE, outdir / DEFAULT_GD_CFG_FILE.name)
    shutil.copyfile(
        DEFAULT_SAM_FILTER_CFG_FILE,
        outdir / DEFAULT_SAM_FILTER_CFG_FILE.name,
    )
    shutil.copyfile(
        DEFAULT_THR_RANGES_CFG_FILE,
        outdir / DEFAULT_THR_RANGES_CFG_FILE.name,
    )
