"""Gene density application."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gene_density.create as gd_create
import pangebin.gene_density.input_output as gd_io
import pangebin.logging as common_log

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class SequenceGeneDensityArguments:
    """Sequence gene density arguments."""

    FASTA_FILE = typer.Argument(
        help="FASTA file",
    )

    GENE_MAPPING_SAM = typer.Argument(
        help="Gene mapping SAM file",
    )

    OUTPUT_FILE = typer.Argument(
        help="Output file",
    )


@APP.command()
def fasta(
    fasta_file: Annotated[Path, SequenceGeneDensityArguments.FASTA_FILE],
    gene_mapping_sam: Annotated[Path, SequenceGeneDensityArguments.GENE_MAPPING_SAM],
    output_file: Annotated[Path, SequenceGeneDensityArguments.OUTPUT_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Compute the gene densities from the mapping of genes against the sequences."""
    common_log.init_logger(_LOGGER, "Computing sequence gene densities.", debug)
    sequence_gene_densities, sequence_intervals = gd_create.sequence_gene_density(
        fasta_file,
        gene_mapping_sam,
    )
    gd_io.to_file_with_intervals(
        output_file,
        sequence_gene_densities,
        gene_mapping_intervals=sequence_intervals,
    )


class FragmentGeneDensityArguments:
    """Fragment gene density arguments."""

    PANASSEMBLY_GFA = typer.Argument(
        help="Pan-assembly GFA file",
    )

    GENE_MAPPING_TO_CONTIGS_SAM = typer.Argument(
        help="Gene mapping to contigs SAM file",
    )

    OUTPUT_FILE = typer.Argument(
        help="Output file",
    )


@APP.command()
def frag(
    panassembly_gfa: Annotated[Path, FragmentGeneDensityArguments.PANASSEMBLY_GFA],
    gene_mapping_to_contigs_sam: Annotated[
        Path,
        FragmentGeneDensityArguments.GENE_MAPPING_TO_CONTIGS_SAM,
    ],
    output_file: Annotated[Path, FragmentGeneDensityArguments.OUTPUT_FILE],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Compute fragment gene densities from the mapping of genes against the contigs."""
    common_log.init_logger(_LOGGER, "Computing fragment gene densities.", debug)
    fragment_gene_densities, fragment_intervals = gd_create.fragment_gene_density(
        panassembly_gfa,
        gene_mapping_to_contigs_sam,
    )
    gd_io.to_file_with_intervals(
        output_file,
        fragment_gene_densities,
        gene_mapping_intervals=fragment_intervals,
    )
