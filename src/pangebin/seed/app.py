"""Seed contig applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gfa.input_output as gfa_io
import pangebin.logging as common_log
import pangebin.seed.create as seed_create
import pangebin.seed.input_output as seed_io
import pangebin.seed.thresholds.app as seed_thr_app

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")

APP.add_typer(seed_thr_app.APP)


class FromGeneDensityArguments:
    """From gene density arguments."""

    GENE_DENSITY_FILE = typer.Argument(
        help="Gene density file",
    )

    OUTPUT_SEED_TSV = typer.Argument(
        help="Seed sequence output file",
    )


@APP.command(name="pos-gd")
def from_positive_gene_densities(
    gene_density_file: Annotated[Path, FromGeneDensityArguments.GENE_DENSITY_FILE],
    output_seed_tsv: Annotated[Path, FromGeneDensityArguments.OUTPUT_SEED_TSV],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Extract seed sequences with positive gene density."""
    common_log.init_logger(
        _LOGGER,
        "Extracting seed sequences with positive gene density.",
        debug,
    )
    with seed_io.Writer.open(output_seed_tsv) as writer:
        for sequence_id in seed_create.from_gene_density(gene_density_file):
            writer.write_sequence(sequence_id)
    _LOGGER.info("Write seed sequence identifiers in file: %s", output_seed_tsv)
    return output_seed_tsv


class FromContigsToFragmentSeedsArguments:
    """From contigs to fragment seeds arguments."""

    SKESA_SEEDS_TSV = typer.Argument(
        help="SKESA seeds TSV file",
    )

    UNICYCLER_SEEDS_TSV = typer.Argument(
        help="Unicycler seeds TSV file",
    )

    PANASSEMBLY_GFA = typer.Argument(
        help="Pan-assembly GFA file",
    )

    OUTPUT_SEED_TSV = typer.Argument(
        help="Seed sequence output file",
    )


@APP.command(name="ctgs-to-frags")
def from_contigs_to_fragment_seeds(
    skesa_seeds_tsv: Annotated[
        Path,
        FromContigsToFragmentSeedsArguments.SKESA_SEEDS_TSV,
    ],
    unicycler_seeds_tsv: Annotated[
        Path,
        FromContigsToFragmentSeedsArguments.UNICYCLER_SEEDS_TSV,
    ],
    panassembly_gfa: Annotated[
        Path,
        FromContigsToFragmentSeedsArguments.PANASSEMBLY_GFA,
    ],
    output_seed_tsv: Annotated[
        Path,
        FromContigsToFragmentSeedsArguments.OUTPUT_SEED_TSV,
    ],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> Path:
    """Extract seed fragments from seed contigs."""
    common_log.init_logger(
        _LOGGER,
        "Extracting seed fragments from seed contigs.",
        debug,
    )
    with (
        seed_io.Reader.open(skesa_seeds_tsv) as skesa_seeds_fin,
        seed_io.Reader.open(unicycler_seeds_tsv) as uni_seeds_fin,
        seed_io.Writer.open(output_seed_tsv) as writer,
    ):
        for seed_fragment in seed_create.from_contigs_to_fragment_seeds(
            iter(skesa_seeds_fin),
            iter(uni_seeds_fin),
            gfa_io.from_file(panassembly_gfa),
        ):
            writer.write_sequence(seed_fragment)
    _LOGGER.info("Write fragment seeds in file: %s", output_seed_tsv)
    return output_seed_tsv
