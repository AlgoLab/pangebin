"""PlasBin-flow compatibility application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.logging as common_log
import pangebin.pbf_comp.input_output as comp_io
import pangebin.plasmidness.input_output as plm_io
import pangebin.seed.input_output as seed_io

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class PlasmidnessArguments:
    """Convert plasmidness file arguments."""

    PBF_PLASMIDNESS_FILE = typer.Argument(
        help="PlasBin-flow plasmidness file",
    )

    PG_PLASMIDNESS_TSV = typer.Argument(
        help="Pangebin plasmidness TSV file",
    )


@APP.command("plm")
def plasmidness(
    pbf_plasmidness_file: Annotated[
        Path,
        PlasmidnessArguments.PBF_PLASMIDNESS_FILE,
    ],
    pg_plasmidness_tsv: Annotated[
        Path,
        PlasmidnessArguments.PG_PLASMIDNESS_TSV,
    ],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Convert PlasBin-flow plasmidness file to PangeBin plasmidness TSV file."""
    common_log.init_logger(_LOGGER, "Converting PlasBin-flow plasmidness file.", debug)
    with (
        comp_io.PBFPLMReader.open(pbf_plasmidness_file) as pbf_plasmidness_fin,
        plm_io.Writer.open(
            pg_plasmidness_tsv,
        ) as pg_plasmidness_fout,
    ):
        for sequence_id, plasmidness in pbf_plasmidness_fin:
            pg_plasmidness_fout.write_sequence_plasmidness(
                sequence_id,
                2 * plasmidness - 1,
            )
    _LOGGER.info(
        "PlasBin-flow plasmidness file %s converted to Pangebin plasmidness file %s",
        pbf_plasmidness_file,
        pg_plasmidness_tsv,
    )


class SeedArguments:
    """Convert seed sequences file arguments."""

    PBF_SEEDS_FILE = typer.Argument(
        help="PlasBin-flow seed sequences file",
    )

    PG_SEEDS_TSV = typer.Argument(
        help="Pangebin seed sequences TSV file",
    )


@APP.command()
def seed(
    pbf_seeds_file: Annotated[
        Path,
        SeedArguments.PBF_SEEDS_FILE,
    ],
    pg_seeds_tsv: Annotated[
        Path,
        SeedArguments.PG_SEEDS_TSV,
    ],
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Convert PlasBin-flow seed sequences file to PangeBin seed sequences TSV file."""
    common_log.init_logger(
        _LOGGER,
        "Converting PlasBin-flow seed sequences file.",
        debug,
    )
    with (
        comp_io.PBFSeedReader.open(pbf_seeds_file) as pbf_seeds_fin,
        seed_io.Writer.open(
            pg_seeds_tsv,
        ) as pg_seeds_fout,
    ):
        for sequence_id in pbf_seeds_fin:
            pg_seeds_fout.write_sequence(sequence_id)
    _LOGGER.info(
        "PlasBin-flow seed sequences file %s"
        " converted to Pangebin seed sequences file %s",
        pbf_seeds_file,
        pg_seeds_tsv,
    )
