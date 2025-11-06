"""PlasBin-flow compatibility application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.pbf_comp.input_output as pbf_io
import pangebin.pbf_comp.ops as comp_ops
import pangebin.plasbin.input_output as pb_io
import pangebin.plasmidness.input_output as plm_io
import pangebin.seed.input_output as seed_io
from pangebin import pblog

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
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> None:
    """Convert PlasBin-flow plasmidness file to PangeBin plasmidness TSV file."""
    pblog.init_logger(_LOGGER, "Converting PlasBin-flow plasmidness file.", debug)
    with (
        pbf_io.PlmReader.open(pbf_plasmidness_file) as pbf_plasmidness_fin,
        plm_io.Writer.open(
            pg_plasmidness_tsv,
        ) as pg_plasmidness_fout,
    ):
        for sequence_id, plasmidness in pbf_plasmidness_fin:
            pg_plasmidness_fout.write_sequence_plasmidness(
                sequence_id,
                comp_ops.pbf_to_pg_plasmidness(plasmidness),
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
def seeds(
    pbf_seeds_file: Annotated[
        Path,
        SeedArguments.PBF_SEEDS_FILE,
    ],
    pg_seeds_tsv: Annotated[
        Path,
        SeedArguments.PG_SEEDS_TSV,
    ],
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> None:
    """Convert PlasBin-flow seed sequences file to PangeBin seed sequences TSV file."""
    pblog.init_logger(
        _LOGGER,
        "Converting PlasBin-flow seed sequences file.",
        debug,
    )
    with (
        pbf_io.SeedReader.open(pbf_seeds_file) as pbf_seeds_fin,
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


class BinArguments:
    """Convert bins file arguments."""

    PG_BINS_PARENT_DIR = typer.Argument(
        help="Pangebin bins parent directory",
    )

    PBF_BIN_INFO_TSV = typer.Argument(
        help="PlasBin-flow bin info TSV file",
    )


@APP.command()
def bins(
    pg_bins_parent_dir: Annotated[
        Path,
        BinArguments.PG_BINS_PARENT_DIR,
    ],
    pbf_bin_info_tsv: Annotated[
        Path,
        BinArguments.PBF_BIN_INFO_TSV,
    ],
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> None:
    """Convert PangeBin-flow bins into PlasBin-flow bin info TSV file."""
    pblog.init_logger(
        _LOGGER,
        "Converting PangeBin-flow bins into PlasBin-flow bin info TSV file.",
        debug,
    )
    io_manager = pb_io.Manager(pb_io.Config(output_directory=pg_bins_parent_dir))
    with pbf_io.BinsWriter.open(pbf_bin_info_tsv) as pbf_bin_info_fout:
        for k in range(io_manager.number_of_bins()):
            pbf_bin_info_fout.write_bin_line(
                comp_ops.pg_bin_dir_to_pbf_bininfo(io_manager, k),
            )
    _LOGGER.info("PlasBin-flow bin info TSV file: %s", pbf_bin_info_tsv)
