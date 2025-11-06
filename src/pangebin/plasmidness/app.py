"""Seed contig applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.gfa.input_output as gfa_io
import pangebin.plasmidness.create as plm_create
import pangebin.plasmidness.input_output as plm_io
from pangebin import pblog

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class FromContigsToFragmentPlasmidnessArguments:
    """From contigs to fragment plasmidness arguments."""

    SKESA_PLASMIDNESS_TSV = typer.Argument(
        help="SKESA plasmidness TSV file",
    )

    UNICYCLER_PLASMIDNESS_TSV = typer.Argument(
        help="Unicycler plasmidness TSV file",
    )

    PANASSEMBLY_GFA = typer.Argument(
        help="Pan-assembly GFA file",
    )

    OUTPUT_PLASMIDNESS_TSV = typer.Argument(
        help="Fragments plasmidness TSV output file",
    )


@APP.command(name="ctgs-to-frags")
def from_contigs_to_fragment_plasmidness(
    skesa_plasmidness_tsv: Annotated[
        Path,
        FromContigsToFragmentPlasmidnessArguments.SKESA_PLASMIDNESS_TSV,
    ],
    unicycler_plasmidness_tsv: Annotated[
        Path,
        FromContigsToFragmentPlasmidnessArguments.UNICYCLER_PLASMIDNESS_TSV,
    ],
    panassembly_gfa: Annotated[
        Path,
        FromContigsToFragmentPlasmidnessArguments.PANASSEMBLY_GFA,
    ],
    output_plasmidness_tsv: Annotated[
        Path,
        FromContigsToFragmentPlasmidnessArguments.OUTPUT_PLASMIDNESS_TSV,
    ],
    debug: Annotated[bool, pblog.OPT_DEBUG] = False,
) -> Path:
    """Get fragments plasmidness from contigs plasmidness."""
    pblog.init_logger(
        _LOGGER,
        "Converting contigs plasmidness to fragments plasmidness.",
        debug,
    )
    with (
        plm_io.Reader.open(skesa_plasmidness_tsv) as skesa_plasmidness_fin,
        plm_io.Reader.open(unicycler_plasmidness_tsv) as uni_plasmidness_fin,
        plm_io.Writer.open(output_plasmidness_tsv) as writer,
    ):
        for (
            frag_id,
            frag_plasmidness,
        ) in plm_create.from_contigs_to_fragments_plasmidness(
            iter(skesa_plasmidness_fin),
            iter(uni_plasmidness_fin),
            gfa_io.from_file(panassembly_gfa),
        ):
            writer.write_sequence_plasmidness(frag_id, frag_plasmidness)
    _LOGGER.info("Write fragments plasmidness in file: %s", output_plasmidness_tsv)
    return output_plasmidness_tsv


# TODO from gene density
