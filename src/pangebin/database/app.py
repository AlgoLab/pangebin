"""Database applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

import pangebin.database.create as db_create
import pangebin.database.input_output as db_io
import pangebin.entrez as pg_entrez
import pangebin.logging as common_log

_LOGGER = logging.getLogger(__name__)

APP = typer.Typer(rich_markup_mode="rich")


class CreateArguments:
    """Arguments for create database application."""

    ACCESSIONS_FILE = typer.Argument(
        help="File with Entrez accessions",
    )


class CreateIOOptions:
    """Input-output options for create database application."""

    __RICH_HELP_PANEL = "Input/Output options"

    OUTPUT_DIR = typer.Option(
        "-o",
        "--outdir",
        help="Output directory",
        rich_help_panel=__RICH_HELP_PANEL,
    )


@APP.command()
def create(
    accessions_file: Annotated[
        Path,
        CreateArguments.ACCESSIONS_FILE,
    ],
    output_dir: Annotated[
        Path,
        CreateIOOptions.OUTPUT_DIR,
    ] = db_io.Config.DEFAULT_OUTPUT_DIR,
    # Entrez configuration
    entrez_email: Annotated[
        str | None,
        pg_entrez.AppOptions.EMAIL,
    ] = pg_entrez.Config.DEFAULT_EMAIL,
    entrez_tool: Annotated[
        str,
        pg_entrez.AppOptions.TOOL,
    ] = pg_entrez.Config.DEFAULT_TOOL,
    entrez_api_key: Annotated[
        str | None,
        pg_entrez.AppOptions.API_KEY,
    ] = pg_entrez.Config.DEFAULT_API_KEY,
    entrez_max_tries: Annotated[
        int,
        pg_entrez.AppOptions.MAX_TRIES,
    ] = pg_entrez.Config.DEFAULT_MAX_TRIES,
    entrez_sleep_between_tries: Annotated[
        int,
        pg_entrez.AppOptions.SLEEP_BETWEEN_TRIES,
    ] = pg_entrez.Config.DEFAULT_SLEEP_BETWEEN_TRIES,
    entrez_config_file: Annotated[
        Path | None,
        pg_entrez.AppOptions.CONFIG_FILE,
    ] = None,
    debug: Annotated[bool, common_log.OPT_DEBUG] = False,
) -> None:
    """Create plasmid database (biosamples, plasmids and paired Illumina SRA ids)."""
    common_log.init_logger(_LOGGER, "Creating database.", debug)

    io_manager = db_io.Manager(db_io.Config(output_directory=output_dir))
    io_manager.config().output_directory().mkdir(parents=True, exist_ok=True)

    config = (
        pg_entrez.Config.from_yaml(entrez_config_file)
        if entrez_config_file
        else pg_entrez.Config(
            email=entrez_email,
            tool=entrez_tool,
            api_key=entrez_api_key,
            max_tries=entrez_max_tries,
            sleep_between_tries=entrez_sleep_between_tries,
        )
    )

    illumina_biosamples, non_illumina_biosamples = (
        db_create.illumina_exp_sra_ids_from_plasmid_accessions(
            db_io.accessions_from_file(accessions_file),
            config,
        )
    )

    db_io.illumina_biosamples_to_file(
        illumina_biosamples,
        io_manager.illumina_biosamples_yaml(),
    )
    db_io.non_illumina_biosamples_to_file(
        non_illumina_biosamples,
        io_manager.non_illumina_biosamples_yaml(),
    )
