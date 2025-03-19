"""Pangebin main application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

from enum import StrEnum

import typer

import pangebin.database.app as db_app
import pangebin.gc_content.app as gc_content_app
import pangebin.gene_density.app as gd_app
import pangebin.gfa.app as gfa_app
import pangebin.ground_truth.app as gt_app
import pangebin.mapping.app as mapping_app
import pangebin.panassembly.app as panassembly_app
import pangebin.pangenome.app as pangenome_app
import pangebin.pipeline.app as pipeline_app
import pangebin.pipeline.seed_thresholds.app as pipe_seed_thr_app
import pangebin.plasbin.app as plasbin_app
import pangebin.seed.app as seed_app
import pangebin.std_asm_graph.app as std_asm_graph_app

APP = typer.Typer(rich_markup_mode="rich")


class _TyperRichHelpPanel(StrEnum):
    """Typer rich help panel categories."""

    PIPELINE = "Pipelines"
    SUBCMD = "Pipeline subcommands"


# ------------------------------------------------------------------------------------ #
#                                   Pipeline Commands                                  #
# ------------------------------------------------------------------------------------ #
APP.command(rich_help_panel=_TyperRichHelpPanel.PIPELINE)(pipeline_app.run)
APP.command(rich_help_panel=_TyperRichHelpPanel.PIPELINE)(
    pipe_seed_thr_app.seed_thresholds,
)

APP.add_typer(pipeline_app.CONFIG_APP, rich_help_panel=_TyperRichHelpPanel.PIPELINE)

# ------------------------------------------------------------------------------------ #
#                                     Sub Commands                                     #
# ------------------------------------------------------------------------------------ #
SUB_APP = typer.Typer(
    name="sub",
    help="Subcommands",
    rich_markup_mode="rich",
)
APP.add_typer(SUB_APP, rich_help_panel=_TyperRichHelpPanel.SUBCMD)


class _SubCMDTyperRichHelpPanel(StrEnum):
    """Typer rich help panel categories."""

    MAIN = "Main commands"
    ATTRIBUTES = "Attribute commands"
    TUNING = "Tuning commands"


SUB_APP.command(rich_help_panel=_SubCMDTyperRichHelpPanel.MAIN)(
    std_asm_graph_app.std_asm_graph,
)
SUB_APP.command(rich_help_panel=_SubCMDTyperRichHelpPanel.MAIN)(pangenome_app.pangenome)
SUB_APP.command(rich_help_panel=_SubCMDTyperRichHelpPanel.MAIN)(
    panassembly_app.panassembly,
)
SUB_APP.command(rich_help_panel=_SubCMDTyperRichHelpPanel.MAIN)(plasbin_app.plasbin)

SUB_APP.add_typer(
    gc_content_app.APP,
    name="gc",
    help="GC content operations.",
    rich_help_panel=_SubCMDTyperRichHelpPanel.ATTRIBUTES,
)
SUB_APP.add_typer(
    gd_app.APP,
    name="gd",
    help="Gene density operations.",
    rich_help_panel=_SubCMDTyperRichHelpPanel.ATTRIBUTES,
)
SUB_APP.add_typer(
    seed_app.APP,
    name="seed",
    help="Seed sequences operations.",
    rich_help_panel=_SubCMDTyperRichHelpPanel.ATTRIBUTES,
)

SUB_APP.command(name="database", rich_help_panel=_SubCMDTyperRichHelpPanel.TUNING)(
    db_app.create,
)
SUB_APP.command(name="ground-truth", rich_help_panel=_SubCMDTyperRichHelpPanel.TUNING)(
    gt_app.create,
)


# ------------------------------------------------------------------------------------ #
#                                   Utility Commands                                   #
# ------------------------------------------------------------------------------------ #
UTILS_APP = typer.Typer(name="utils", help="Utility commands", rich_markup_mode="rich")
APP.add_typer(UTILS_APP, rich_help_panel=_TyperRichHelpPanel.SUBCMD)

UTILS_APP.add_typer(
    gfa_app.APP,
    name="gfa",
    help="GFA operations.",
)
UTILS_APP.add_typer(
    mapping_app.APP,
    name="map",
    help="Mapping operations.",
)


if __name__ == "__main__":
    APP()
