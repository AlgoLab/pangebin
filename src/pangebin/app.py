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
import pangebin.seed.app as seed_app
import pangebin.std_asm_graph.app as std_asm_graph_app


class _TyperRichHelpPanel(StrEnum):
    """Typer rich help panel categories."""

    MAIN = "Main commands"
    TUNING = "Tuning commands"
    ATTRIBUTES = "Attribute commands"
    UTILS = "Utility commands"


APP = typer.Typer(rich_markup_mode="rich")

APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(std_asm_graph_app.std_asm_graph)
APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(pangenome_app.pangenome)
APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(panassembly_app.panassembly)
# APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(plasbin_app.plasbin)

APP.command(name="database", rich_help_panel=_TyperRichHelpPanel.TUNING)(
    db_app.create,
)
APP.command(name="ground-truth", rich_help_panel=_TyperRichHelpPanel.TUNING)(
    gt_app.create,
)

APP.add_typer(
    gfa_app.APP,
    name="gfa",
    help="GFA operations.",
    rich_help_panel=_TyperRichHelpPanel.UTILS,
)
APP.add_typer(
    mapping_app.APP,
    name="map",
    help="Mapping operations.",
    rich_help_panel=_TyperRichHelpPanel.UTILS,
)

APP.add_typer(
    gc_content_app.APP,
    name="gc",
    help="GC content operations.",
    rich_help_panel=_TyperRichHelpPanel.ATTRIBUTES,
)
APP.add_typer(
    gd_app.APP,
    name="gd",
    help="Gene density operations.",
    rich_help_panel=_TyperRichHelpPanel.ATTRIBUTES,
)
APP.add_typer(
    seed_app.APP,
    name="seed",
    help="Seed sequences operations.",
    rich_help_panel=_TyperRichHelpPanel.ATTRIBUTES,
)

if __name__ == "__main__":
    APP()
