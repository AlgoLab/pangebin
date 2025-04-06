"""PlasBin application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import typer

import pangebin.plasbin.binlab.app as binlab_app
import pangebin.plasbin.decomp.app as decomp_app
import pangebin.plasbin.once.app as once_app

APP = typer.Typer(rich_markup_mode="rich")
APP.add_typer(decomp_app.APP)
APP.add_typer(binlab_app.APP)
APP.add_typer(once_app.APP)
