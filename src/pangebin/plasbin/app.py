"""PlasBin application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import typer

import pangebin.plasbin.decomp.app as decomp_app

APP = typer.Typer(rich_markup_mode="rich")
APP.add_typer(decomp_app.APP)
