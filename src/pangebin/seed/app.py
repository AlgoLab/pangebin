"""Seed contig applications."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import typer

import pangebin.seed.thresholds.app as seed_thr_app

APP = typer.Typer(rich_markup_mode="rich")

APP.add_typer(seed_thr_app.APP)
