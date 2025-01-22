"""Pangebin main application module."""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated

import gfapy  # type: ignore[import-untyped]
import pandas as pd
import typer

import pangebin.gfa.app as gfa_app
import pangebin.panassembly.app as panassembly_app
import pangebin.pangenome.app as pangenome_app
import pangebin.plasbin.app as plasbin_app
import pangebin.std_asm_graph.app as std_asm_graph_app


class _TyperRichHelpPanel(StrEnum):
    """Typer rich help panel categories."""

    MAIN = "Main commands"
    UTILS = "Utility commands"


APP = typer.Typer(rich_markup_mode="rich")

APP.add_typer(
    gfa_app.APP,
    name="gfa",
    help="GFA operations.",
    rich_help_panel=_TyperRichHelpPanel.UTILS,
)

APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(std_asm_graph_app.std_asm_graph)
APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(pangenome_app.pangenome)
APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(panassembly_app.panassembly)
APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)(plasbin_app.plasbin)


class ModArgs:
    """ModifyBins arguments."""

    ARG_BIN_FILE = typer.Argument(
        help="Bin file",
    )

    ARG_OUTPUT = typer.Argument(
        help="Output filename",
    )

    ARG_MODTYPE = typer.Argument(
        help="Bin Modification type (NVE= naive) (OVL=overlap)",
    )
    ARG_PANGENOME = typer.Argument(
        help="Pangenome GFA file",
    )


@APP.command(rich_help_panel=_TyperRichHelpPanel.MAIN)
def mod_bins(
    pangenome_graph: Annotated[Path, ModArgs.ARG_PANGENOME],
    bin_file: Annotated[Path, ModArgs.ARG_BIN_FILE],
    output: Annotated[Path, ModArgs.ARG_OUTPUT],
    modtype: Annotated[str, ModArgs.ARG_MODTYPE],
):
    """Modify the bins."""
    bins = pd.read_csv(bin_file, sep="\t")
    gfa = gfapy.Gfa.from_file(pangenome_graph)

    bin_set = {}  # {B1: [contig1, contig2, ...], B2: [contig3, contig4, ...]}
    for i, row in bins.iterrows():
        if row["plasmid"] not in bin_set:
            bin_set[row["plasmid"]] = []
            bin_set[row["plasmid"]].append(row["contig"])
        else:
            bin_set[row["plasmid"]].append(row["contig"])

    # REFACTOR StrEnum for binning_mode
    if modtype == "NVE":
        new_bins = set()
        for bin in bin_set:
            candidates = []
            for bin_2 in bin_set:
                if len(set(bin_set[bin]).intersection(set(bin_set[bin_2]))) > 0:
                    candidates.append(bin_2)
            new_bins.add(tuple(candidates))

        new_bin_set = {}
        counter = 1
        for bin in new_bins:
            contigs = []
            for i in bin:
                contigs += [x for x in bin_set[i]]
            contigs = list(set(contigs))
            new_bin_set[f"B{counter}"] = contigs
            counter += 1

        bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
        for bin in new_bin_set:
            bin_len = 0
            for contig in new_bin_set[bin]:
                gfa_contig = gfa.segment(str(contig))
                contig_len = gfa_contig.LN
                if len(contig_len) > 0:
                    contig_len = contig_len.values[0]
                else:
                    contig_len = -1
                bin_out_csv.loc[len(bin_out_csv)] = {
                    "plasmid": str(bin),
                    "contig": str(contig),
                    "contig_len": contig_len,
                }
        bin_out_csv.to_csv(output, sep="\t", index=False)

    elif modtype == "OVL":
        # augment each bin with the fragments of the pangenome that belongs to the same contig
        new_bins = set()
        for bin in bin_set:
            contigs = set()  ## the fragments of contigs to add to the bin
            fragments = set(bin_set[bin])
            for f in fragments:
                seg = gfa.segment(f)
                if seg != None:
                    contig_list = seg.cl.split(",")
                    for c in contig_list:
                        contigs.add(c)

            candidates = set()
            for c in contigs:
                path = gfa.line(str(c))
                if path is not None:
                    segs = path.segment_names
                    for seg in segs:
                        candidates.add(str(seg.name))
            new_bins.add(tuple(candidates))

        new_bin_set = {}
        counter = 1
        for bin in new_bins:
            contigs = list(set(bin))
            new_bin_set[f"B{counter}"] = contigs
            counter += 1

        bin_out_csv = pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
        for bin in new_bin_set:
            bin_len = 0
            for contig in new_bin_set[bin]:
                gfa_contig = gfa.segment(str(contig))
                contig_len = gfa_contig.LN
                contig_len = contig_len.values[0] if len(contig_len) > 0 else -1
                bin_out_csv.loc[len(bin_out_csv)] = {
                    "plasmid": str(bin),
                    "contig": str(contig),
                    "contig_len": contig_len,
                }

        bin_out_csv.to_csv(output, sep="\t", index=False)

    else:
        pass


if __name__ == "__main__":
    APP()
