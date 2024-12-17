"""Pangebin preprocess module."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import gfapy
import pandas as pd
import typer

import pangebin.plasbin.app as plasbin_app
from pangebin.graph_utils import (
    add_gfa_to_pangenome,
    clean_pangenome,
    compute_scores,
    convert_kc_to_dp,
    extract_gfagz,
    gfa_to_fasta,
    mix_fasta,
    remove_nodes,
    rename_contigs,
)

APP = typer.Typer(rich_markup_mode="rich")

APP.command(help="Execute PlasBin-Flow.")(plasbin_app.plasbin)


@dataclass
class PreprocessArgs:
    """GFAUtils arguments."""

    ARG_SAMPLE_NAME = typer.Argument(
        help="Sample ID",
    )

    ARG_INPUT_UNI_GFA = typer.Argument(
        help="Unicler GFA assembly graph file",
    )

    ARG_INPUT_SKE_GFA = typer.Argument(
        help="Unicler GFA assembly graph file",
    )

    ARG_OUTPUT_DIR = typer.Option(
        "--outdir",
        # DOCU: here we will place unicycler.gfa skesa.gfa, unicycler.fasta, skesa.fasta, mixed.fasta.gz
        help="Output folder",
    )

    OPT_THRESHOLD = typer.Option(
        "--thr",
        help="Threshold length for removing nodes",
    )


@dataclass
class PanassemblyArgs:
    """Pangenome assembly arguments."""

    ARG_PANGENOME = typer.Argument(
        help="Pangenome GFA file",
    )

    ARG_UNI_PREPROCESSED = typer.Argument(
        help="Unicycler GFA Preprocessed graph",
    )

    ARG_SKE_PREPROCESSED = typer.Argument(
        help="Skesa GFA Preprocessed graph",
    )

    ARG_OUTPUT_DIR = typer.Argument(
        help="Output folder",
    )


@dataclass
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


@APP.command()
def preprocess(
    outdir: Annotated[Path, PreprocessArgs.ARG_OUTPUT_DIR],
    sample: Annotated[str, PreprocessArgs.ARG_SAMPLE_NAME],
    ugfa: Annotated[Path, PreprocessArgs.ARG_INPUT_UNI_GFA],
    sgfa: Annotated[Path, PreprocessArgs.ARG_INPUT_SKE_GFA],
    threshold: Annotated[int, PreprocessArgs.OPT_THRESHOLD],
):
    """Preprocess GFA Assembly files."""
    typer.echo("Preprocessing GFA files")
    outdir.mkdir(parents=True, exist_ok=True)

    ext_u_gfa = extract_gfagz(ugfa)
    ext_s_gfa = extract_gfagz(sgfa)

    r_u_gfa = rename_contigs(ext_u_gfa, "uni")
    r_s_gfa = convert_kc_to_dp(rename_contigs(ext_s_gfa, "ske"))

    rnode_u_gfa = remove_nodes(
        r_u_gfa,
        threshold,
        f"{outdir}/{sample}.{threshold}.u.gfa",
    )
    rnode_s_gfa = remove_nodes(
        r_s_gfa,
        threshold,
        f"{outdir}/{sample}.{threshold}.s.gfa",
    )

    u_fasta = gfa_to_fasta(rnode_u_gfa, f"{outdir}/{sample}.{threshold}.u.fasta")
    s_fasta = gfa_to_fasta(rnode_s_gfa, f"{outdir}/{sample}.{threshold}.s.fasta")
    mix_fasta([u_fasta, s_fasta], f"{outdir}/{sample}.{threshold}.mix.fasta")


@APP.command()
def panassembly(
    pangenome: Annotated[Path, PanassemblyArgs.ARG_PANGENOME],
    skesa_assembly: Annotated[Path, PanassemblyArgs.ARG_SKE_PREPROCESSED],
    unicycler_assembly: Annotated[Path, PanassemblyArgs.ARG_UNI_PREPROCESSED],
    sample: Annotated[str, PreprocessArgs.ARG_SAMPLE_NAME],
    outdir: Annotated[Path, PanassemblyArgs.ARG_OUTPUT_DIR],
):
    """Make a pangenome assembly (panassembly) from the set of original assemblers plus the pangenome."""
    outdir.mkdir(parents=True, exist_ok=True)
    gfa = gfapy.Gfa.from_file(pangenome)
    cl_pangenome = clean_pangenome(gfa)
    assemblers = [skesa_assembly, unicycler_assembly]
    for a in assemblers:
        gfa = gfapy.Gfa.from_file(f"{a}")
        add_gfa_to_pangenome(gfa, cl_pangenome)
    compute_scores(cl_pangenome)
    filename = f"{outdir}/{sample}.panasm.gfa"
    cl_pangenome.to_file(f"{filename}")


@APP.command()
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
