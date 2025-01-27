"""Pangenome creation."""

import logging
from pathlib import Path

import gfapy  # type: ignore[import-untyped]

import pangebin.input_output as common_io
import pangebin.pangenome.config as pangenome_config
import pangebin.std_asm_graph.fasta as std_asm_graph_fasta
from pangebin import subprocess_lib

_LOGGER = logging.getLogger(__name__)

_NUMBER_OF_HAPLOTYPES = 2

NEXTFLOW_CMD_STR = "nextflow"


def nfcore_pangenome(
    mixed_fasta_gz_path: Path,
    config: pangenome_config.Pangenome,
    out_nfcore_dir_path: Path,
) -> None:
    """Create nf-core/pangenome assembly."""
    cmd_path = subprocess_lib.command_path(NEXTFLOW_CMD_STR)
    cli_line = [
        cmd_path,
        "run",
        "nf-core/pangenome",
        "-r",
        config.release(),
        #
        # nf-core options
        #
        "-profile",
        config.profile(),
        ("-resume" if config.resume() else ""),
        #
        # nf-core/pangenome options
        #
        "--input",
        mixed_fasta_gz_path,
        "--n_haplotypes",
        _NUMBER_OF_HAPLOTYPES,
        "--seqwish_min_match_length",
        pangenome_config.SEQWISH_MIN_MATCH_LENGTH,
        "--smoothxg_poa_params",
        pangenome_config.SMOOTHXG_POA_PARAMS,
        "--wfmash_segment_length",
        pangenome_config.WFMASH_SEGMENT_LENGTH,
        "--wfmash_chunks",
        pangenome_config.WFMASH_CHUNKS,
        "--wfmash_map_pct_id",
        pangenome_config.WFMASH_MAP_PCT_ID,
        "--outdir",
        out_nfcore_dir_path,
    ]
    if config.supp_nfcore_pangenome_config_path() is not None:
        cli_line.extend(("-params-file", config.supp_nfcore_pangenome_config_path()))

    subprocess_lib.run_cmd(cli_line, NEXTFLOW_CMD_STR)


def add_false_sequence(mixed_fasta_gz_path: Path) -> None:
    """Add false sequence."""
    # XXX tmp fix of https://github.com/nf-core/pangenome/issues/215#issuecomment-2593269805
    with common_io.open_file_append(mixed_fasta_gz_path) as f:
        f.write(">none\n")
        f.write("*\n")


def rename_gfa_paths(pangenome_gfa_path: Path) -> None:
    """Rename pangenome GFA paths."""
    pangenome_gfa = gfapy.Gfa.from_file(pangenome_gfa_path)
    for path in pangenome_gfa.paths:
        path.name = std_asm_graph_fasta.pansn_name_to_contig_name(path.name)
    pangenome_gfa.to_file(pangenome_gfa_path)
