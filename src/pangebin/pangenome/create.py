"""Pangenome creation."""

import logging
import subprocess
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
    try:
        subprocess.run(  # noqa: S603
            [str(x) for x in cli_line],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _nfcore_pan_cmd_err = NFCorePangenomeCommandFailedError(exc)
        _LOGGER.critical(str(_nfcore_pan_cmd_err))
        raise _nfcore_pan_cmd_err from exc


class NFCorePangenomeCommandFailedError(Exception):
    """nf-core/pangenome command failed error."""

    def __init__(self, called_proc_exc: subprocess.CalledProcessError) -> None:
        """Initialize."""
        super().__init__()
        self.__called_proc_exc = called_proc_exc

    def called_proc_exc(self) -> subprocess.CalledProcessError:
        """Return the command."""
        return self.__called_proc_exc

    def __str__(self) -> str:
        """Return the error message."""
        return f"nf-core/pangenome command failed: {self.__called_proc_exc.stderr}"


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
