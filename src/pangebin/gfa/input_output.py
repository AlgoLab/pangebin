"""GFA input output module."""

import gzip
import shutil
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.iter as gfa_iter
import pangebin.input_output as io


def from_file(gfa_path: Path) -> gfapy.Gfa:
    """Read a GFA file."""
    if io.is_gz_file(gfa_path):
        return from_gfa_gz(gfa_path)
    return gfapy.Gfa.from_file(gfa_path)


def from_gfa_gz(in_gfa_gz_path: Path) -> gfapy.Gfa:
    """Read a GFA gzip compressed file."""
    with (
        gzip.open(in_gfa_gz_path, "rb") as f_in,
        NamedTemporaryFile("wb") as f_out,
    ):
        shutil.copyfileobj(f_in, f_out)  # type: ignore[misc]
        f_out.flush()
        gfa = gfapy.Gfa.from_file(f_out.name, vlevel=0)
    gfa.validate()
    return gfa


def gfa_to_fasta_file(graph: gfapy.Gfa, fasta_path: Path) -> None:
    """Write a FASTA file from GFA graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph
    fasta_path : Path
        Path to FASTA file

    """
    with fasta_path.open("w") as f_out:
        for seq_record in gfa_iter.iter_gfa_to_fasta(graph):
            f_out.write(seq_record.format("fasta"))


def gfa_file_to_fasta_file(gfa_path: Path, fasta_path: Path) -> Path:
    """Convert GFA to FASTA."""
    # REFACTOR may be deprecated
    command = [
        "awk",
        '/^S/{print ">"$2; print $3}',
        str(gfa_path),
        ">",
        str(fasta_path),
    ]
    subprocess.run(command, check=True)  # noqa: S603
    return Path(fasta_path)
