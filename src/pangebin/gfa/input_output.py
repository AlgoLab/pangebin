"""GFA input output module."""

import gzip
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.iter as gfa_iter
import pangebin.gfa.segment as gfa_segment
import pangebin.input_output as io


def from_file(gfa_path: Path) -> gfapy.Gfa:
    """Read a GFA file (compressed or not)."""
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


def gfa_file_to_fasta_file(
    gfa_file: Path,
    fasta_path: Path,
    sep: str = gfa_segment.DEFAULT_ATTRIBUTE_STR_SEP,
) -> None:
    """Write a FASTA file from GFA file.

    Parameters
    ----------
    gfa_file : Path
        GFA graph file
    fasta_path : Path
        Path to FASTA file
    sep : str, optional
        string for separating GFA attributes, default is space

    """
    with fasta_path.open("w") as f_out:
        for seq_record in gfa_iter.sequence_records(gfa_file, sep=sep):
            f_out.write(seq_record.format("fasta"))
