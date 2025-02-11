"""GFA input output module."""

import gzip
import shutil
import subprocess
from collections.abc import Iterator
from pathlib import Path
from tempfile import NamedTemporaryFile

import gfapy  # type: ignore[import-untyped]
from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]

import pangebin.gfa.iter as gfa_iter
import pangebin.gfa.line as gfa_line
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


def iter_segment(gfa_file: Path) -> Iterator[GfaSegment]:
    """Get a segment line iterator."""
    with io.open_file_read(gfa_file) as f_in:
        yield from (line for line in f_in if line[0] == gfa_line.Type.SEGMENT)


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
        for seq_record in gfa_iter.sequence_records(graph):
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
