"""FASTA module."""

import shutil
from pathlib import Path

from Bio import SeqIO


def remove_small_sequences(
    in_fasta_path: Path,
    min_contig_length: int,
    out_fasta_path: Path,
) -> None:
    """Remove small sequences from a FASTA file."""
    with Path(in_fasta_path).open() as spe, Path(out_fasta_path).open("w") as out:
        for record in SeqIO.parse(spe, "fasta"):
            if len(record) >= min_contig_length:
                out.write(f">{record.id}\n")
                out.write(f"{record.seq!s}\n")


def merge_two_fastas(
    first_fasta_path: Path,
    second_fasta_path: Path,
    merged_fasta_path: Path,
) -> None:
    """Merge two FASTA files into one.

    Parameters
    ----------
    first_fasta_path : Path
        First FASTA file.
    second_fasta_path : Path
        Second FASTA file.
    merged_fasta_path : Path
        Merged FASTA file.

    """
    with merged_fasta_path.open("w") as f_out:
        for fasta in [first_fasta_path, second_fasta_path]:
            with fasta.open() as f_in:
                shutil.copyfileobj(f_in, f_out)
