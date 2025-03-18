"""Seed sequence input/output module."""

from collections.abc import Iterable
from pathlib import Path


def to_tsv(sequence_ids: Iterable[str], out_file: Path) -> None:
    """Write seed sequences to file.

    Parameters
    ----------
    sequence_ids : iterable of str
        Iterable of sequence ids
    out_file : Path
        Path to output file
    """
    with out_file.open("w") as f_out:
        for sequence_id in sequence_ids:
            f_out.write(sequence_id)
            f_out.write("\n")
