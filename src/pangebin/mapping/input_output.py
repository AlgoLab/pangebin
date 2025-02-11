"""Mapping input-output."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd
from Bio import SeqIO

import pangebin.mapping.filter as map_filter
import pangebin.mapping.items as map_items

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

_LOGGER = logging.getLogger(__name__)


def iter_mapping_in_file(
    mapping_file: Path,
) -> Iterator[map_items.Mapping]:
    """Iterate over a Blast6 mapping file.

    Parameters
    ----------
    mapping_file : Path
        Path to the Blast6 mapping file.

    Yields
    ------
    Mapping
        A Mapping object.

    """
    with mapping_file.open("r") as f:
        for line in f:
            yield map_items.Mapping.from_string(line)


def to_filtered_sam_file(
    input_sam: Path,
    filtered_sam: Path,
    config: map_filter.Config,
    query_fasta: Path | None = None,
    subject_fasta: Path | None = None,
) -> None:
    """Convert a SAM file to a filtered SAM file.

    Parameters
    ----------
    input_sam : Path
        Path to the input SAM file.
    filtered_sam : Path
        Path to the filtered SAM file.
    config : map_filter.Config
        Config object.
    query_fasta : Path, optional
        Path to the query FASTA file.
    subject_fasta : Path, optional
        Path to the subject FASTA file.

    """
    q_len_dict: dict[str, int] | None = (
        None
        if query_fasta is None
        else {record.id: len(record) for record in SeqIO.parse(query_fasta, "fasta")}
    )
    s_len_dict: dict[str, int] | None = (
        None
        if subject_fasta is None
        else {record.id: len(record) for record in SeqIO.parse(subject_fasta, "fasta")}
    )
    _number_of_original_mapping = 0
    _number_of_kept_mapping = 0
    with input_sam.open("r") as f_in, filtered_sam.open("w") as f_out:
        for line in f_in:
            _number_of_original_mapping += 1
            mapping = map_items.Mapping.from_string(line)
            if map_filter.accept_mapping(mapping, config, q_len_dict, s_len_dict):
                _number_of_kept_mapping += 1
                f_out.write(line)

    _LOGGER.info(
        "Number of kept mappings: %i/%i (%.2f%%)",
        _number_of_kept_mapping,
        _number_of_original_mapping,
        _number_of_kept_mapping / _number_of_original_mapping * 100,
    )
    _LOGGER.info(
        "Filtered SAM file: %s",
        filtered_sam,
    )


# REFACTOR do I need dataframe?
# FIXME no reorder coordinates
def mapping_file_to_dataframe(
    mapping_file: Path,
) -> pd.DataFrame:
    """Convert a Blast6 mapping file to a DataFrame.

    Parameters
    ----------
    mapping_file : Path
        Path to the Blast6 mapping file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the Blast6 mapping data.

    """
    mappings_df = pd.DataFrame(
        columns=list(map_items.BLAST6_COL_TYPES.keys()),
        dtype=map_items.BLAST6_COL_TYPES,
    )
    try:
        for mapping in iter_mapping_in_file(mapping_file):
            mappings_df.loc[len(mappings_df)] = mapping.to_list()
    except Exception:
        _LOGGER.exception("Error reading mapping file %s", mapping_file)
        raise

    return mappings_df


class ReadMappingFileError(Exception):
    """Read mapping file error."""

    def __init__(self, mapping_file: Path) -> None:
        """Initialize object."""
        super().__init__(f"Error reading mapping file {mapping_file}")
        self.__mapping_file = mapping_file

    def mapping_file(self) -> Path:
        """Get mapping file."""
        return self.__mapping_file
