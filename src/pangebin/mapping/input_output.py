"""Mapping input-output."""

import logging
from pathlib import Path

import pandas as pd

import pangebin.mapping.items as mapping_items

_LOGGER = logging.getLogger(__name__)


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

    Note
    ----
    The pident value is normalized between 0 and 1.

    """
    try:
        mappings_df = pd.read_csv(
            mapping_file,
            sep="\t",
            names=list(mapping_items.BLAST6_COL_TYPES.keys()),
            dtype=mapping_items.BLAST6_COL_TYPES,
        )
    except Exception:
        _LOGGER.exception("Error reading mapping file %s", mapping_file)
        raise
    #
    # Normalize pident
    #
    for idx, row in mappings_df.iterrows():
        prev_val = row[mapping_items.BlastColumns.PIDENT]
        mappings_df[mapping_items.BlastColumns.PIDENT][idx] = prev_val / 100.0

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
