"""Mapping iterator module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pangebin.mapping.items as map_items

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


def mapping_from_sam(
    sam_file: Path,
) -> Iterator[map_items.Mapping]:
    """Iterate over a Blast6 mapping file.

    Parameters
    ----------
    sam_file : Path
        Path to the Blast6 SAM mapping file.

    Yields
    ------
    Mapping
        A Mapping object.

    """
    with sam_file.open("r") as f:
        for line in f:
            yield map_items.Mapping.from_string(line)
