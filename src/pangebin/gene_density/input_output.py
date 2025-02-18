"""Gene density input-output module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pangebin.mapping.items as map_items

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


# REFACTOR ouput only non-0 gene density (i.e. with mapping intervals)
# It will enable to reduce the size of the dictionaries(?)
def to_file_with_intervals(
    file: Path,
    gene_densities: dict[str, float],
    gene_mapping_intervals: dict[str, list[map_items.MappingInterval]] | None = None,
) -> None:
    """Write gene density (and intervals if provided) to file.

    Parameters
    ----------
    file : Path
        Path to the output file.
    gene_densities : dict[str, float]
        Gene density for each sequence.
    gene_mapping_intervals : dict[str, list[MappingInterval]], optional
        Gene intervals for each sequence, by default None
    """
    with file.open("w") as f_out:
        if gene_mapping_intervals is not None:
            for sequence_id, gene_density in gene_densities.items():
                interval_str = (
                    ""
                    if sequence_id not in gene_mapping_intervals
                    else "\t"
                    + ",".join(
                        str(interval)
                        for interval in gene_mapping_intervals[sequence_id]
                    )
                )
                f_out.write(
                    f"{sequence_id}\t{gene_density}" + interval_str + "\n",
                )
        else:
            for sequence_id, gene_density in gene_densities.items():
                f_out.write(f"{sequence_id}\t{gene_density}\n")


def from_file(file: Path) -> Iterator[tuple[str, float]]:
    """Read gene densities from file.

    Yields
    ------
    str
        Sequence name
    float
        Its gene density
    """
    with file.open() as f_in:
        for line in f_in:
            str_items: list[str] = line.split()
            yield str_items[0], float(str_items[1])


def from_file_with_intervals(
    file: Path,
) -> Iterator[
    tuple[
        str,
        float,
        list[map_items.MappingInterval],
    ]
]:
    """Read gene density and intervals from file.

    Yields
    ------
    str
        Sequence name
    float
        Its gene density
    list[MappingInterval]
        Its intervals
    """
    with file.open() as f_in:
        for line in f_in:
            str_items: list[str] = line.split()
            yield (
                str_items[0],
                float(str_items[1]),
                (
                    [
                        map_items.MappingInterval.from_string(interval_str)
                        for interval_str in str_items[2].split(",")
                    ]
                    if len(str_items) > 2  # noqa: PLR2004
                    else []
                ),
            )
