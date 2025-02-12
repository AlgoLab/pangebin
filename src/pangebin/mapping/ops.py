"""Mapping operations."""

from pathlib import Path

import pangebin.mapping.items as map_items
import pangebin.mapping.iter as map_iter


def subject_intervals_from_sam(
    sam_file: Path,
) -> dict[str, list[map_items.MappingInterval]]:
    """Return the ammping intervals for each subject in a SAM file.

    Parameters
    ----------
    sam_file : Path
        Path to the Blast6 SAM mapping file.

    Returns
    -------
    dict[str, list[MappingInterval]]
        Subject ID (str): List of MappingIntervals

    """
    subject_intervals: dict[str, list[map_items.MappingInterval]] = {}
    for mapping in map_iter.mapping_from_sam(sam_file):
        interval = map_items.MappingInterval(
            mapping.qseqid(),
            mapping.sstart(),
            mapping.send(),
        )
        if mapping.sseqid() not in subject_intervals:
            subject_intervals[mapping.sseqid()] = [interval]
        else:
            subject_intervals[mapping.sseqid()].append(interval)
    return subject_intervals
