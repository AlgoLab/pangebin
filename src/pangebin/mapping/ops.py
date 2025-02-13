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
        interval_start, interval_end = (
            (mapping.sstart(), mapping.send())
            if mapping.sstart() < mapping.send()
            else (mapping.send(), mapping.sstart())
        )
        interval = map_items.MappingInterval(
            interval_start,
            interval_end,
        )
        if mapping.sseqid() not in subject_intervals:
            subject_intervals[mapping.sseqid()] = [interval]
        else:
            subject_intervals[mapping.sseqid()].append(interval)
    return subject_intervals


def interval_union(
    intervals: list[map_items.MappingInterval],
) -> list[map_items.MappingInterval]:
    """Return the union of a list of intervals.

    Parameters
    ----------
    intervals : list[MappingInterval]
        List of MappingIntervals

    Returns
    -------
    list[MappingInterval]
        List of MappingIntervals

    """
    union_intervals: list[map_items.MappingInterval] = []
    start_sorted_intervals = sorted(intervals, key=lambda x: x.start())
    if start_sorted_intervals:
        union_intervals.append(start_sorted_intervals[0])
    for interval in start_sorted_intervals[1:]:
        if union_intervals[-1].end() >= interval.start() - 1:
            union_intervals[-1].set_end(max(union_intervals[-1].end(), interval.end()))
        else:
            union_intervals.append(interval)
    return union_intervals
