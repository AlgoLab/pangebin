"""GFA iterator module."""

from __future__ import annotations

from typing import TYPE_CHECKING

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gfapy  # type: ignore[import-untyped]
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


def __gfa_segment_attribute_to_string(segment: GfaSegment, sep: str = " ") -> str:
    """Format GFA segment attributes as a string.

    Parameters
    ----------
    segment: GfaSegment
        GFA segment
    sep: str, optional
        string for separating GFA attributes, default is space

    Returns
    -------
    str
        GFA segment attributes in format sep.join(key:value)

    """
    return sep.join(
        [f"{tag_name}:{segment.get(tag_name)}" for tag_name in segment.tagnames],
    )


def iter_gfa_to_fasta(
    graph: gfapy.Gfa,
    sep: str = " ",
) -> Iterator[SeqRecord]:
    """Iterate over contigs and their attributes from a GFA files.

    Parameters
    ----------
    graph: gfapy.Gfa
        GFA graph
    sep: str, optional
        string for separating GFA attributes, default is space

    Yields
    ------
    SeqRecord
        FASTA record

    Notes
    -----
    The FASTA header format is:
    <contig name> <contig name>.GFA <attributes string=sep.join(key:value)>

    """
    seg: GfaSegment
    for seg in graph.segments:
        yield SeqRecord(
            Seq(seg.sequence),
            id=seg.name,
            name=seg.name,
            description=f"{seg.name}.GFA {__gfa_segment_attribute_to_string(seg, sep=sep)}",  # noqa: E501
        )
