"""GFA operations."""

from __future__ import annotations

import datetime
import logging
import shutil
from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]

import pangebin.gfa.header as gfa_header
import pangebin.gfa.line as gfa_line
import pangebin.gfa.tag as gfa_tag
import pangebin.input_output as io
from pangebin.gfa import segment as gfa_segment

if TYPE_CHECKING:
    from pathlib import Path

    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


_LOGGER = logging.getLogger(__name__)


def is_skesa_gfa_fixed(gfa_path: Path) -> bool:
    """Check ig the Skeza GFA file is fixed."""
    yes_fix_tag = (
        f"{gfa_header.Tag.SKESA_FIX}"
        f":{gfa_header.TagType.SKESA_FIX}"
        f":{gfa_header.SkesaFixTagValue.YES}"
    )
    with io.open_file_read(gfa_path) as f_in:
        for line in f_in:
            if line.startswith(
                str(gfa_line.Type.HEADER),
            ) and yes_fix_tag in line.rstrip(
                "\n",
            ).split("\t"):
                _LOGGER.debug("Skesa GFA file is fixed.")
                return True
    _LOGGER.debug("Skesa GFA file is not fixed.")
    return False


def fix_skesa_gfa(
    in_gfa_path: Path,
    out_gfa_path: Path | None = None,
) -> None:
    """Fix a Skeza GFA file.

    Parameters
    ----------
    in_gfa_path : Path
        Path to input GFA file
    out_gfa_path : Path | None, optional
        Path to output GFA file, must be different from input if used, by default None


    Warnings
    --------
    This function modifies the input GFA file if out_gfa_path is None.

    Raises
    ------
    ValueError
        If input and output files are the same

    """
    replace_file = out_gfa_path is None
    if replace_file:
        _LOGGER.debug("Fixing Skeza GFA file %s", in_gfa_path)
    else:
        _LOGGER.debug("Fixing Skeza GFA file %s to %s.", in_gfa_path, out_gfa_path)

    if out_gfa_path is None:
        out_gfa_path = in_gfa_path.parent / (
            f"{datetime.datetime.now(tz=datetime.UTC).isoformat()}.gfa"
        )
        _LOGGER.debug("Temporary output file: %s", out_gfa_path)
    elif in_gfa_path == out_gfa_path:
        _err_msg = f"Input and output files are the same: {in_gfa_path}"
        _LOGGER.error(_err_msg)
        raise ValueError(_err_msg)

    yes_fix_tag = (
        f"{gfa_header.Tag.SKESA_FIX}"
        f":{gfa_header.TagType.SKESA_FIX}"
        f":{gfa_header.SkesaFixTagValue.YES}"
    )
    with (
        io.open_file_read(in_gfa_path) as f_in,
        io.open_file_write(out_gfa_path) as f_out,
    ):
        f_out.write(f"{gfa_line.Type.HEADER}\t{yes_fix_tag}\n")
        for line in f_in:
            if line.startswith(gfa_line.Type.SEGMENT):
                split_line = line.split("\t")
                f_out.write(split_line[0])  # S
                f_out.write("\t")
                f_out.write(split_line[1])  # segment name
                f_out.write("\t")
                f_out.write(split_line[2])  # sequence
                for optional_field in split_line[3:]:
                    tag_name, tag_type, value = optional_field.split(":")
                    if tag_type == gfa_tag.FieldType.SIGNED_INT:
                        f_out.write(f"\t{tag_name}:{tag_type}:{int(float(value))}")
                    else:
                        f_out.write(f"\t{tag_name}:{tag_type}:{value}")
                f_out.write("\n")
            else:
                f_out.write(line)

    if replace_file:
        _LOGGER.debug(
            "Replacing temporary output file %s with input file %s.",
            out_gfa_path,
            in_gfa_path,
        )
        in_gfa_path.unlink()
        with (
            io.open_file_read(out_gfa_path) as f_in,
            io.open_file_write(in_gfa_path) as f_out,
        ):
            shutil.copyfileobj(f_in, f_out)
        out_gfa_path.unlink()


def set_segment_length_tags(graph: gfapy.Gfa) -> None:
    """Set the segment length attribute in a GFA graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting segment length tags.")
    segment: GfaSegment
    for segment in graph.segments:
        segment.set_datatype(
            gfa_segment.Tag.LENGTH,
            gfa_segment.TagType.LENGTH,
        )
        if len(segment.sequence) == 0:
            segment.sequence = gfapy.Placeholder()
            gfa_segment.set_length(segment, 0)
        else:
            gfa_segment.set_length(segment)


def rename_segments(
    graph: gfapy.Gfa,
    segment_prefix_name: str,
) -> None:
    """Rename segments in a GFA graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph
    segment_prefix_name : str
        Prefix for segment names

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Renaming segments.")
    segment: GfaSegment
    for counter, segment in enumerate(graph.segments):
        segment.name = gfa_segment.format_name(segment_prefix_name, counter + 1)


def convert_kmer_coverage_to_normalized_coverage(
    graph: gfapy.Gfa,
) -> None:
    """Convert k-mer coverage to normalized coverage.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Converting k-mer coverage to normalized coverage.")
    total_coverage = sum(gfa_segment.kmer_coverage(seg) for seg in graph.segments)
    total_length = sum(gfa_segment.length(seg) for seg in graph.segments)

    seg: GfaSegment
    for seg in graph.segments:
        seg.set_datatype(
            gfa_segment.Tag.NORMALIZED_COVERAGE,
            gfa_segment.TagType.NORMALIZED_COVERAGE,
        )
        gfa_segment.set_normalized_coverage(
            seg,
            float(
                (gfa_segment.kmer_coverage(seg) * total_length)
                / (gfa_segment.length(seg) * total_coverage),
            ),
        )
        seg.delete(gfa_segment.Tag.KMER_COVERAGE)


def set_standardized_header_tag(gfa: gfapy.Gfa) -> None:
    """Set the standardized header tag in the GFA graph."""
    _LOGGER.info("Setting standardized header tag.")
    gfa.header.set_datatype(
        gfa_header.Tag.STANDARDIZED,
        gfa_header.TagType.STANDARDIZED,
    )
    gfa.header.add(gfa_header.Tag.STANDARDIZED, gfa_header.StandardizedTagValue.YES)


def is_standardized(gfa: gfapy.Gfa) -> bool:
    """Check if a GFA graph is standardized."""
    return (
        gfa.header.get(gfa_header.Tag.STANDARDIZED)
        == gfa_header.StandardizedTagValue.YES
    )
