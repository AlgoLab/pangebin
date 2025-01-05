"""GFA operations."""

from __future__ import annotations

import datetime
import logging
import shutil
from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]
from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]

import pangebin.input_output as io
from pangebin.assembler import ContigPrefix
from pangebin.gfa.items import (
    SKESA_FIX_HEADER_TAG,
    SKESA_FIX_HEADER_TAG_TYPE,
    GFAFieldType,
    GFALineType,
    SequenceTag,
    SequenceTagType,
    SkesaFixHeaderTagValue,
)

if TYPE_CHECKING:
    from pathlib import Path

_LOGGER = logging.getLogger(__name__)


def is_skesa_gfa_fixed(gfa_path: Path) -> bool:
    """Check ig the Skeza GFA file is fixed."""
    yes_fix_tag = (
        f"{SKESA_FIX_HEADER_TAG}"
        f":{SKESA_FIX_HEADER_TAG_TYPE}"
        f":{SkesaFixHeaderTagValue.YES}"
    )
    with io.open_file_read(gfa_path) as f_in:
        for line in f_in:
            if line.startswith(str(GFALineType.HEADER)) and yes_fix_tag in line.rstrip(
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
        f"{SKESA_FIX_HEADER_TAG}"
        f":{SKESA_FIX_HEADER_TAG_TYPE}"
        f":{SkesaFixHeaderTagValue.YES}"
    )
    with (
        io.open_file_read(in_gfa_path) as f_in,
        io.open_file_write(out_gfa_path) as f_out,
    ):
        f_out.write(f"{GFALineType.HEADER}\t{yes_fix_tag}\n")
        for line in f_in:
            if line.startswith(GFALineType.SEGMENT):
                split_line = line.split("\t")
                f_out.write(split_line[0])  # S
                f_out.write("\t")
                f_out.write(split_line[1])  # segment name
                f_out.write("\t")
                f_out.write(split_line[2])  # sequence
                for optional_field in split_line[3:]:
                    tag_name, tag_type, value = optional_field.split(":")
                    if tag_type == GFAFieldType.SIGNED_INT:
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


def rename_contigs(
    graph: gfapy.Gfa,
    segment_name_prefix: ContigPrefix,
) -> None:
    """Rename contigs in a GFA graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph
    segment_name_prefix : ContigPrefix
        Prefix for contig names

    Warnings
    --------
    This function mutates the GFA graph

    """
    _old_level = graph.vlevel
    graph.vlevel = 3

    seg: GfaSegment
    for counter, seg in enumerate(graph.segments):
        if seg.LN is not None:
            if seg.LN == 0:
                seg.sequence = gfapy.Placeholder()
                seg.LN = 0
        elif len(seg.sequence) == 0:
            seg.sequence = gfapy.Placeholder()
            seg.LN = 0
        seg.name = f"{segment_name_prefix}{counter + 1}"

    graph.vlevel = _old_level


def convert_kc_to_dp(
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
    total_coverage = 0
    total_length = 0

    _old_level = graph.vlevel
    graph.vlevel = 3

    seg: GfaSegment
    for seg in graph.segments:
        total_coverage += seg.KC
        if seg.LN is None:
            seg.set_datatype(SequenceTag.LENGTH, SequenceTagType.LENGTH)
            seg.LN = len(seg.sequence)
        total_length += seg.LN

    for seg in graph.segments:
        seg.set_datatype(
            SequenceTag.NORMALIZED_COVERAGE,
            SequenceTagType.NORMALIZED_COVERAGE,
        )
        seg.dp = float((seg.KC * total_length) / (seg.LN * total_coverage))
        seg.KC = None

    graph.vlevel = _old_level
