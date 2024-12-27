"""GFA operations."""

from __future__ import annotations

import datetime
import logging
from typing import TYPE_CHECKING

import pangebin.input_output as io
from pangebin.gfa.items import (
    SKESA_FIX_HEADER_TAG,
    SKESA_FIX_HEADER_TAG_TYPE,
    GFAFieldType,
    GFALineType,
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
            f"{datetime.datetime.now(tz=datetime.UTC).isoformat()}"
            f"_{in_gfa_path.name}"
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
        in_gfa_path.unlink()
        out_gfa_path.replace(in_gfa_path)
