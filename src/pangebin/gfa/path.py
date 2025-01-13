"""Path GFA API wrapper."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pangebin.gfa.line as gfa_line

if TYPE_CHECKING:
    import gfapy  # type: ignore[import-untyped]
    from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
        Path as GfaPath,
    )

_LOGGER = logging.getLogger(__name__)


def get_path_line_by_name(gfa: gfapy.Gfa, name: str) -> GfaPath:
    """Get path by name.

    Parameters
    ----------
    gfa : gfapy.Gfa
        GFA graph
    name : str
        Path name

    Returns
    -------
    GfaPath
        GFA path line

    Raises
    ------
    ValueError
        Invalid path name

    """
    line: gfapy.Line | None = gfa.line(str(name))
    if line is None or line.record_type != gfa_line.Type.PATH:
        _err_msg = f"Invalid path name: {name}, line is {line}"
        _LOGGER.error(_err_msg)
        raise ValueError(_err_msg)
    return line
