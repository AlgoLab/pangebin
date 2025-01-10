"""Path GFA API wrapper."""

from __future__ import annotations

import logging
from enum import StrEnum
from typing import TYPE_CHECKING

import pangebin.gfa.line as gfa_line
import pangebin.gfa.tag as gfa_tag

if TYPE_CHECKING:
    import gfapy  # type: ignore[import-untyped]
    from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
        Path as GfaPath,
    )

_LOGGER = logging.getLogger(__name__)


class Tag(StrEnum):
    """Path tags."""

    LENGTH = "LN"
    FROM_ASSEMBLER = "aa"
    NORMALIZED_COVERAGE = "dp"


class TagType(StrEnum):
    """Path tag types."""

    LENGTH = gfa_tag.FieldType.SIGNED_INT
    FROM_ASSEMBLER = gfa_tag.FieldType.CHAR
    NORMALIZED_COVERAGE = gfa_tag.FieldType.FLOAT


class FromAssemblerTagValue(StrEnum):
    """From assembler tag values."""

    UNICYCLER_CONTIG = "u"
    SKESA_CONTIG = "s"


def length(path: GfaPath) -> int:
    """Get path length."""
    return path.get(Tag.LENGTH)


def set_length(path: GfaPath, length: int) -> None:
    """Set path length."""
    path.set(Tag.LENGTH, length)


def from_assembler(path: GfaPath) -> FromAssemblerTagValue:
    """Get from assembler tag value."""
    return FromAssemblerTagValue(path.get(Tag.FROM_ASSEMBLER))


def set_from_assembler(path: GfaPath, value: FromAssemblerTagValue) -> None:
    """Set from assembler tag value."""
    path.set(Tag.FROM_ASSEMBLER, value)


def normalized_coverage(path: GfaPath) -> float:
    """Get normalized coverage."""
    return path.get(Tag.NORMALIZED_COVERAGE)


def set_normalized_coverage(path: GfaPath, coverage: float) -> None:
    """Set normalized coverage."""
    path.set(Tag.NORMALIZED_COVERAGE, coverage)


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
