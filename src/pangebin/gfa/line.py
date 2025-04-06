"""GFA line API wrapper."""

from __future__ import annotations

from enum import StrEnum


class Type(StrEnum):
    """GFA line types."""

    COMMENT = "#"
    HEADER = "H"
    SEGMENT = "S"
    LINK = "L"
    JUMP = "J"
    CONTAINMENT = "C"
    PATH = "P"
    WALK = "W"
