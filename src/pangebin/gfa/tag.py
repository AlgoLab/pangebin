"""GFA items module."""

from __future__ import annotations

from enum import StrEnum


class FieldType(StrEnum):
    """GFA field types."""

    CHAR = "A"
    SIGNED_INT = "i"
    FLOAT = "f"
    STRING = "Z"
    JSON = "J"
    BYTE_ARRAY = "H"
    INT_OR_FLOAT_ARRAY = "B"


class IntOrFloatArrayType(StrEnum):
    """Type letter for array of integers or floats."""

    INT8_T = "c"
    UINT8_T = "C"
    INT16_T = "s"
    UINT16_T = "S"
    INT32_T = "i"
    UINT32_T = "I"
    FLOAT = "f"
