"""GFA items module."""

from enum import StrEnum


class GFALineType(StrEnum):
    """GFA line types."""

    COMMENT = "#"
    HEADER = "H"
    SEGMENT = "S"
    LINK = "L"
    JUMP = "J"
    CONTAINMENT = "C"
    PATH = "P"
    WALK = "W"


class GFAFieldType(StrEnum):
    """GFA field types."""

    CHAR = "A"
    SIGNED_INT = "i"
    FLOAT = "f"
    STRING = "Z"
    JSON = "J"
    BYTE_ARRAY = "H"
    INT_OR_FLOAT_ARRAY = "B"


SKESA_FIX_HEADER_TAG = "FX"
SKESA_FIX_HEADER_TAG_TYPE = GFAFieldType.CHAR


class SkesaFixHeaderTagValue(StrEnum):
    """Skesa fix header tag values."""

    YES = "Y"
    NO = "N"
