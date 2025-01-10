"""GFA header API wrapper."""

from enum import StrEnum

from pangebin.gfa.tag import FieldType


class Tag(StrEnum):
    """GFA header tag types."""

    SKESA_FIX = "FX"
    STANDARDIZED = "Sd"


class TagType(StrEnum):
    """GFA header tag types."""

    SKESA_FIX = FieldType.CHAR
    PREPROCESSED = FieldType.CHAR


class SkesaFixTagValue(StrEnum):
    """Skesa fix header tag values."""

    YES = "Y"
    NO = "N"


class PreprocessedTagValue(StrEnum):
    """Preprocessed header tag values."""

    YES = "Y"
    NO = "N"
