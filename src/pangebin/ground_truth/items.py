"""Ground truth items."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


class PlasmidBioSampleNCBI:
    """Plasmid bio sample NCBI."""

    KEY_BIOSAMPLE_ID = "BioSample_ID"
    KEY_SRA_ID = "SRA_ID"
    KEY_GENBANK_IDS = "Plasmid_GenBank_IDs"

    @classmethod
    def from_dict(cls, dict_from_yaml: dict) -> PlasmidBioSampleNCBI:
        """Create object from dict."""
        return cls(
            dict_from_yaml[cls.KEY_BIOSAMPLE_ID],
            dict_from_yaml[cls.KEY_SRA_ID],
            dict_from_yaml[cls.KEY_GENBANK_IDS],
        )

    def __init__(
        self,
        biosample_id: str,
        sra_id: str,
        genbank_ids: Iterable[str],
    ) -> None:
        self.__biosample_id = biosample_id
        self.__sra_id = sra_id
        self.__genbank_ids = list(genbank_ids)

    def biosample_id(self) -> str:
        """Biosample ID."""
        return self.__biosample_id

    def sra_id(self) -> str:
        """SRA ID."""
        return self.__sra_id

    def genbank_ids(self) -> list[str]:
        """GenBank IDs."""
        return self.__genbank_ids

    def to_dict(self) -> dict[str, str | list[str]]:
        """Convert to dict."""
        return {
            self.KEY_BIOSAMPLE_ID: self.__biosample_id,
            self.KEY_SRA_ID: self.__sra_id,
            self.KEY_GENBANK_IDS: self.__genbank_ids,
        }


class PlasmidContig:
    """Plasmid contig."""

    @classmethod
    def from_tsv_row(cls, tsv_row: str) -> PlasmidContig:
        """Create object from TSV row."""
        items_str: list[str] = tsv_row.split()
        return cls(
            items_str[0],
            items_str[1],
            int(items_str[2]),
            int(items_str[3]),
            float(items_str[4]),
        )

    def __init__(
        self,
        plasmid_id: str,
        contig_id: str,
        plasmid_length: int,
        contig_length: int,
        contig_coverage: float,
    ) -> None:
        self.__data = (
            plasmid_id,
            contig_id,
            plasmid_length,
            contig_length,
            contig_coverage,
        )

    def plasmid_id(self) -> str:
        """Plasmid ID."""
        return self.__data[0]

    def contig_id(self) -> str:
        """Contig ID."""
        return self.__data[1]

    def plasmid_length(self) -> int:
        """Plasmid length."""
        return self.__data[2]

    def contig_length(self) -> int:
        """Contig length."""
        return self.__data[3]

    def contig_coverage(self) -> float:
        """Contig coverage."""
        return self.__data[4]

    def to_tsv_row(self) -> str:
        """Convert to TSV row."""
        return "\t".join(
            (str(item) for item in self.__data),
        )


class NonPlasmidContig:
    """Non plasmid contig."""

    @classmethod
    def from_tsv_row(cls, tsv_row: str) -> NonPlasmidContig:
        """Create object from TSV row."""
        items_str: list[str] = tsv_row.split()
        return cls(items_str[0], int(items_str[1]))

    def __init__(self, contig_id: str, contig_length: int) -> None:
        """Initialize object."""
        self.__data = (contig_id, contig_length)

    def contig_id(self) -> str:
        """Contig ID."""
        return self.__data[0]

    def contig_length(self) -> int:
        """Contig length."""
        return self.__data[1]

    def to_tsv_row(self) -> str:
        """Convert to TSV row."""
        return "\t".join(
            (str(item) for item in self.__data),
        )
