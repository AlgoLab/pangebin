"""Plasmid database items."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


class IlluminaBioSamples:
    """BioSample, their plasmids and reads."""

    KEY_BIOSAMPLE_ID = "BioSampleID"
    KEY_SAMPLE_SRA_ID = "SampleSRAID"
    KEY_PLASMID_GENBANK_IDS = "PlasmidGenBankIDs"
    KEY_RUN_SRA_IDS = "RunSRAIDs"

    @classmethod
    def from_dict(cls, dict_from_yaml: dict) -> IlluminaBioSamples:
        """Create object from dict."""
        return cls(
            dict_from_yaml[cls.KEY_BIOSAMPLE_ID],
            dict_from_yaml[cls.KEY_SAMPLE_SRA_ID],
            dict_from_yaml[cls.KEY_PLASMID_GENBANK_IDS],
            dict_from_yaml[cls.KEY_RUN_SRA_IDS],
        )

    def __init__(
        self,
        biosample_id: str,
        sample_sra_id: str,
        plasmids: Iterable[str],
        run_sra_ids: Iterable[str],
    ) -> None:
        """Initialize object.

        Parameters
        ----------
        biosample_id : str
            BioSample ID
        sample_sra_id : str
            Sample SRA ID
        plasmids : iterable of str
            Set of plasmid Genbank IDs
        run_sra_ids : iterable of str
            Set of reads run SRA IDs
        """
        self.__biosample_id = biosample_id
        self.__sample_sra_id = sample_sra_id
        self.__plasmids = list(plasmids)
        self.__run_sra_ids = list(run_sra_ids)

    def biosample_id(self) -> str:
        """Get BioSample ID."""
        return self.__biosample_id

    def sample_sra_id(self) -> str:
        """Get Sample SRA ID."""
        return self.__sample_sra_id

    def plasmids(self) -> list[str]:
        """Get plasmid Genbank IDs."""
        return self.__plasmids

    def run_sra_ids(self) -> list[str]:
        """Get run SRA reads IDs."""
        return self.__run_sra_ids

    def to_dict(self) -> dict[str, str | list[str]]:
        """Convert to dict."""
        return {
            self.KEY_BIOSAMPLE_ID: self.__biosample_id,
            self.KEY_SAMPLE_SRA_ID: self.__sample_sra_id,
            self.KEY_PLASMID_GENBANK_IDS: self.__plasmids,
            self.KEY_RUN_SRA_IDS: self.__run_sra_ids,
        }


class NonIlluminaBioSamples:
    """BioSample and their plasmids without Illumina reads."""

    KEY_BIOSAMPLE_ID = "BioSampleID"
    KEY_SAMPLE_SRA_ID = "SampleSRAID"
    KEY_PLASMID_GENBANK_IDS = "PlasmidGenBankIDs"

    @classmethod
    def from_dict(cls, dict_from_yaml: dict) -> NonIlluminaBioSamples:
        """Create object from dict."""
        return cls(
            dict_from_yaml[cls.KEY_BIOSAMPLE_ID],
            dict_from_yaml[cls.KEY_SAMPLE_SRA_ID],
            dict_from_yaml[cls.KEY_PLASMID_GENBANK_IDS],
        )

    def __init__(
        self,
        biosample_id: str,
        sample_sra_id: str,
        plasmids: Iterable[str],
    ) -> None:
        """Initialize object.

        Parameters
        ----------
        biosample_id : str
            BioSample ID
        sample_sra_id : str
            Sample SRA ID
        plasmids : iterable of str
            Set of plasmid Genbank IDs
        """
        self.__biosample_id = biosample_id
        self.__sample_sra_id = sample_sra_id
        self.__plasmids = list(plasmids)

    def biosample_id(self) -> str:
        """Get BioSample ID."""
        return self.__biosample_id

    def sample_sra_id(self) -> str:
        """Get Sample SRA ID."""
        return self.__sample_sra_id

    def plasmids(self) -> list[str]:
        """Get plasmid Genbank IDs."""
        return self.__plasmids

    def to_dict(self) -> dict[str, str | list[str]]:
        """Convert to dict."""
        return {
            self.KEY_BIOSAMPLE_ID: self.__biosample_id,
            self.KEY_SAMPLE_SRA_ID: self.__sample_sra_id,
            self.KEY_PLASMID_GENBANK_IDS: self.__plasmids,
        }
