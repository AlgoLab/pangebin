"""Ground truth input-output module."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper
import logging

import pangebin.ground_truth.items as gt_items

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

_LOGGER = logging.getLogger(__name__)


# REFACTOR unused function
def plasmid_biosamples_from_file(
    yaml_file: Path,
) -> Iterator[gt_items.PlasmidBioSampleNCBI]:
    """Get plasmid biosamples from file.

    Parameters
    ----------
    yaml_file : Path
        Path to the YAML file

    Yield
    -----
    PlasmidBioSampleNCBI
        Plasmid BioSample from NCBI
    """
    with yaml_file.open() as f_in:
        for dict_from_yaml in yaml.safe_load_all(f_in):
            yield gt_items.PlasmidBioSampleNCBI.from_dict(dict_from_yaml)


def plasmid_contigs_to_file(
    plasmid_contigs: Iterable[gt_items.PlasmidContig],
    out_file: Path,
) -> None:
    """Write plasmid contigs to file.

    Parameters
    ----------
    plasmid_contigs : iterable of PlasmidContig
        Iterable of plasmid contigs
    out_file : Path
        Path to output file
    """
    _LOGGER.info("Write plasmid contigs file: %s", out_file)
    with out_file.open("w") as f_out:
        for plasmid_contig in plasmid_contigs:
            f_out.write(plasmid_contig.to_tsv_row())
            f_out.write("\n")


def plasmid_contigs_from_file(file: Path) -> Iterator[gt_items.PlasmidContig]:
    """Get plasmid contigs from file.

    Parameters
    ----------
    file : Path
        Path to the file

    Yield
    -----
    PlasmidContig
        Plasmid contig
    """
    with file.open() as f_in:
        for line in f_in:
            yield gt_items.PlasmidContig.from_tsv_row(line)


def non_plasmid_contigs_to_file(
    non_plasmid_contigs: Iterable[gt_items.NonPlasmidContig],
    out_file: Path,
) -> None:
    """Write non plasmid contigs to file.

    Parameters
    ----------
    non_plasmid_contigs : iterable of NonPlasmidContig
        Iterable of non plasmid contigs
    out_file : Path
        Path to output file
    """
    with out_file.open("w") as f_out:
        for non_plasmid_contig in non_plasmid_contigs:
            f_out.write(non_plasmid_contig.to_tsv_row())
            f_out.write("\n")


def non_plasmid_contigs_from_file(file: Path) -> Iterator[gt_items.NonPlasmidContig]:
    """Get non plasmid contigs from file.

    Parameters
    ----------
    file : Path
        Path to the file

    Yield
    -----
    NonPlasmidContig
        Non plasmid contig
    """
    with file.open() as f_in:
        for line in f_in:
            yield gt_items.NonPlasmidContig.from_tsv_row(line)


class Manager:
    """Ground truth IO manager."""

    MERGED_PLASMID_FASTANAME = Path("merged_plasmid_sequences_fasta")
    MAPPING_SAMNAME = Path("contigs_vs_plamids.sam")
    FILTERED_MAPPING_SAMNAME = Path("filtered_contigs_vs_plamids.sam")

    PLASMID_CONTIGS_FILENAME = Path("plasmid_contigs.tsv")
    NON_PLASMID_CONTIGS_FILENAME = Path("non_plasmid_contigs.tsv")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def merged_plasmid_fasta(self) -> Path:
        """Get merged plasmid FASTA path."""
        return self.__config.output_directory() / self.MERGED_PLASMID_FASTANAME

    def mapping_sam(self) -> Path:
        """Get mapping file path."""
        return self.__config.output_directory() / self.MAPPING_SAMNAME

    def filtered_mapping_sam(self) -> Path:
        """Get filtered mapping file path."""
        return self.__config.output_directory() / self.FILTERED_MAPPING_SAMNAME

    def plasmid_contigs_file(self) -> Path:
        """Get plasmid contigs file path."""
        return self.__config.output_directory() / self.PLASMID_CONTIGS_FILENAME

    def non_plasmid_contigs_file(self) -> Path:
        """Get non plasmid contigs file path."""
        return self.__config.output_directory() / self.NON_PLASMID_CONTIGS_FILENAME

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config:
    """Ground truth IO configuration."""

    DEFAULT_OUTPUT_DIR = Path("./ground_truth")

    KEY_OUTPUT_DIR = "output_directory"

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> Config:
        """Create config instance from a YAML file."""
        with Path(yaml_filepath).open("r") as file:
            config_data = yaml.safe_load(file)
        return cls.from_dict(config_data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_OUTPUT_DIR, cls.DEFAULT_OUTPUT_DIR),
        )

    def __init__(self, output_directory: Path = DEFAULT_OUTPUT_DIR) -> None:
        """Initialize object."""
        self.__output_directory = output_directory

    def output_directory(self) -> Path:
        """Get output directory."""
        return self.__output_directory

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_OUTPUT_DIR: self.__output_directory,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        yaml_filepath = Path(yaml_filepath)
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
