"""Input/output module for pangenome."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pangebin.yaml_interface import YAMLInterface


class Manager:
    """Standardize input/output manager."""

    PANGENOME = Path("pangenome.gfa")

    NFCORE_SUBDIR = Path("nfcore")
    MIXED_FASTA_GFA_EXTENSIONS = ".gz.gfaffix.unchop.Ygs.view.gfa"

    FINAL_GFA_DIR = Path("FINAL_GFA")

    def __init__(self, config: Config) -> None:
        """Initialize object."""
        self.__config = config

    def pangenome_gfa_path(self) -> Path:
        """Get Pangenome GFA path."""
        return self.__config.output_directory() / self.PANGENOME

    def nfcore_pangenome_directory(self) -> Path:
        """Get Nextflow core Pangenome directory path."""
        return self.config().output_directory() / self.NFCORE_SUBDIR

    def nfcore_pangenome_gfa_path(self) -> Path:
        """Get Nextflow core Pangenome GFA path."""
        return (
            self.nfcore_pangenome_directory()
            / self.FINAL_GFA_DIR
            / (self.__config.mixed_fasta_path().name + self.MIXED_FASTA_GFA_EXTENSIONS)
        )

    def config(self) -> Config:
        """Get config."""
        return self.__config


class Config(YAMLInterface):
    """Pangenome IO config class."""

    DEFAULT_DIR = Path("./pangenome")

    KEY_MIXED_FASTA_PATH = "mixed_fasta_path"
    KEY_OUTPUT_DIR = "output_directory"

    NAME = "Pangenome IO config"

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            Path(config_dict[cls.KEY_MIXED_FASTA_PATH]),
            config_dict.get(cls.KEY_OUTPUT_DIR, cls.DEFAULT_DIR),
        )

    def __init__(
        self,
        mixed_fasta_path: Path,
        output_directory: Path = DEFAULT_DIR,
    ) -> None:
        """Initialize object."""
        self.__mixed_fasta_path = mixed_fasta_path
        self.__output_directory = output_directory

    def mixed_fasta_path(self) -> Path:
        """Get mixed FASTA path."""
        return self.__mixed_fasta_path

    def output_directory(self) -> Path:
        """Get output directory."""
        return self.__output_directory

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MIXED_FASTA_PATH: self.__mixed_fasta_path,
            self.KEY_OUTPUT_DIR: self.__output_directory,
        }
