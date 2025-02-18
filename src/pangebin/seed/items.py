"""Seed sequence items."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class SeedContigThresholdTestItem:
    """Seed contig threshold test item."""

    @classmethod
    def from_datatest_line(cls, line: str) -> SeedContigThresholdTestItem:
        """Create object from datatest line."""
        items_str: list[str] = line.split()
        return cls(
            Path(items_str[0]),
            Path(items_str[1]),
            Path(items_str[2]),
        )

    def __init__(
        self,
        plasmid_contigs_file: Path,
        non_plasmid_contigs_file: Path,
        contig_gene_densities_file: Path,
    ) -> None:
        self.__plasmid_contigs_file = plasmid_contigs_file
        self.__non_plasmid_contigs_file = non_plasmid_contigs_file
        self.__contig_gene_densities_file = contig_gene_densities_file

    def plasmid_contigs_file(self) -> Path:
        """Plasmid contigs file."""
        return self.__plasmid_contigs_file

    def non_plasmid_contigs_file(self) -> Path:
        """Non plasmid contigs file."""
        return self.__non_plasmid_contigs_file

    def contig_gene_densities_file(self) -> Path:
        """Contig gene densities file."""
        return self.__contig_gene_densities_file

    def to_dataset_line(self) -> str:
        """Convert to dataset line."""
        return (
            f"{self.__plasmid_contigs_file}\t"
            f"{self.__non_plasmid_contigs_file}\t"
            f"{self.__contig_gene_densities_file}"
        )


class SeedContigThresholds:
    """Seed contig thresholds."""

    KEY_MEANS = "means"
    KEY_BEST_SP_NPS = "best_sp_nps"
    KEY_THRESHOLDS = "thresholds"

    KEY_LENGHT = "length"
    KEY_GENE_DENSITY = "gene_density"

    @classmethod
    def from_yaml(cls, file: Path) -> SeedContigThresholds:
        """Create object from yaml file."""
        with file.open() as f_in:
            dict_from_yaml = yaml.safe_load(f_in)
        return cls.from_dict(dict_from_yaml)

    @classmethod
    def from_dict(cls, dict_from_yaml: dict) -> SeedContigThresholds:
        """Create object from dict."""
        return cls(
            dict_from_yaml[cls.KEY_MEANS][cls.KEY_LENGHT],
            dict_from_yaml[cls.KEY_MEANS][cls.KEY_GENE_DENSITY],
            dict_from_yaml[cls.KEY_THRESHOLDS],
            (
                (
                    int(threshold_pair[cls.KEY_LENGHT]),
                    float(threshold_pair[cls.KEY_GENE_DENSITY]),
                )
                for threshold_pair in dict_from_yaml[cls.KEY_THRESHOLDS]
            ),
        )

    def __init__(
        self,
        mean_length: float,
        mean_gene_density: float,
        best_sp_nps: int,
        threshold_pairs: Iterable[tuple[int, float]],
    ) -> None:
        self.__mean_length = mean_length
        self.__mean_gene_density = mean_gene_density
        self.__best_sp_nps = best_sp_nps
        self.__threshold_pairs = list(threshold_pairs)

    def mean_length(self) -> float:
        """Mean sequence length."""
        return self.__mean_length

    def mean_gene_density(self) -> float:
        """Mean gene density."""
        return self.__mean_gene_density

    def best_sp_nps(self) -> int:
        """Best SP - NPS."""
        return self.__best_sp_nps

    def threshold_pairs(self) -> list[tuple[int, float]]:
        """Threshold pairs."""
        return self.__threshold_pairs

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MEANS: {
                self.KEY_LENGHT: self.__mean_length,
                self.KEY_GENE_DENSITY: self.__mean_gene_density,
            },
            self.KEY_BEST_SP_NPS: self.__best_sp_nps,
            self.KEY_THRESHOLDS: [
                {
                    self.KEY_LENGHT: length,
                    self.KEY_GENE_DENSITY: gene_density,
                }
                for length, gene_density in self.__threshold_pairs
            ],
        }

    def to_yaml(self, file: Path) -> Path:
        """Convert to yaml file."""
        with file.open("w") as f_out:
            yaml.dump(self.to_dict(), f_out, Dumper, sort_keys=False)
        return file
