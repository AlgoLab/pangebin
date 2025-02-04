"""Blast mapping operations."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml  # type: ignore[import-untyped]

from pangebin.mapping import items

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


def filter_mappings(
    mappings_df: pd.DataFrame,
    q_len_dict: dict[str, int] | None = None,
    s_len_dict: dict[str, int] | None = None,
    config: Config | None = None,
) -> pd.DataFrame:
    """Return a subset of mappings passing the thresholds.

    Parameters
    ----------
    mappings_df : DataFrame
        DataFrame obtained from read_blast_outfmt6_file
    q_len_dict : dict[str, int], optional
        query id (str): query length (int)
    s_len_dict : dict[str, int], optional
        subject id (str): subject length (int)
    config : Config, optional
        Filter config

    Returns
    -------
    DataFrame
        DataFrame with mappings that pass the thresholds

    """
    if config is None:
        config = Config()
    filtered_mappings_df = pd.DataFrame(columns=mappings_df.columns)
    for _, mapping in mappings_df.iterrows():
        length_ok = mapping[items.BlastColumns.LENGTH] >= config.min_length()
        pident_ok = mapping[items.BlastColumns.PIDENT] >= config.min_pident()

        q_len = (
            abs(
                mapping[items.BlastColumns.QEND] - mapping[items.BlastColumns.QSTART],
            )
            + 1
        )
        q_cov_ok = q_len_dict is None or (
            q_len / q_len_dict[mapping[items.BlastColumns.QSEQID]] >= config.min_q_cov()
        )

        s_len = (
            abs(
                mapping[items.BlastColumns.SEND] - mapping[items.BlastColumns.SSTART],
            )
            + 1
        )
        s_cov_ok = s_len_dict is None or (
            s_len / s_len_dict[mapping[items.BlastColumns.SSEQID]] >= config.min_s_cov()
        )

        if length_ok and pident_ok and q_cov_ok and s_cov_ok:
            filtered_mappings_df.loc[len(filtered_mappings_df)] = mapping

    return filtered_mappings_df


class Config:
    """Standardize config class."""

    DEFAULT_MIN_LENGTH = 1
    DEFAULT_MIN_PIDENT = 0
    DEFAULT_MIN_Q_COV = 0
    DEFAULT_MIN_S_COV = 0

    KEY_MIN_LENGTH = "min_length"
    KEY_MIN_PIDENT = "min_pident"
    KEY_MIN_Q_COV = "min_q_cov"
    KEY_MIN_S_COV = "min_s_cov"

    DEFAULT_YAML_FILE = Path("mapping_filter_config.yaml")

    NAME = "Mapping filter config"

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
            config_dict.get(cls.KEY_MIN_LENGTH, cls.DEFAULT_MIN_LENGTH),
            config_dict.get(cls.KEY_MIN_PIDENT, cls.DEFAULT_MIN_PIDENT),
            config_dict.get(cls.KEY_MIN_Q_COV, cls.DEFAULT_MIN_Q_COV),
            config_dict.get(cls.KEY_MIN_S_COV, cls.DEFAULT_MIN_S_COV),
        )

    def __init__(
        self,
        min_length: int = DEFAULT_MIN_LENGTH,
        min_pident: float = DEFAULT_MIN_PIDENT,
        min_q_cov: float = DEFAULT_MIN_Q_COV,
        min_s_cov: float = DEFAULT_MIN_S_COV,
    ) -> None:
        """Initialize object."""
        self.__min_length = min_length
        self.__min_pident = min_pident
        self.__min_q_cov = min_q_cov
        self.__min_s_cov = min_s_cov

    def min_length(self) -> int:
        """Get min length."""
        return self.__min_length

    def min_pident(self) -> float:
        """Get min pident."""
        return self.__min_pident

    def min_q_cov(self) -> float:
        """Get min q cov."""
        return self.__min_q_cov

    def min_s_cov(self) -> float:
        """Get min s cov."""
        return self.__min_s_cov

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MIN_LENGTH: self.__min_length,
            self.KEY_MIN_PIDENT: self.__min_pident,
            self.KEY_MIN_Q_COV: self.__min_q_cov,
            self.KEY_MIN_S_COV: self.__min_s_cov,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        yaml_filepath = Path(yaml_filepath)
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
