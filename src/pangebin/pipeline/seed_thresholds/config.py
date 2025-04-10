"""Seed thresholds pipeline configuration module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pangebin.entrez as pg_entrez
import pangebin.ground_truth.config as gt_cfg
import pangebin.mapping.filter as map_filter
import pangebin.seed.thresholds.config as seed_thr_cfg
from pangebin import subprocess_lib
from pangebin.yaml_interface import YAMLInterface


class Config(YAMLInterface):
    """Seed thresholds pipeline configuration."""

    KEY_GROUND_TRUTH_CFG = "ground_truth_cfg"
    KEY_SAM_FILTER_CFG = "sam_filter_cfg"
    KEY_THR_RANGES_CFG = "threshold_ranges_cfg"
    KEY_ENTREZ_CFG = "entrez_cfg"
    KEY_RESSOURCES_CFG = "ressources_cfg"

    DEFAULT_FILENAME = Path("pipeline_seed_thresholds_cfg.yaml")

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            gt_cfg.Config.from_dict(config_dict[cls.KEY_GROUND_TRUTH_CFG])
            if cls.KEY_GROUND_TRUTH_CFG in config_dict
            else gt_cfg.Config(),
            map_filter.Config.from_dict(config_dict[cls.KEY_SAM_FILTER_CFG])
            if cls.KEY_SAM_FILTER_CFG in config_dict
            else map_filter.Config(),
            seed_thr_cfg.ThresholdRanges.from_dict(config_dict[cls.KEY_THR_RANGES_CFG])
            if cls.KEY_THR_RANGES_CFG in config_dict
            else seed_thr_cfg.ThresholdRanges(),
            pg_entrez.Config.from_dict(config_dict[cls.KEY_ENTREZ_CFG])
            if cls.KEY_ENTREZ_CFG in config_dict
            else pg_entrez.Config(),
            subprocess_lib.RessourcesConfig.from_dict(
                config_dict[cls.KEY_RESSOURCES_CFG],
            )
            if cls.KEY_RESSOURCES_CFG in config_dict
            else subprocess_lib.RessourcesConfig(),
        )

    def __init__(
        self,
        ground_truth_config: gt_cfg.Config | None = None,
        sam_filter_config: map_filter.Config | None = None,
        threshold_ranges: seed_thr_cfg.ThresholdRanges | None = None,
        entrez_config: pg_entrez.Config | None = None,
        ressources_config: subprocess_lib.RessourcesConfig | None = None,
    ) -> None:
        self.__gt_cfg = (
            ground_truth_config if ground_truth_config is not None else gt_cfg.Config()
        )
        self.__sam_filter_cgf = (
            sam_filter_config if sam_filter_config is not None else map_filter.Config()
        )
        self.__threshold_ranges = (
            threshold_ranges
            if threshold_ranges is not None
            else seed_thr_cfg.ThresholdRanges()
        )
        self.__entrez_cfg = (
            entrez_config if entrez_config is not None else pg_entrez.Config()
        )
        self.__ressource_cfg = (
            ressources_config
            if ressources_config is not None
            else subprocess_lib.RessourcesConfig()
        )

    def ground_truth_config(self) -> gt_cfg.Config:
        """Get ground truth config."""
        return self.__gt_cfg

    def sam_filter_config(self) -> map_filter.Config:
        """Get sam filter config."""
        return self.__sam_filter_cgf

    def threshold_ranges(self) -> seed_thr_cfg.ThresholdRanges:
        """Get threshold ranges."""
        return self.__threshold_ranges

    def entrez_config(self) -> pg_entrez.Config:
        """Get entrez config."""
        return self.__entrez_cfg

    def ressources_config(self) -> subprocess_lib.RessourcesConfig:
        """Get ressources config."""
        return self.__ressource_cfg

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_GROUND_TRUTH_CFG: self.__gt_cfg.to_dict(),
            self.KEY_SAM_FILTER_CFG: self.__sam_filter_cgf.to_dict(),
            self.KEY_THR_RANGES_CFG: self.__threshold_ranges.to_dict(),
            self.KEY_ENTREZ_CFG: self.__entrez_cfg.to_dict(),
            self.KEY_RESSOURCES_CFG: self.__ressource_cfg.to_dict(),
        }


if __name__ == "__main__":
    Config().to_yaml(Config.DEFAULT_FILENAME)
