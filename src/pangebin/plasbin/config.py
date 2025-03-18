"""PangeBin-flow config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


class Config:
    """PangeBin-flow config class."""

    DEFAULT_GAMMA_MCF = 0.9
    DEFAULT_GAMMA_MGC = 0.9

    KEY_GAMMA_MCF = "gamma_mcf"
    KEY_GAMMA_MGC = "gamma_mgc"

    DEFAULT_YAML_FILE = Path("pangebin_flow_config.yaml")

    NAME = "PangeBin-flow config"

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
            config_dict.get(cls.KEY_GAMMA_MCF, cls.DEFAULT_GAMMA_MCF),
            config_dict.get(cls.KEY_GAMMA_MGC, cls.DEFAULT_GAMMA_MGC),
        )

    def __init__(
        self,
        gamma_mcf: float = DEFAULT_GAMMA_MCF,
        gamma_mgc: float = DEFAULT_GAMMA_MGC,
    ) -> None:
        """Initialize object."""
        self.__gamma_mcf = gamma_mcf
        self.__gamma_mgc = gamma_mgc

    def gamma_mcf(self) -> float:
        """Get gamma mcf."""
        return self.__gamma_mcf

    def gamma_mgc(self) -> float:
        """Get gamma mgc."""
        return self.__gamma_mgc

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_GAMMA_MCF: self.__gamma_mcf,
            self.KEY_GAMMA_MGC: self.__gamma_mgc,
        }

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
