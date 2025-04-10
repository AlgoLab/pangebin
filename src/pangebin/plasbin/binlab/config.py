"""PangeBin-flow binlab config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import typer

from pangebin.yaml_interface import YAMLInterface


class Binlab(YAMLInterface):
    """PangeBin-flow binlab model config class."""

    DEFAULT_GAMMA_MBS = 0.9

    KEY_GAMMA_MBS = "gamma_mbs"

    DEFAULT_YAML_FILE = Path("plasbin_binlab_config.yaml")

    NAME = "PangeBin-flow binlab binning config"

    @classmethod
    def default(cls) -> Binlab:
        """Get default config."""
        return cls(cls.DEFAULT_GAMMA_MBS)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Binlab:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_GAMMA_MBS, cls.DEFAULT_GAMMA_MBS),
        )

    def __init__(
        self,
        gamma_mbs: float,
    ) -> None:
        """Initialize object."""
        self.__gamma_mbs = gamma_mbs

    def gamma_mbs(self) -> float:
        """Get gamma mbs."""
        return self.__gamma_mbs

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_GAMMA_MBS: self.__gamma_mbs,
        }


class BinlabOptions:
    """Plasbin binlab options."""

    _RICH_HELP_PANEL = "Plasbin binlab model options"

    GAMMA_MBS = typer.Option(
        help="Gamma MBS coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--binlab-cfg",
        help="The configuration file path",
        rich_help_panel=_RICH_HELP_PANEL,
    )


if __name__ == "__main__":
    default_config = Binlab.default()
    default_config.to_yaml(Binlab.DEFAULT_YAML_FILE)
