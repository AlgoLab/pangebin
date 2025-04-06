"""PangeBin-flow config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import typer

from pangebin.yaml import YAMLInterface


class Decomp(YAMLInterface):
    """PangeBin-flow decomp model config class."""

    DEFAULT_GAMMA_MCF = 0.9
    DEFAULT_GAMMA_MGC = 0.9

    KEY_GAMMA_MCF = "gamma_mcf"
    KEY_GAMMA_MGC = "gamma_mgc"

    DEFAULT_YAML_FILE = Path("plasbin_decomp_config.yaml")

    NAME = "PangeBin-flow decomp binning config"

    @classmethod
    def default(cls) -> Decomp:
        """Get default config."""
        return cls(cls.DEFAULT_GAMMA_MCF, cls.DEFAULT_GAMMA_MGC)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Decomp:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_GAMMA_MCF, cls.DEFAULT_GAMMA_MCF),
            config_dict.get(cls.KEY_GAMMA_MGC, cls.DEFAULT_GAMMA_MGC),
        )

    def __init__(
        self,
        gamma_mcf: float,
        gamma_mgc: float,
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


class DecompOptions:
    """Plasbin decomp options."""

    _RICH_HELP_PANEL = "Plasbin decomp model options"

    GAMMA_MCF = typer.Option(
        help="Gamma MCF coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    GAMMA_MGC = typer.Option(
        help="Gamma MGC coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--decomp-cfg",
        help="The configuration file path",
        rich_help_panel=_RICH_HELP_PANEL,
    )


if __name__ == "__main__":
    default_config = Decomp.default()
    default_config.to_yaml(Decomp.DEFAULT_YAML_FILE)
