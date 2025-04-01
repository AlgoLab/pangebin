"""PangeBin-flow config module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import typer

import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.network as net
from pangebin.yaml import YAMLInterface


class Binning(YAMLInterface):
    """PangeBin-flow binning config class."""

    DEFAULT_SINK_ARCS_DOMAIN = net.SinkArcsDomain.ALL
    DEFAULT_MIN_FLOW = 0.0001
    DEFAULT_MIN_CUMULATIVE_LENGTH = 1000
    DEFAULT_CIRCULAR = False
    DEFAULT_OBJ_FUN_DOMAIN = pb_lp_obj.ObjectiveFunctionDomain.ALL
    DEFAULT_GAMMA_MCF = 0.9
    DEFAULT_GAMMA_MGC = 0.9

    KEY_SINK_ARC_DEFINITION = "sink_arc_definition"
    KEY_MIN_FLOW = "min_flow"
    KEY_MIN_CUMULATIVE_LENGTH = "min_cumulative_len"
    KEY_CIRCULAR = "circular"
    KEY_OBJ_FUN_DOMAIN = "obj_fun_domain"
    KEY_GAMMA_MCF = "gamma_mcf"
    KEY_GAMMA_MGC = "gamma_mgc"

    DEFAULT_YAML_FILE = Path("binning_config.yaml")

    NAME = "PangeBin-flow binning config"

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Binning:
        """Convert dict to object."""
        return cls(
            config_dict.get(
                cls.KEY_SINK_ARC_DEFINITION,
                cls.DEFAULT_SINK_ARCS_DOMAIN,
            ),
            config_dict.get(cls.KEY_MIN_FLOW, cls.DEFAULT_MIN_FLOW),
            config_dict.get(
                cls.KEY_MIN_CUMULATIVE_LENGTH,
                cls.DEFAULT_MIN_CUMULATIVE_LENGTH,
            ),
            config_dict.get(cls.KEY_CIRCULAR, cls.DEFAULT_CIRCULAR),
            config_dict.get(cls.KEY_OBJ_FUN_DOMAIN, cls.DEFAULT_OBJ_FUN_DOMAIN),
            config_dict.get(cls.KEY_GAMMA_MCF, cls.DEFAULT_GAMMA_MCF),
            config_dict.get(cls.KEY_GAMMA_MGC, cls.DEFAULT_GAMMA_MGC),
        )

    def __init__(  # noqa: PLR0913
        self,
        sink_arcs_domain: net.SinkArcsDomain = DEFAULT_SINK_ARCS_DOMAIN,
        min_flow: float = DEFAULT_MIN_FLOW,
        min_cumulative_len: int = DEFAULT_MIN_CUMULATIVE_LENGTH,
        circular: bool = DEFAULT_CIRCULAR,  # noqa: FBT001
        obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain = DEFAULT_OBJ_FUN_DOMAIN,
        gamma_mcf: float = DEFAULT_GAMMA_MCF,
        gamma_mgc: float = DEFAULT_GAMMA_MGC,
    ) -> None:
        """Initialize object."""
        self.__sink_arcs_domain = sink_arcs_domain
        self.__min_flow = min_flow
        self.__min_cumulative_len = min_cumulative_len
        self.__circular = circular
        self.__obj_fun_domain = obj_fun_domain
        self.__gamma_mcf = gamma_mcf
        self.__gamma_mgc = gamma_mgc

    def sink_arcs_domain(self) -> net.SinkArcsDomain:
        """Get sink-arcs domain."""
        return self.__sink_arcs_domain

    def min_flow(self) -> float:
        """Get min flow."""
        return self.__min_flow

    def min_cumulative_len(self) -> int:
        """Get min cumulative length."""
        return self.__min_cumulative_len

    def circular(self) -> bool:
        """Get circular."""
        return self.__circular

    def obj_fun_domain(self) -> pb_lp_obj.ObjectiveFunctionDomain:
        """Get objective function domain."""
        return self.__obj_fun_domain

    def gamma_mcf(self) -> float:
        """Get gamma mcf."""
        return self.__gamma_mcf

    def gamma_mgc(self) -> float:
        """Get gamma mgc."""
        return self.__gamma_mgc

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_SINK_ARC_DEFINITION: str(self.__sink_arcs_domain),
            self.KEY_MIN_FLOW: self.__min_flow,
            self.KEY_MIN_CUMULATIVE_LENGTH: self.__min_cumulative_len,
            self.KEY_CIRCULAR: self.__circular,
            self.KEY_OBJ_FUN_DOMAIN: str(self.__obj_fun_domain),
            self.KEY_GAMMA_MCF: self.__gamma_mcf,
            self.KEY_GAMMA_MGC: self.__gamma_mgc,
        }


class BinningOptions:
    """Binning options."""

    _RICH_HELP_PANEL = "Binning options"

    SINK_ARCS_DOMAIN = typer.Option(
        help="Sink-arcs domain",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    MIN_FLOW = typer.Option(
        help="Minimum flow",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    MIN_CUMULATIVE_LENGTH = typer.Option(
        help="Minimum cumulative length",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    CIRCULAR = typer.Option(
        help="The flow is circular",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    OBJ_FUN_DOMAIN = typer.Option(
        help="Objective function domain",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    GAMMA_MCF = typer.Option(
        help="Gamma MCF coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    GAMMA_MGC = typer.Option(
        help="Gamma MGC coefficient",
        rich_help_panel=_RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--bin-cfg",
        help="The configuration file path",
        rich_help_panel=_RICH_HELP_PANEL,
    )


if __name__ == "__main__":
    default_config = Binning()
    default_config.to_yaml(Binning.DEFAULT_YAML_FILE)
