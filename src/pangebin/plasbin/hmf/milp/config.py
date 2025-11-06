"""HMF MILP config module."""

from __future__ import annotations

from typing import Any

import pangebin.plasbin.milp.objectives as cmn_lp_objs
from pangebin.yaml_interface import YAMLInterface


# REFACTOR split config API from IO reader/writer
# Config
# ConfigReader
# ConfigWriter
class Config(YAMLInterface):
    """HMF MILP configurations."""

    DEFAULT_MIN_FLOW = 0.0001
    DEFAULT_MIN_CUMULATIVE_LENGTH = 1000
    DEFAULT_PLASMIDNESS_COEFFICIENT = 0
    DEFAULT_OBJ_FUN_DOMAIN = cmn_lp_objs.ObjectiveFunctionDomain.ALL

    KEY_MIN_FLOW = "min_flow"
    KEY_MIN_CUMULATIVE_LENGTH = "min_cumulative_len"
    KEY_PLASMIDNESS_COEFFICIENT = "plasmidness_coefficient"
    KEY_OBJ_FUN_DOMAIN = "obj_fun_domain"

    @classmethod
    def default(cls) -> Config:
        """Get default config."""
        return cls(
            cls.DEFAULT_MIN_FLOW,
            cls.DEFAULT_MIN_CUMULATIVE_LENGTH,
            cls.DEFAULT_PLASMIDNESS_COEFFICIENT,
            cls.DEFAULT_OBJ_FUN_DOMAIN,
        )

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Config:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_MIN_FLOW, cls.DEFAULT_MIN_FLOW),
            config_dict.get(
                cls.KEY_MIN_CUMULATIVE_LENGTH,
                cls.DEFAULT_MIN_CUMULATIVE_LENGTH,
            ),
            config_dict.get(
                cls.KEY_PLASMIDNESS_COEFFICIENT,
                cls.DEFAULT_PLASMIDNESS_COEFFICIENT,
            ),
            cmn_lp_objs.ObjectiveFunctionDomain(
                config_dict.get(cls.KEY_OBJ_FUN_DOMAIN, cls.DEFAULT_OBJ_FUN_DOMAIN),
            ),
        )

    def __init__(
        self,
        min_flow: float,
        min_cumulative_len: int,
        plasmidness_coefficient: float,
        obj_fun_domain: cmn_lp_objs.ObjectiveFunctionDomain,
    ) -> None:
        """Initialize object."""
        self.__min_flow = min_flow
        self.__min_cumulative_len = min_cumulative_len
        self.__plasmidness_coefficient = plasmidness_coefficient
        self.__obj_fun_domain = obj_fun_domain

    def min_flow(self) -> float:
        """Get min flow."""
        return self.__min_flow

    def min_cumulative_len(self) -> int:
        """Get min cumulative length."""
        return self.__min_cumulative_len

    def plasmidness_coefficient(self) -> float:
        """Get plasmidness coefficient."""
        return self.__plasmidness_coefficient

    def obj_fun_domain(self) -> cmn_lp_objs.ObjectiveFunctionDomain:
        """Get objective function domain."""
        return self.__obj_fun_domain

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        return {
            self.KEY_MIN_FLOW: self.__min_flow,
            self.KEY_MIN_CUMULATIVE_LENGTH: self.__min_cumulative_len,
            self.KEY_PLASMIDNESS_COEFFICIENT: self.__plasmidness_coefficient,
            self.KEY_OBJ_FUN_DOMAIN: str(self.__obj_fun_domain),
        }
