"""YAML interface module."""

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from pathlib import Path

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from typing import Any, Self


class YAMLInterface(metaclass=ABCMeta):
    """YAML interface."""

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> Self:
        """Get object from YAML file."""
        with yaml_filepath.open("r") as file:
            obj_dict = yaml.safe_load(file)
        return cls.from_dict(obj_dict)

    @classmethod
    @abstractmethod
    def from_dict(cls, obj_dict: dict[str, Any]) -> Self:
        """Convert dict to object."""
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        raise NotImplementedError

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        yaml_filepath = Path(yaml_filepath)
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath
