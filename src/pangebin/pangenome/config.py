"""Configuration for pangenome."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from enum import StrEnum

_PARENT_DIR = Path(__file__).resolve().parent

BASE_NFCORE_PANGENOME_CONFIG_PATH = _PARENT_DIR / "nf-core_pangenome_config.yaml"

SEQWISH_MIN_MATCH_LENGTH = 1000
SMOOTHXG_POA_PARAMS = "asm5"
WFMASH_SEGMENT_LENGTH = 200
WFMASH_CHUNKS = 1
WFMASH_MAP_PCT_ID = 95.0


class NextflowProfile(StrEnum):
    """Supported Nextflow profiles."""

    TEST = "test"
    DOCKER = "docker"
    SINGULARITY = "singularity"
    PODMAN = "podman"
    SHIFTER = "shifter"
    CHARLIECLOUD = "charliecloud"
    APPTAINER = "apptainer"
    CONDA = "conda"


class Pangenome:
    """Configuration for pangenome."""

    DEFAULT_RELEASE = "1.1.2"
    DEFAULT_PROFILE = NextflowProfile.PODMAN
    DEFAULT_RESUME = True
    DEFAULT_SUPP_NFCORE_PANGENOME_CONFIG_PATH = None

    KEY_RELEASE = "release"
    KEY_PROFILE = "profile"
    KEY_RESUME = "resume"
    KEY_NFCORE_PANGENOME_CONFIG_PATH = "supp_nfcore_pangenome_config_path"

    DEFAULT_YAML_FILE = Path("pangenome_config.yaml")

    NAME = "Pangenome config"

    @classmethod
    def from_yaml(cls, yaml_filepath: Path) -> Pangenome:
        """Create config instance from a YAML file."""
        with Path(yaml_filepath).open("r") as file:
            config_data = yaml.safe_load(file)
        return cls.from_dict(config_data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> Pangenome:
        """Convert dict to object."""
        return cls(
            config_dict.get(cls.KEY_RELEASE, cls.DEFAULT_RELEASE),
            NextflowProfile(config_dict.get(cls.KEY_PROFILE, cls.DEFAULT_PROFILE)),
            config_dict.get(cls.KEY_RESUME, cls.DEFAULT_RESUME),
            Path(
                config_dict.get(
                    cls.KEY_NFCORE_PANGENOME_CONFIG_PATH,
                    cls.DEFAULT_SUPP_NFCORE_PANGENOME_CONFIG_PATH,
                ),
            ),
        )

    def __init__(
        self,
        release: str = DEFAULT_RELEASE,
        profile: NextflowProfile = DEFAULT_PROFILE,
        resume: bool = DEFAULT_RESUME,  # noqa: FBT001
        supp_nfcore_pangenome_config_path: Path
        | None = DEFAULT_SUPP_NFCORE_PANGENOME_CONFIG_PATH,
    ) -> None:
        """Initialize object."""
        self.__release = release
        self.__profile = profile
        self.__resume = resume
        self.__supp_nfcore_pangenome_config_path = supp_nfcore_pangenome_config_path

    def release(self) -> str:
        """Get release."""
        return self.__release

    def resume(self) -> bool:
        """Get resume."""
        return self.__resume

    def profile(self) -> NextflowProfile:
        """Get profile."""
        return self.__profile

    def supp_nfcore_pangenome_config_path(self) -> Path | None:
        """Get nf-core/pangenome config path."""
        return self.__supp_nfcore_pangenome_config_path

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict."""
        cfg_dict = {
            self.KEY_RELEASE: self.__release,
            self.KEY_PROFILE: str(self.__profile),
            self.KEY_RESUME: self.__resume,
        }
        if self.__supp_nfcore_pangenome_config_path is not None:
            cfg_dict[self.KEY_NFCORE_PANGENOME_CONFIG_PATH] = str(
                self.__supp_nfcore_pangenome_config_path,
            )
        return cfg_dict

    def to_yaml(self, yaml_filepath: Path) -> Path:
        """Write to yaml."""
        yaml_filepath = Path(yaml_filepath)
        with yaml_filepath.open("w") as file:
            yaml.dump(self.to_dict(), file, Dumper=Dumper, sort_keys=False)
        return yaml_filepath


def fmt_nfcore_pan_cli(
    key: str,
    value: Any,  # noqa: ANN401
) -> str:
    """Format nf-core/pangenome CLI option."""
    return f"--{key} {value}"
