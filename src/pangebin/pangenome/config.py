"""Configuration for pangenome."""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Any

from pangebin.yaml_interface import YAMLInterface

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


class Pangenome(YAMLInterface):
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


def fmt_nfcore_pan_cli(
    key: str,
    value: Any,  # noqa: ANN401
) -> str:
    """Format nf-core/pangenome CLI option."""
    return f"--{key} {value}"
