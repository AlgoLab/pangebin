"""NCBI Entrez module."""

from __future__ import annotations

import logging

import typer
from Bio import Entrez

from pangebin.yaml_interface import YAMLInterface

_LOGGER = logging.getLogger(__name__)


class Config(YAMLInterface):
    """NCBI Entrez configurations."""

    KEY_EMAIL = "email"
    KEY_TOOL = "tool"
    KEY_API_KEY = "api_key"
    KEY_MAX_TRIES = "max_tries"
    KEY_SLEEP_BETWEEN_TRIES = "sleep_between_tries"

    DEFAULT_EMAIL = None
    DEFAULT_TOOL = "biopython"
    DEFAULT_API_KEY = None
    DEFAULT_MAX_TRIES = 3
    DEFAULT_SLEEP_BETWEEN_TRIES = 15

    @classmethod
    def from_dict(cls, config: dict) -> Config:
        """Create Entrez object from dictionary."""
        return cls(
            email=config.get(cls.KEY_EMAIL, cls.DEFAULT_EMAIL),
            tool=config.get(cls.KEY_TOOL, cls.DEFAULT_TOOL),
            api_key=config.get(cls.KEY_API_KEY, cls.DEFAULT_API_KEY),
            max_tries=config.get(cls.KEY_MAX_TRIES, cls.DEFAULT_MAX_TRIES),
            sleep_between_tries=config.get(
                cls.KEY_SLEEP_BETWEEN_TRIES,
                cls.DEFAULT_SLEEP_BETWEEN_TRIES,
            ),
        )

    def __init__(
        self,
        email: str | None = DEFAULT_EMAIL,
        tool: str = DEFAULT_TOOL,
        api_key: str | None = DEFAULT_API_KEY,
        max_tries: int = DEFAULT_MAX_TRIES,
        sleep_between_tries: int = DEFAULT_SLEEP_BETWEEN_TRIES,
    ) -> None:
        """Initialize object."""
        self.__email = email
        self.__tool = tool
        self.__api_key = api_key
        self.__max_tries = max_tries
        self.__sleep_between_tries = sleep_between_tries

    def email(self) -> str | None:
        """Return email."""
        return self.__email

    def tool(self) -> str:
        """Return tool."""
        return self.__tool

    def api_key(self) -> str | None:
        """Return API key."""
        return self.__api_key

    def max_tries(self) -> int:
        """Return max tries."""
        return self.__max_tries

    def sleep_between_tries(self) -> int:
        """Return sleep between tries."""
        return self.__sleep_between_tries

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_EMAIL: self.email(),
            self.KEY_TOOL: self.tool(),
            self.KEY_API_KEY: self.api_key(),
            self.KEY_MAX_TRIES: self.max_tries(),
            self.KEY_SLEEP_BETWEEN_TRIES: self.sleep_between_tries(),
        }


def set_entrez_config(config: Config) -> None:
    """Set Entrez configuration."""
    _LOGGER.debug("Setting Entrez configuration")
    if config.email() is not None:
        Entrez.email = config.email()  # type: ignore[assignment]
    Entrez.tool = config.tool()
    if config.api_key() is not None:
        Entrez.api_key = config.api_key()  # type: ignore[assignment]
    Entrez.max_tries = config.max_tries()
    Entrez.sleep_between_tries = config.sleep_between_tries()


class AppOptions:
    """Entrez app options."""

    __RICH_HELP_PANEL = "Entrez options"

    EMAIL = typer.Option(
        "--email",
        help="Email address to fetch NCBI database",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    TOOL = typer.Option(
        "--tool",
        help="Tool name",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    API_KEY = typer.Option(
        "--api-key",
        help="API key",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    MAX_TRIES = typer.Option(
        "--max-tries",
        help="Max tries",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    SLEEP_BETWEEN_TRIES = typer.Option(
        "--sleep-between-tries",
        help="Sleep between tries",
        rich_help_panel=__RICH_HELP_PANEL,
    )

    CONFIG_FILE = typer.Option(
        "--entrez-cfg",
        help="The configuration file path",
        rich_help_panel=__RICH_HELP_PANEL,
    )
