"""Script to help to publish a new release of the python package.

Usage
-----
```sh
python3.10 config/new_release.py --help
```

Warning:
-------
This script is user-interactive.

"""

# Due to typer usage:
# ruff: noqa: TC001, TC003, UP007, FBT001, FBT002, PLR0913

from __future__ import annotations

import datetime
import logging
import shlex
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Annotated

import tomli
import typer
from packaging.version import InvalidVersion, Version

APP = typer.Typer(rich_markup_mode="rich")

_LOGGER = logging.getLogger("new_release")
_LOGGER.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
formatter = logging.Formatter(
    "%(asctime)s - %(name)s [%(levelname)s] %(message)s",
)
console_handler.setFormatter(formatter)
_LOGGER.addHandler(console_handler)


PYPROJECT: Path = Path("pyproject.toml")

DEFAULT_EXECUTE = False
#
# Git
#
DEFAULT_MAIN_BRANCH = "master"
DEFAULT_DEV_BRANCH = "develop"


class ChangesLevel(str, Enum):
    """Changes level."""

    PATCH = "patch"
    MINOR = "minor"
    MAJOR = "major"


@dataclass
class HelpReleaseArgs:
    """Help to release a new package version."""

    ARG_LEVEL = typer.Argument(
        help="The changes level",
    )

    OPT_EXECUTE = typer.Option(
        help="Execute the script",
    )


@APP.command()
def help_release(
    level: Annotated[ChangesLevel, HelpReleaseArgs.ARG_LEVEL],
    execute: Annotated[bool, HelpReleaseArgs.OPT_EXECUTE] = DEFAULT_EXECUTE,
) -> None:
    """Help to release a new package version."""
    current_version = get_current_version()
    new_version = get_new_version(current_version, level)
    new_release_branch(new_version)
    make_release_changes(new_version)
    commit_release_changes(new_version)
    merge_release_main(new_version)
    merge_release_develop(new_version)
    remove_release_branch(new_version)


def get_current_version() -> Version:
    """Get current version from toml.

    Returns
    -------
    Version
        Current version number

    """
    with PYPROJECT.open(mode="rb") as f_in:
        config = tomli.load(f_in)
    version_str: str = config["project"]["version"]
    try:
        current_version: Version = Version(version_str)
    except InvalidVersion:
        _LOGGER.critical("The current version `%s` is not valid.", version_str)
        sys.exit(1)
    _LOGGER.info("The current version is: `%s`", current_version)
    return current_version


def get_new_version(current_version: Version, changes_level: ChangesLevel) -> Version:
    """Get the new version number.

    Parameters
    ----------
    current_version : Version
        Current version number to compare with
    changes_level : ChangesLevel
        The changes level

    Returns
    -------
    Version
        New correct version number

    """
    match changes_level:
        case ChangesLevel.PATCH:
            return Version(
                f"{current_version.major}"
                f".{current_version.minor}"
                f".{current_version.micro + 1}",
            )
        case ChangesLevel.MINOR:
            return Version(
                f"{current_version.major}.{current_version.minor + 1}.0",
            )
        case ChangesLevel.MAJOR:
            return Version(
                f"{current_version.major + 1}.0.0",
            )


def new_release_branch(new_version: Version) -> None:
    """Return the new release branch message.

    Parameters
    ----------
    new_version : Version
        Version number

    """
    cmd = [
        "git",
        "checkout",
        "-b",
        f"release-{new_version}",
        f"{DEFAULT_DEV_BRANCH}",
    ]
    _LOGGER.info("Creating a new release branch:\n> %s\n", " ".join(cmd))
    try:
        subprocess.run(  # noqa: S603
            cmd,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while creating a new release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)


def make_release_changes(new_version: Version) -> None:
    """Make the release changes message.

    Parameters
    ----------
    new_version : Version
        Version number

    """
    __wait_done_request(
        f"[REQUEST] Changes for the release:\n"
        "* Change the version number in the `pyproject.toml` file"
        f" to {new_version}\n"
        "* Update the content of the `CHANGELOG.md` file\n"
        "    * Release version and date"
        f" `[{new_version}] -"
        f" {datetime.datetime.now(tz=datetime.timezone.utc).date()}`\n"
        "    * Please be carefull of the content\n"
        "* Then enter `y`\n",
    )


def commit_release_changes(new_version: Version) -> None:
    """Commit the release changes message.

    Parameters
    ----------
    new_version : Version
        Version number

    """
    cmd = [
        "git",
        "commit",
        "-a",
        "-m",
        f"Bumped version number to {new_version}",
    ]
    _LOGGER.info(
        "Commiting the release changes:\n> %s\n",
        " ".join(shlex.quote(s) for s in cmd),
    )
    try:
        subprocess.run(  # noqa: S603
            cmd,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while commiting the release changes.")
        _LOGGER.critical(exc)
        sys.exit(1)


def merge_release_main(new_version: Version) -> None:
    """Merge the release and the main branches message.

    Parameters
    ----------
    new_version : Version
        Version number

    """
    cmd_checkout = [
        "git",
        "checkout",
        f"{DEFAULT_MAIN_BRANCH}",
    ]
    cmd_merge = [
        "git",
        "merge",
        "--no-ff",
        f"release-{new_version}",
        "-m",
        "Merge release",
    ]
    cmd_tag = [
        "git",
        "tag",
        "-a",
        f"v{new_version}",
        "-m",
        f"Release v{new_version}",
    ]
    cmd_push = [
        "git",
        "push",
        "origin",
        f"{DEFAULT_MAIN_BRANCH}",
        "--follow-tags",
    ]
    _LOGGER.info(
        "Merging the release branch to the main branch:\n> %s\n> %s\n> %s\n> %s\n",
        " ".join(cmd_checkout),
        " ".join(shlex.quote(s) for s in cmd_merge),
        " ".join(shlex.quote(s) for s in cmd_tag),
        " ".join(cmd_push),
    )
    try:
        subprocess.run(  # noqa: S603
            cmd_checkout,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while checkouting the main branch.")
        _LOGGER.critical(exc)
        sys.exit(1)

    try:
        subprocess.run(  # noqa: S603
            cmd_merge,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while merging the release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)
    try:
        subprocess.run(  # noqa: S603
            cmd_tag,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while tagging the release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)

    try:
        subprocess.run(  # noqa: S603
            cmd_push,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while pushing the release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)


def merge_release_develop(new_version: Version) -> None:
    """Merge the release and the develop branches message.

    Parameters
    ----------
    new_version : Version
        Version number

    """
    cmd_checkout = [
        "git",
        "checkout",
        f"{DEFAULT_DEV_BRANCH}",
    ]
    cmd_merge = [
        "git",
        "merge",
        "--no-ff",
        f"release-{new_version}",
        "-m",
        "Merge release",
    ]
    cmd_push = [
        "git",
        "push",
        "origin",
        f"{DEFAULT_DEV_BRANCH}",
        "--follow-tags",
    ]
    _LOGGER.info(
        "Merging the release branch to the develop branch:\n> %s\n> %s\n> %s\n",
        " ".join(cmd_checkout),
        " ".join(shlex.quote(s) for s in cmd_merge),
        " ".join(cmd_push),
    )
    try:
        subprocess.run(  # noqa: S603
            cmd_checkout,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while checkouting the develop branch.")
        _LOGGER.critical(exc)
        sys.exit(1)
    try:
        subprocess.run(  # noqa: S603
            cmd_merge,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while merging the release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)
    try:
        subprocess.run(  # noqa: S603
            cmd_push,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while pushing the release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)


def remove_release_branch(new_version: Version) -> None:
    """Remove the release branch message.

    Parameters
    ----------
    new_version : Version
        Version number

    """
    cmd = [
        "git",
        "branch",
        "-d",
        f"release-{new_version}",
    ]
    _LOGGER.info("Removing the release branch:\n> %s\n", " ".join(cmd))
    try:
        subprocess.run(  # noqa: S603
            cmd,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        _LOGGER.critical("An error occurred while removing the release branch.")
        _LOGGER.critical(exc)
        sys.exit(1)


def __wait_done_request(msg: str) -> None:
    done_answer = "n"
    while done_answer != "y":
        print()  # noqa: T201
        done_answer = input(msg)
    print()  # noqa: T201


if __name__ == "__main__":
    APP()
