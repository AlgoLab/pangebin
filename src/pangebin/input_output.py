"""Input-Output PlasBin-Flow module."""

import os


def is_gz_file(filepath: str | os.PathLike) -> bool:
    """Check if a file is gzipped.

    Parameters
    ----------
    filepath : str | os.PathLike
        Path of file to check.

    Returns
    -------
    bool
        True if file is gzipped.

    """
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"
