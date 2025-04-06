"""GFA views module."""

from pathlib import Path

import rich

import pangebin.gfa.line as gfa_line
import pangebin.input_output as io


def print_stats(gfa_path: Path) -> None:
    """Print stats of a GFA graph."""
    number_of_segments = 0
    number_of_links = 0
    number_of_paths = 0
    with io.open_file_read(gfa_path) as f_in:
        for line in f_in:
            if line[0] == gfa_line.Type.SEGMENT:
                number_of_segments += 1
            elif line[0] == gfa_line.Type.LINK:
                number_of_links += 1
            elif line[0] == gfa_line.Type.PATH:
                number_of_paths += 1
    rich.print(f"Number of segments: {number_of_segments}")
    rich.print(f"Number of links: {number_of_links}")
    rich.print(f"Number of paths: {number_of_paths}")
