"""Assemble reads."""

from pathlib import Path

import pangebin.assembly.config as asm_config
from pangebin import subprocess_lib

UNICYCLER_CMD_STR = "unicycler"


def unicycler_paired_end_reads(
    reads_1_fastq: Path,
    reads_2_fastq: Path,
    out_directory: Path,
    config: asm_config.Unicycler,
) -> None:
    """Assemble paired short reads with Unicycler."""
    cmd_path = subprocess_lib.command_path(UNICYCLER_CMD_STR)
    cli_line = [
        cmd_path,
        "-1",
        reads_1_fastq,
        "-2",
        reads_2_fastq,
        "-o",
        out_directory,
        "--threads",
        config.number_threads(),
    ]
    subprocess_lib.run_cmd(cli_line, UNICYCLER_CMD_STR)
