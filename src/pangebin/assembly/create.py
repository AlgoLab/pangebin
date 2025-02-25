"""Assemble reads."""

from pathlib import Path

from pangebin import subprocess_lib

UNICYCLER_CMD_STR = "unicycler"


def unicycler_paired_end_reads(
    reads_1_fastq: Path,
    reads_2_fastq: Path,
    out_directory: Path,
    ressources_config: subprocess_lib.RessourcesConfig,
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
        ressources_config.max_cores(),
        "--spades_options",
        f"-m {ressources_config.max_memory()}",
    ]
    subprocess_lib.run_cmd(cli_line, UNICYCLER_CMD_STR)
