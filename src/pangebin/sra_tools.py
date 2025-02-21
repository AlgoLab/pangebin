"""SRA tools module."""

import shutil
from pathlib import Path

from pangebin import subprocess_lib

PREFETCH_CMD = "prefetch"


def run_prefetch(srr_id: str) -> Path:
    """Run prefetch.

    Parameters
    ----------
    srr_id : str
        SRR ID.

    Returns
    -------
    Path
        Path of the downloaded directory.
    """
    cmd_path = subprocess_lib.command_path(PREFETCH_CMD)
    cli_line = [cmd_path, srr_id]
    subprocess_lib.run_cmd(cli_line, PREFETCH_CMD)
    return Path(f"{srr_id}")


FASTQ_DUMP_CMD = "fastq-dump"


def run_fastq_dump(input_dir: Path) -> Path:
    """Run fastq-dump.

    Parameters
    ----------
    input_dir : Path
        Input directory.

    Returns
    -------
    Path
        Path of the FASTQ file.
    """
    cmd_path = subprocess_lib.command_path(FASTQ_DUMP_CMD)
    cli_line = [cmd_path, input_dir]
    subprocess_lib.run_cmd(cli_line, FASTQ_DUMP_CMD)
    return Path(f"{input_dir}.fastq")


def run_fastq_dump_paired(input_dir: Path) -> tuple[Path, Path]:
    """Run fastq-dump for paired-end reads.

    Parameters
    ----------
    input_dir : Path
        Input directory.

    Returns
    -------
    Path
        Path of the first FASTQ file.
    Path
        Path of the second FASTQ file.
    """
    cmd_path = subprocess_lib.command_path(FASTQ_DUMP_CMD)
    cli_line = [cmd_path, "--split-3", input_dir]
    subprocess_lib.run_cmd(cli_line, FASTQ_DUMP_CMD)
    return Path(f"{input_dir}_1.fastq"), Path(f"{input_dir}_2.fastq")


def download_sra_fastq(
    srr_id: str,
    output_dir: Path,
    remove_dir: bool = True,  # noqa: FBT001, FBT002
) -> Path:
    """Download a SRA file.

    Parameters
    ----------
    srr_id : str
        SRR ID.
    output_dir : Path
        Output directory where to save the FASTQ.
    remove_dir : bool, optional
        Remove the downloaded directory, by default True

    Returns
    -------
    Path
        Path of the FASTQ file.
    """
    input_dir = run_prefetch(srr_id)
    fastq_file = run_fastq_dump(input_dir)
    if remove_dir:
        shutil.rmtree(input_dir, ignore_errors=True)
    return Path(shutil.move(fastq_file, output_dir))


def download_sra_fastq_paired(
    srr_id: str,
    output_dir: Path,
    remove_dir: bool = True,  # noqa: FBT001, FBT002
) -> tuple[Path, Path]:
    """Download a SRA file.

    Parameters
    ----------
    srr_id : str
        SRR ID.
    output_dir : Path
        Output directory where to save the FASTQ.
    remove_dir : bool, optional
        Remove the downloaded prefetch directory, by default True

    Returns
    -------
    Path
        Path of the first FASTQ file.
    Path
        Path of the second FASTQ file.
    """
    input_dir = run_prefetch(srr_id)
    fastq_file_1, fastq_file_2 = run_fastq_dump_paired(input_dir)
    if remove_dir:
        shutil.rmtree(input_dir, ignore_errors=True)
    return Path(shutil.move(fastq_file_1, output_dir)), Path(
        shutil.move(fastq_file_2, output_dir),
    )
