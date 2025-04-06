"""Run Blast."""

from pathlib import Path

from pangebin import subprocess_lib

MAKEBLASTDB_CMD_STR = "makeblastdb"
BLAST_CMD_STR = "blastn"


def blast_map(
    query_fasta_file: Path,
    subject_fasta_file: Path,
    out_mapping_file: Path,
) -> None:
    """Blast the query against the subject."""
    subject_db_file = subject_fasta_file.with_suffix(".db")

    cmd_path = subprocess_lib.command_path(MAKEBLASTDB_CMD_STR)
    cli_line = [
        cmd_path,
        "-in",
        subject_fasta_file,
        "-dbtype",
        "nucl",
        "-out",
        subject_db_file,
    ]
    subprocess_lib.run_cmd(cli_line, MAKEBLASTDB_CMD_STR)

    cmd_path = subprocess_lib.command_path(BLAST_CMD_STR)
    cli_line = [
        "blastn",
        "-task",
        "megablast",
        "-query",
        query_fasta_file,
        "-db",
        subject_db_file,
        "-out",
        out_mapping_file,
        "-outfmt",
        "6",
    ]
    subprocess_lib.run_cmd(cli_line, BLAST_CMD_STR)

    for file in subject_db_file.parent.glob("*.db.n*"):
        file.unlink()
