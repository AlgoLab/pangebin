"""Standard assembly FASTA file."""

from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO

from pangebin import assembler

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

_PANSN_SAMPLE = "StandardMix"
_PANSN_SEP = "#"


def fastas_to_pansn_mixed_fasta(
    skesa_fasta_path: Path,
    unicycler_fasta_path: Path,
    mixed_fasta_path: Path,
) -> None:
    """Merge FASTA files according to the PanSN standard.

    See https://github.com/pangenome/PanSN-spec
    """
    with mixed_fasta_path.open("w") as mixed_fasta_file:
        for fasta_in, haplotype_id in (
            (skesa_fasta_path, assembler.HaplotypeID.SKESA),
            (unicycler_fasta_path, assembler.HaplotypeID.UNICYCLER),
        ):
            seq_record: SeqRecord
            for seq_record in SeqIO.parse(fasta_in, "fasta"):
                fmt_name = fmt_seq_name_to_pansn_spec(seq_record.name, haplotype_id)
                mixed_fasta_file.write(
                    f">{fmt_name}\n",
                )
                mixed_fasta_file.write(f"{seq_record.seq!s}\n")


def fmt_seq_name_to_pansn_spec(
    sequence_name: str,
    haplotype_id: assembler.HaplotypeID,
) -> str:
    """Format sequence name according to the PanSN standard.

    See https://github.com/pangenome/PanSN-spec
    """
    return _PANSN_SEP.join(
        [
            _PANSN_SAMPLE,
            str(haplotype_id),
            sequence_name,
        ],
    )


def pansn_name_to_contig_name(pansn_name: str) -> str:
    """Get contig name from PanSN name."""
    return pansn_name.split(_PANSN_SEP)[2]
