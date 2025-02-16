"""Ground truth create module."""

from __future__ import annotations

from typing import TYPE_CHECKING

from Bio import Entrez, SeqIO

import pangebin.ground_truth.config as gt_config
import pangebin.ground_truth.items as gt_items
import pangebin.mapping.create as map_create
import pangebin.mapping.filter as map_filter
import pangebin.mapping.input_output as map_io
import pangebin.mapping.ops as map_ops

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

import logging

_LOGGER = logging.getLogger(__name__)


def separate_plasmid_and_other_contigs(
    contigs_fasta: Path,
    plasmid_genbank_ids: Iterable[str],
    output_dir: Path,
    config: gt_config.Config | None = None,
    email_address: str | None = None,
) -> tuple[Path, Path]:
    """Create files of plasmid and non-plasmid contigs.

    Parameters
    ----------
    contigs_fasta : Path
        Path to contigs FASTA file.
    plasmid_genbank_ids : Iterable[str]
        GenBank IDs of plasmid sequences.
    output_dir : Path
        Path to output directory.
    config : gd_config.Config, optional
        Configuration object.
    email_address : str | None, optional
        Email address to fetch NCBI database, by default None

    Returns
    -------
    Path
        Plasmid contigs file
    Path
        Non-plasmid contigs file
    """
    if config is None:
        config = gt_config.Config()

    if email_address is not None:
        Entrez.email = email_address  # type: ignore[assignment]

    plasmid_contigs_file = output_dir / "plasmid_contigs.tsv"
    non_plasmid_contigs_file = output_dir / "non_plasmid_contigs.tsv"

    merged_fasta, plasmid_lengths = __merged_plasmid_sequences_fasta(
        plasmid_genbank_ids,
        output_dir,
    )
    mapping_file = output_dir / "contigs_vs_plasmids.sam"
    map_create.blast_map(contigs_fasta, merged_fasta, mapping_file)
    filtered_mapping_file = output_dir / "filtered_contigs_vs_plasmids.sam"
    map_io.to_filtered_sam_file(
        mapping_file,
        filtered_mapping_file,
        map_filter.Config(min_pident=config.min_pident()),
    )
    q_intervals_for_each_s = map_ops.queries_intervals_for_each_subject_from_sam(
        filtered_mapping_file,
    )

    contig_lengths: dict[str, int] = {
        seq_record.id: len(seq_record.seq)
        for seq_record in SeqIO.parse(contigs_fasta, "fasta")
    }

    plasmid_contigs: list[gt_items.PlasmidContig] = []
    non_plasmid_contigs: list[gt_items.NonPlasmidContig] = []

    plasmid_contig_id: set[str] = set()

    for plasmid_id, all_contig_intervals in q_intervals_for_each_s.items():
        for contig_id, contig_intervals in all_contig_intervals.items():
            disjoint_contig_intervals = map_ops.interval_union(contig_intervals)
            contig_coverage = (
                sum(len(interval) for interval in disjoint_contig_intervals)
                / contig_lengths[contig_id]
            )
            if contig_coverage >= config.min_contig_coverage():
                plasmid_contig_id.add(contig_id)
                plasmid_contigs.append(
                    gt_items.PlasmidContig(
                        plasmid_id=plasmid_id,
                        contig_id=contig_id,
                        plasmid_length=plasmid_lengths[plasmid_id],
                        contig_length=contig_lengths[contig_id],
                        contig_coverage=contig_coverage,
                    ),
                )

    for contig_id, contig_length in contig_lengths.items():
        if contig_id not in plasmid_contig_id:
            non_plasmid_contigs.append(
                gt_items.NonPlasmidContig(
                    contig_id=contig_id,
                    contig_length=contig_length,
                ),
            )

    _LOGGER.info("Write plasmid contigs file: %s", plasmid_contigs_file)
    with plasmid_contigs_file.open("w") as f_out:
        for plasmid_contig in plasmid_contigs:
            f_out.write(plasmid_contig.to_tsv_row() + "\n")

    _LOGGER.info("Write non plasmid contigs file: %s", non_plasmid_contigs_file)
    with non_plasmid_contigs_file.open("w") as f_out:
        for non_plasmid_contig in non_plasmid_contigs:
            f_out.write(non_plasmid_contig.to_tsv_row() + "\n")

    return (plasmid_contigs_file, non_plasmid_contigs_file)


def __merged_plasmid_sequences_fasta(
    genbank_ids: Iterable[str],
    output_dir: Path,
) -> tuple[Path, dict[str, int]]:
    _LOGGER.info("Merge plasmid sequences.")
    merged_fasta = output_dir / "merged_plasmid_sequences.fasta"

    plasmid_lengths: dict[str, int] = {}
    with merged_fasta.open("w") as f_out:
        for genbank_id in genbank_ids:
            stream = Entrez.efetch(
                db="nuccore",
                id=genbank_id,
                rettype="gb",
                retmode="text",
            )

            seq_record = SeqIO.read(stream, "genbank")
            f_out.write(f">{seq_record.id}\n")
            f_out.write(f"{seq_record.seq}\n")
            plasmid_lengths[seq_record.id] = len(seq_record.seq)

    return merged_fasta, plasmid_lengths
