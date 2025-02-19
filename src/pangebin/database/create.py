"""Create database module."""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING
from urllib.error import HTTPError, URLError

from Bio import Entrez
from defusedxml import ElementTree as DET  # type: ignore[import-untyped] # noqa: N814
from rich import progress

import pangebin.database.items as db_items
import pangebin.entrez as pg_entrez
import pangebin.logging as common_log

if TYPE_CHECKING:
    from collections.abc import Iterable

_LOGGER = logging.getLogger(__name__)


def illumina_exp_sra_ids_from_plasmid_accessions(
    accessions: Iterable[str],
    entrez_config: pg_entrez.Config,
) -> tuple[
    list[db_items.IlluminaBioSamples],
    list[db_items.NonIlluminaBioSamples],
]:
    """Retrieve paired Illumina experiment SRA ids from plasmid accessions."""
    pg_entrez.set_entrez_config(entrez_config)

    illumina_biosamples: dict[str, db_items.IlluminaBioSamples] = {}
    non_illumina_biosamples: dict[str, db_items.NonIlluminaBioSamples] = {}

    _nb_plasmids_without_biosample = 0
    _nb_biosamples_without_illumina_paired_reads = 0

    accessions = list(accessions)
    try:
        handle = Entrez.efetch(
            db="nuccore",
            id=accessions,
            rettype="gb",
            retmode="xml",
            seq_start=1,
            seq_stop=1,
        )
    except URLError as e:
        _err_msg = f"Failed to fetch GenBank record for {accessions}"
        _LOGGER.exception(_err_msg)
        raise ValueError(_err_msg) from e

    records = Entrez.read(handle)
    handle.close()

    for record in progress.track(
        records,
        description="Processing plasmid accessions...",
        console=common_log.CONSOLE,
    ):
        genblank_primary_accession = str(record["GBSeq_primary-accession"])

        biosample_id = biosample_ids_from_genbank_record(record)

        if biosample_id is None:
            _nb_plasmids_without_biosample += 1
            _LOGGER.warning(
                "No biosample id for genbank %s",
                genblank_primary_accession,
            )
        elif biosample_id in illumina_biosamples:
            illumina_biosamples[biosample_id].plasmids().append(
                genblank_primary_accession,
            )
        elif biosample_id in non_illumina_biosamples:
            non_illumina_biosamples[biosample_id].plasmids().append(
                genblank_primary_accession,
            )
        else:
            sample_sra_id = sample_sra_id_from_biosample(biosample_id, entrez_config)

            illumina_exp_sra_ids = illumina_exp_sra_ids_from_sample_sra_id(
                sample_sra_id,
            )

            if not illumina_exp_sra_ids:
                _nb_biosamples_without_illumina_paired_reads += 1
                _LOGGER.warning(
                    "No Illumina paired reads for biosample %s",
                    biosample_id,
                )
                non_illumina_biosamples[biosample_id] = db_items.NonIlluminaBioSamples(
                    biosample_id,
                    sample_sra_id,
                    [genblank_primary_accession],
                )
            else:
                illumina_biosamples[biosample_id] = db_items.IlluminaBioSamples(
                    biosample_id,
                    sample_sra_id,
                    [genblank_primary_accession],
                    illumina_exp_sra_ids,
                )

    if _nb_plasmids_without_biosample > 0:
        _LOGGER.warning(
            "Plasmids without biosample: %s",
            _nb_plasmids_without_biosample,
        )
    if _nb_biosamples_without_illumina_paired_reads > 0:
        _LOGGER.warning(
            "Biosamples without Illumina paired reads: %s",
            _nb_biosamples_without_illumina_paired_reads,
        )

    return (
        list(illumina_biosamples.values()),
        list(non_illumina_biosamples.values()),
    )


def biosample_ids_from_genbank_record(record: dict) -> str | None:
    """Extract biosample ids from genbank record."""
    if "GBSeq_xrefs" in record:
        for xref in record["GBSeq_xrefs"]:
            if xref["GBXref_dbname"] == "BioSample":
                return str(xref["GBXref_id"])
    return None


def sample_sra_id_from_biosample(
    biosample_id: str,
    entrez_config: pg_entrez.Config,
) -> str:
    """Extract sample SRA id from biosample id."""
    for i in range(entrez_config.max_tries()):
        try:
            handle = Entrez.efetch(
                db="biosample",
                id=biosample_id,
                rettype="full",
                retmode="xml",
            )
        except HTTPError as e:  # noqa: PERF203
            # XXX efetch biosample db can send unexpected 400 error, while id exists
            if e.code != 400 or i == entrez_config.max_tries() - 1:  # noqa: PLR2004
                raise
            _LOGGER.debug(
                "Fetching Biosample %s returns a 400 HTTP error, retrying",
                biosample_id,
            )
            time.sleep(entrez_config.sleep_between_tries())
        except URLError as e:
            _err_msg = f"Failed to fetch Biosample {biosample_id}"
            _LOGGER.exception(_err_msg)
            raise ValueError(_err_msg) from e
        else:
            break

    xml_data = handle.read()
    handle.close()

    root = DET.fromstring(xml_data)

    # Extract the SRA identifiers (which are inside the <Ids> section with db="SRA")
    return str(root.find(".//Id[@db='SRA']").text)


def illumina_exp_sra_ids_from_sample_sra_id(sample_sra_id: str) -> list[str]:
    """Extract illumina experiment SRA ids linked to the sample SRA id."""
    srx_sra_ids: list[str] = []
    handle = Entrez.esearch(db="SRA", term=sample_sra_id)
    record = Entrez.read(handle, validate=False)
    handle.close()
    sra_summaries = summaries("sra", record["IdList"])
    for sra_summary in sra_summaries:
        srx_sra = illumina_exp_sra_from_summary(sra_summary)
        if srx_sra is not None:
            srx_sra_ids.append(srx_sra)
    return srx_sra_ids


def summaries(
    database: str,
    identifier: str | list[str],
) -> dict:
    """Retrieve Entrez summaries."""
    handle = Entrez.esummary(db=database, id=identifier)
    record = Entrez.read(handle)
    handle.close()
    return record


def illumina_exp_sra_from_summary(sra_summary: dict) -> str | None:
    """Extract paired illumina experiment SRA from summary."""
    sra_info = f"<ExpXml>{sra_summary['ExpXml']}</ExpXml>"

    root = DET.fromstring(sra_info)

    platform = root.find("Summary").find("Platform").text
    if "illumina" not in platform.lower():
        return None

    child = root
    for tag in ("Library_descriptor", "LIBRARY_LAYOUT", "PAIRED"):
        child = child.find(tag)
        if child is None:
            return None

    return str(root.find("Experiment").get("acc"))
