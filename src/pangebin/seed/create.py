"""Seed contig creation module."""

import logging
from itertools import chain
from pathlib import Path
from statistics import mean

import pangebin.gene_density.input_output as gd_io
import pangebin.ground_truth.input_output as gt_io
import pangebin.seed.config as seed_config
import pangebin.seed.input_output as seed_io
import pangebin.seed.items as seed_items

_LOGGER = logging.getLogger(__name__)


def thresholds_from_datatest(
    datatest: Path,
    threshold_ranges: seed_config.ThresholdRanges,
) -> seed_items.SeedContigThresholds:
    """Compute threshold."""
    plasmid_ctg_len_gd_plasmid, nonplasmid_ctg_len_gd = (
        dispath_plasmid_and_nonplasmid_attributes(datatest)
    )
    mean_length: float = mean(  # type: ignore[assignment]
        (
            value[0]  # type: ignore[type-var]
            for value in chain(
                plasmid_ctg_len_gd_plasmid.values(),
                nonplasmid_ctg_len_gd.values(),
            )
        ),
    )
    mean_gene_density: float = mean(  # type: ignore[assignment]
        (
            value[1]  # type: ignore[type-var]
            for value in chain(
                plasmid_ctg_len_gd_plasmid.values(),
                nonplasmid_ctg_len_gd.values(),
            )
        ),
    )
    best_thresholds, best_sp_nps = get_best_thresholds(
        plasmid_ctg_len_gd_plasmid,
        nonplasmid_ctg_len_gd,
        threshold_ranges,
    )

    return seed_items.SeedContigThresholds(
        mean_length,
        mean_gene_density,
        best_sp_nps,
        best_thresholds,
    )


def dispath_plasmid_and_nonplasmid_attributes(
    datatest: Path,
) -> tuple[
    dict[str, tuple[int, float, str]],
    dict[str, tuple[int, float]],
]:
    """Dispatch plasmid and nonplasmid attributes."""
    plasmid_ctg_len_gd_plasmid: dict[str, tuple[int, float, str]] = {}
    nonplasmid_ctg_len_gd: dict[str, tuple[int, float]] = {}

    for k, test_item in enumerate(seed_io.test_items_from_file(datatest)):
        for pls_ctg in gt_io.plasmid_contigs_from_file(
            test_item.plasmid_contigs_file(),
        ):
            ctg_id = f"{k}_{pls_ctg.contig_id()}"
            plasmid_ctg_len_gd_plasmid[ctg_id] = (
                pls_ctg.contig_length(),
                0.0,
                pls_ctg.plasmid_id(),
            )

        for nonpls_ctg in gt_io.non_plasmid_contigs_from_file(
            test_item.non_plasmid_contigs_file(),
        ):
            ctg_id = f"{k}_{nonpls_ctg.contig_id()}"
            nonplasmid_ctg_len_gd[ctg_id] = (nonpls_ctg.contig_length(), 0.0)

        for original_ctg_id, gene_density in gd_io.from_file(
            test_item.contig_gene_densities_file(),
        ):
            ctg_id = f"{k}_{original_ctg_id}"
            if ctg_id in plasmid_ctg_len_gd_plasmid:
                plasmid_ctg_len_gd_plasmid[ctg_id] = (
                    plasmid_ctg_len_gd_plasmid[ctg_id][0],
                    gene_density,
                    plasmid_ctg_len_gd_plasmid[ctg_id][2],
                )
            else:
                nonplasmid_ctg_len_gd[ctg_id] = (
                    nonplasmid_ctg_len_gd[ctg_id][0],
                    gene_density,
                )

    return plasmid_ctg_len_gd_plasmid, nonplasmid_ctg_len_gd


def get_best_thresholds(
    plasmid_ctg_len_gd: dict[str, tuple[int, float, str]],
    nonplasmid_ctg_len_gd: dict[str, tuple[int, float]],
    threshold_ranges: seed_config.ThresholdRanges,
) -> tuple[list[tuple[int, float]], int]:
    """Get best thresholds that maximizes SP - NPS."""
    best_thresholds: list[tuple[int, float]] = []
    best_sp_nps = -len(nonplasmid_ctg_len_gd)
    min_length = threshold_ranges.min_length()
    while min_length <= threshold_ranges.max_length():
        min_gene_density = threshold_ranges.min_gene_density()
        while min_gene_density <= threshold_ranges.max_gene_density():
            sp_nps = sp_minus_nps(
                min_length,
                min_gene_density,
                plasmid_ctg_len_gd,
                nonplasmid_ctg_len_gd,
            )
            if sp_nps > best_sp_nps:
                best_sp_nps = sp_nps
            elif sp_nps == best_sp_nps:
                best_thresholds.append((min_length, min_gene_density))
            min_gene_density += threshold_ranges.step_gene_density()
        min_length += threshold_ranges.step_length()

    _LOGGER.info(
        "Found %i best pairs of thresholds with SP - NPS = %i",
        len(best_thresholds),
        best_sp_nps,
    )
    return best_thresholds, best_sp_nps


def sp_minus_nps(
    min_length: int,
    min_gene_density: float,
    plasmid_ctg_len_gd: dict[str, tuple[int, float, str]],
    nonplasmid_ctg_len_gd: dict[str, tuple[int, float]],
) -> int:
    """Compute SP - NPS."""
    seed_plasmids = set()
    for ctg_len, ctg_gd, pls_id in plasmid_ctg_len_gd.values():
        if ctg_len >= min_length and ctg_gd >= min_gene_density:
            seed_plasmids.add(pls_id)

    nps = sum(
        1
        for ctg_len, ctg_gd in nonplasmid_ctg_len_gd.values()
        if ctg_len >= min_length and ctg_gd >= min_gene_density
    )

    return len(seed_plasmids) - nps
