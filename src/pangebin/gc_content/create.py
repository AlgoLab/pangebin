"""GC content create module."""

import logging
import math
from collections.abc import Iterator
from pathlib import Path

import scipy.special as sc_spe  # type: ignore[import-untyped]
from Bio.SeqRecord import SeqRecord
from scipy import integrate

import pangebin.gfa.iter as gfa_iter
from pangebin.gc_content import items

_LOGGER = logging.getLogger(__name__)


DEFAULT_PSEUDO_COUNT = 10


def gfa_file_to_gc_scores(
    gfa_file: Path,
    gc_content_intervals: items.Intervals,
    pseudo_count: int = DEFAULT_PSEUDO_COUNT,
) -> items.IntervalsAndScores:
    """Compute GC scores for a GFA graph.

    Parameters
    ----------
    gfa_file : Path
        GFA file
    gc_content_intervals : items.Intervals
        GC content intervals
    pseudo_count : int, optional
        Pseudocount, by default DEFAULT_PSEUDO_COUNT

    Returns
    -------
    items.IntervalAndScores
        GC scores for each GC interval

    """
    intervals_and_scores = items.IntervalsAndScores.from_intervals(gc_content_intervals)
    for seq_record in gfa_iter.sequence_records(gfa_file):
        intervals_and_scores.add_sequence_scores(
            items.SequenceProbasAndScores(
                seq_record.name,
                (
                    sequence_gc_proba_and_score(
                        seq_record,
                        gc_content_intervals,
                        list(gc_content_intervals.interval_equiprobabilities()),
                        pseudo_count=pseudo_count,
                    )
                ),
            ),
        )
    return intervals_and_scores


def sequence_gc_proba_and_score(
    seq_record: SeqRecord,
    gc_content_intervals: items.Intervals,
    all_prob_b: list[float],
    pseudo_count: int,
) -> Iterator[tuple[float, float]]:
    """Compute sequence GC probabilities and scores for each GC interval.

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence record
    gc_content_intervals : items.Intervals
        GC content intervals
    all_prob_b : list of float
        Uniform probabilities of each GC interval
    pseudo_count : int
        Pseudocount

    Yield
    -----
    tuple[float, float]
        GC probability and score for the sequence record

    """

    def proba_to_count_n_in_a_plasmid_with_gc_ratio(
        seq_length: int,
        seq_gc_count: int,
        plasmid_gc_ratio: float,
        pseudo_count: int,
    ) -> float:
        """Give the probability of the contig to be in a plasmid with a given GC ratio.

        Compute probability of observing the number of GC nucleotides in the sequence
        within a molecule of GC content <plasmid_gc_ratio>
        using pseducount <pseudo_count>.

        Parameters
        ----------
        seq_length : int
            Sequence length
        seq_gc_count : int
            Number of GC nucleotides
        plasmid_gc_ratio : float
            GC ratio of the molecule
        pseudo_count : int
            Pseudocount

        Note
        ----
        Done via logarithm to avoid overflow.

        """

        def ln_n_choose_k(n: int, k: int) -> float:
            """Compute ln of n choose k.

            Note than n! = gamma(n+1)
            """
            return (
                sc_spe.gammaln(n + 1)
                - sc_spe.gammaln(k + 1)
                - sc_spe.gammaln(n - k + 1)
            )

        alpha = pseudo_count * plasmid_gc_ratio
        beta = pseudo_count * (1 - plasmid_gc_ratio)

        ln_proba_to_be_plasmid_with_gc_ratio = (
            ln_n_choose_k(seq_length, seq_gc_count)
            + sc_spe.betaln(seq_gc_count + alpha, seq_length - seq_gc_count + beta)
            - sc_spe.betaln(alpha, beta)
        )
        return math.exp(ln_proba_to_be_plasmid_with_gc_ratio)

    all_prob_n_knw_b_x_prob_b = []  # P(n|b, l) x P(b)

    seq_length = len(seq_record)
    seq_gc_count = sum(1 for nt in seq_record.seq if nt in ("G", "C"))
    _LOGGER.debug(
        "GC ratio of %s is %s (len = %s)",
        seq_record.name,
        seq_gc_count / seq_length,
        seq_length,
    )
    for k, gc_interval in enumerate(gc_content_intervals):
        prob_n_knw_b = integrate.quad(
            lambda plasmid_gc_ratio: proba_to_count_n_in_a_plasmid_with_gc_ratio(
                seq_length,
                seq_gc_count,
                plasmid_gc_ratio,
                pseudo_count,
            ),
            *gc_interval,
        )
        prob_n_knw_b_x_prob_b = prob_n_knw_b[0] / all_prob_b[k]
        all_prob_n_knw_b_x_prob_b.append(prob_n_knw_b_x_prob_b)

    max_prob_n_knw_b_x_prob_b = max(all_prob_n_knw_b_x_prob_b)

    normalized_probs = [
        prob_n_knw_b_x_prob_b / max_prob_n_knw_b_x_prob_b
        for prob_n_knw_b_x_prob_b in all_prob_n_knw_b_x_prob_b
    ]
    for normalized_prob in normalized_probs:
        yield (normalized_prob, 2 * normalized_prob - 1)
