"""PlasBin-flow conversion items module."""

from collections.abc import Iterable


# REFACTOR repetition of SequenceNormCoverage
class ContigMult:
    """PlasBin-flow contig multiplicity."""

    def __init__(self, contig_id: str, contig_mult: float) -> None:
        self.__contig_id = contig_id
        self.__contig_mult = contig_mult

    def identifier(self) -> str:
        """Contig ID."""
        return self.__contig_id

    def multiplicity(self) -> float:
        """Contig multiplicity."""
        return self.__contig_mult


class PBFBinInfo:
    """PlasBin-flow bin info."""

    def __init__(
        self,
        bin_id: str,
        flow: float,
        gc_interval: tuple[float, float],
        contigs_mults: Iterable[ContigMult],
    ) -> None:
        self.__bin_id = bin_id
        self.__flow = flow
        self.__gc_interval = gc_interval
        self.__contigs_mults = list(contigs_mults)

    def identifier(self) -> str:
        """Bin ID."""
        return self.__bin_id

    def flow(self) -> float:
        """Flow."""
        return self.__flow

    def gc_interval(self) -> tuple[float, float]:
        """GC interval."""
        return self.__gc_interval

    def contigs_mults(self) -> Iterable[ContigMult]:
        """Contigs multiplicity."""
        return self.__contigs_mults
