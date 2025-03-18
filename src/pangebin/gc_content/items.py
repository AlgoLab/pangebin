"""GC content items."""

from __future__ import annotations

import yaml  # type: ignore[import-untyped]

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

import logging

_LOGGER = logging.getLogger(__name__)


class SequenceProbasAndScores(list[tuple[float, float]]):
    """GC content score for a sequence."""

    def __init__(
        self,
        seq_id: str,
        gc_prob_scores: Iterable[tuple[float, float]],
    ) -> None:
        """Initialize object."""
        super().__init__(gc_prob_scores)
        self.__seq_id = seq_id

    def sequence_id(self) -> str:
        """Sequence ID."""
        return self.__seq_id

    def probas(self) -> list[float]:
        """Probas."""
        return [proba for proba, _ in self]

    def scores(self) -> list[float]:
        """Scores."""
        return [score for _, score in self]

    def to_dict(self) -> dict[str, list[list[float]]]:
        """Convert to dictionary."""
        return {self.__seq_id: [list(interval) for interval in self]}


class Intervals:
    """GC content intervals."""

    __MIN_FLOAT = 0.0
    __MAX_FLOAT = 1.0

    __step_list: list[float]

    @classmethod
    def from_file(cls, file: Path) -> Intervals:
        """Read intervals from file."""
        intervals = []
        with file.open() as f_in:
            for line in f_in:
                intervals.append(float(line))  # noqa: PERF401
        return cls(intervals)

    def __init__(self, step_list: Iterable[float] | None = None) -> None:
        """Initialize object."""
        self.__step_list = []
        if step_list is not None:
            for step in step_list:
                self.append(step)

    def append(self, step: float) -> None:
        """Append interval."""
        if len(self.__step_list) == 0:
            if step != self.__MIN_FLOAT:
                _err_msg = (
                    f"Expected {self.__MIN_FLOAT} for the first interval, got {step}"
                )
                _LOGGER.error(_err_msg)
                raise ValueError(_err_msg)
        elif step < self.__step_list[-1]:
            _err_msg = f"Expected increasing steps, got {step} < {self.__step_list[-1]}"
            _LOGGER.error(_err_msg)
            raise ValueError(_err_msg)
        elif step > self.__MAX_FLOAT:
            _err_msg = f"Expected {self.__MAX_FLOAT} for the last interval, got {step}"
            _LOGGER.error(_err_msg)
            raise ValueError(_err_msg)
        self.__step_list.append(step)

    def __iter__(self) -> Iterator[tuple[float, float]]:
        """Get interval iterator."""
        yield from zip(self.__step_list[:-1], self.__step_list[1:])

    def interval_equiprobabilities(self) -> Iterator[float]:
        """Get interval equiprobability iterator."""
        for interval in self:
            yield (interval[1] - interval[0])

    def __len__(self) -> int:
        """Get the number of intervals."""
        return len(self.__step_list) - 1

    def to_list(self) -> list[float]:
        """Convert to list."""
        return self.__step_list

    def to_file(self, file: Path) -> None:
        """Write intervals to file."""
        with file.open("w") as f_out:
            for increasing_steps in self.__step_list:
                f_out.write(f"{increasing_steps}\n")


class IntervalsAndScores:
    """GC content intervals and scores."""

    KEY_INTERVALS = "intervals"
    KEY_SEQUENCES = "sequences"

    @classmethod
    def from_yaml(cls, file: Path) -> IntervalsAndScores:
        """Create object from yaml file."""
        with file.open() as f_in:
            dict_from_yaml = yaml.safe_load(f_in)
        return cls.from_dict(dict_from_yaml)

    @classmethod
    def from_dict(cls, dict_from_yaml: dict) -> IntervalsAndScores:
        """Create object from dict."""
        return cls(
            Intervals(dict_from_yaml[cls.KEY_INTERVALS]),
            (
                SequenceProbasAndScores(seq_id, seq_scores)
                for seq_id, seq_scores in dict_from_yaml[cls.KEY_SEQUENCES].items()
            ),
        )

    @classmethod
    def from_intervals(cls, intervals: Intervals) -> IntervalsAndScores:
        """Create object from intervals."""
        return cls(intervals, [])

    def __init__(
        self,
        intervals: Intervals,
        scores: Iterable[SequenceProbasAndScores],
    ) -> None:
        """Initialize object."""
        self.__intervals = intervals
        self.__scores = list(scores)

    def intervals(self) -> Intervals:
        """Intervals."""
        return self.__intervals

    def scores(self) -> Iterator[SequenceProbasAndScores]:
        """Get sequence probas and scores iterator."""
        yield from self.__scores

    def add_sequence_scores(self, score: SequenceProbasAndScores) -> None:
        """Add sequence scores."""
        if len(score) != len(self.intervals()):
            _err_msg = f"Expected {len(self.intervals())} scores, got {len(score)}"
            _LOGGER.error(_err_msg)
            raise ValueError(_err_msg)
        self.__scores.append(score)

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            self.KEY_INTERVALS: self.intervals().to_list(),
            self.KEY_SEQUENCES: {
                seq_id: seq_scores
                for score in self.scores()
                for seq_id, seq_scores in score.to_dict().items()
            },
        }

    def to_yaml(self, file: Path) -> Path:
        """Convert to yaml file."""
        with file.open("w") as f_out:
            yaml.dump(self.to_dict(), f_out, Dumper)
        return file
