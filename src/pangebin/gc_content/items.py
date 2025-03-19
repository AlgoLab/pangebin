"""GC content items."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

import logging

_LOGGER = logging.getLogger(__name__)


class Intervals:
    """GC content intervals."""

    __MIN_FLOAT = 0.0
    __MAX_FLOAT = 1.0

    __step_list: list[float]

    @classmethod
    def from_file(cls, file: Path) -> Intervals:
        """Read intervals from file."""
        intervals = cls()
        with file.open() as f_in:
            for line in f_in:
                intervals.append(float(line))
        return intervals

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
        yield from zip(self.__step_list[:-1], self.__step_list[1:], strict=False)

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


class SequenceGCScores:
    """Sequence probability scores."""

    def __init__(
        self,
        sequence_id: str,
        prob_score: Iterable[float],
    ) -> None:
        """Initialize object."""
        self.__sequence_id = sequence_id
        self.__prob_score = list(prob_score)

    def sequence_id(self) -> str:
        """Sequence ID."""
        return self.__sequence_id

    def probability_scores(self) -> list[float]:
        """Probability score."""
        return self.__prob_score
