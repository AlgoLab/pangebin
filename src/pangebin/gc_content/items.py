"""GC content items."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

import logging

_LOGGER = logging.getLogger(__name__)


class Intervals:
    """GC content intervals."""

    __MIN_FLOAT = 0.0
    __MAX_FLOAT = 1.0

    __step_list: list[float]

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


class IntervalFormatter:
    """Interval formatter."""

    SEP = "_"

    @classmethod
    def to_str(cls, interval: tuple[float, float]) -> str:
        """Format interval."""
        return f"{interval[0]}{cls.SEP}{interval[1]}"

    @classmethod
    def from_str(cls, interval_str: str) -> tuple[float, float]:
        """Parse interval."""
        l_split = interval_str.split(cls.SEP)
        return (float(l_split[0]), float(l_split[1]))


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
