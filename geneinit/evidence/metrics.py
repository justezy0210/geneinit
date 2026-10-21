"""Small deterministic interval metrics shared by evidence integrations."""

from __future__ import annotations

from collections import defaultdict
from typing import Callable, Generic, Iterable, Iterator, Sequence, TypeVar


Interval = tuple[int, int]
T = TypeVar("T")


def merge_intervals(intervals: Iterable[Interval]) -> tuple[Interval, ...]:
    ordered = sorted((int(start), int(end)) for start, end in intervals)
    if any(start < 1 or end < start for start, end in ordered):
        raise ValueError("invalid 1-based inclusive interval")
    if not ordered:
        return ()
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def interval_bp(intervals: Iterable[Interval]) -> int:
    return sum(end - start + 1 for start, end in merge_intervals(intervals))


def overlap_bp(left: Iterable[Interval], right: Iterable[Interval]) -> int:
    left_items = merge_intervals(left)
    right_items = merge_intervals(right)
    left_index = right_index = overlap = 0
    while left_index < len(left_items) and right_index < len(right_items):
        left_start, left_end = left_items[left_index]
        right_start, right_end = right_items[right_index]
        overlap += max(0, min(left_end, right_end) - max(left_start, right_start) + 1)
        if left_end < right_end:
            left_index += 1
        else:
            right_index += 1
    return overlap


def precision_recall_f1(
    matches: int, predicted_total: int, evidence_total: int
) -> tuple[float | None, float | None, float | None]:
    precision = matches / predicted_total if predicted_total else None
    recall = matches / evidence_total if evidence_total else None
    if precision is None or recall is None:
        return precision, recall, None
    if precision + recall == 0:
        return precision, recall, 0.0
    return precision, recall, 2 * precision * recall / (precision + recall)


def interval_metrics(
    predicted: Iterable[Interval], evidence: Iterable[Interval]
) -> tuple[int, float | None, float | None, float | None]:
    predicted_items = merge_intervals(predicted)
    evidence_items = merge_intervals(evidence)
    shared = overlap_bp(predicted_items, evidence_items)
    precision, recall, f1 = precision_recall_f1(
        shared, interval_bp(predicted_items), interval_bp(evidence_items)
    )
    return shared, precision, recall, f1


def introns_from_blocks(blocks: Iterable[Interval]) -> tuple[Interval, ...]:
    merged = merge_intervals(blocks)
    return tuple(
        (left_end + 1, right_start - 1)
        for (_, left_end), (right_start, _) in zip(merged, merged[1:])
        if left_end + 1 <= right_start - 1
    )


def intron_metrics(
    predicted: Sequence[Interval], evidence: Sequence[Interval]
) -> tuple[int, float | None, float | None, float | None]:
    matches = len(set(predicted) & set(evidence))
    precision, recall, f1 = precision_recall_f1(
        matches, len(predicted), len(evidence)
    )
    return matches, precision, recall, f1


def clip_intervals(
    intervals: Iterable[Interval], start: int, end: int
) -> tuple[Interval, ...]:
    return merge_intervals(
        (max(start, left), min(end, right))
        for left, right in intervals
        if max(start, left) <= min(end, right)
    )


class IntervalIndex(Generic[T]):
    """Deterministic coarse-bin interval index with identity de-duplication."""

    def __init__(
        self,
        items: Iterable[T],
        *,
        seqid: Callable[[T], str],
        start: Callable[[T], int],
        end: Callable[[T], int],
        strand: Callable[[T], str] | None = None,
        bin_size: int = 100_000,
    ) -> None:
        if bin_size < 1:
            raise ValueError("interval-index bin size must be positive")
        self.bin_size = bin_size
        self._seqid = seqid
        self._start = start
        self._end = end
        self._strand = strand
        self._bins: dict[tuple[str, str | None, int], list[T]] = defaultdict(list)
        for item in items:
            left, right = start(item), end(item)
            if left < 1 or right < left:
                raise ValueError("invalid indexed interval")
            strand_value = strand(item) if strand else None
            for bin_id in range((left - 1) // bin_size, (right - 1) // bin_size + 1):
                self._bins[(seqid(item), strand_value, bin_id)].append(item)

    def query(
        self, seqid: str, start: int, end: int, strand: str | None = None
    ) -> Iterator[T]:
        if start < 1 or end < start:
            raise ValueError("invalid query interval")
        seen: set[int] = set()
        for bin_id in range((start - 1) // self.bin_size, (end - 1) // self.bin_size + 1):
            for item in self._bins.get((seqid, strand if self._strand else None, bin_id), ()):
                identity = id(item)
                if identity in seen:
                    continue
                seen.add(identity)
                if self._start(item) <= end and self._end(item) >= start:
                    yield item
