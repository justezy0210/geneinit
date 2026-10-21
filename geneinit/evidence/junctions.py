"""Strand-aware exact splice-junction counting and recurrence aggregation."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import csv
import gzip
import heapq
import os
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence


MATCH_CIGAR_OPS = {0, 7, 8}
REFERENCE_CIGAR_OPS = {0, 2, 3, 7, 8}
JUNCTION_COLUMNS = (
    "seqid",
    "intron_start",
    "intron_end",
    "strand",
    "primary_reads",
    "unique_reads",
    "hq_reads",
    "max_min_anchor",
    "max_mapq",
)
AGGREGATE_COLUMNS = (
    "seqid",
    "intron_start",
    "intron_end",
    "strand",
    "total_primary_reads",
    "total_unique_reads",
    "total_hq_reads",
    "max_sample_hq_reads",
    "sample_count_hq_ge1",
    "sample_count_hq_ge3",
    "sample_count_hq_ge5",
    "condition_count_hq_ge1",
    "condition_count_hq_ge3",
    "condition_count_ge2rep_hq_ge1",
    "condition_count_ge2rep_hq_ge3",
    "condition_count_ge3rep_hq_ge1",
    "condition_count_ge3rep_hq_ge3",
    "max_min_anchor",
    "max_mapq",
    "hq_sample_counts",
)


@dataclass(slots=True)
class JunctionCount:
    primary_reads: int = 0
    unique_reads: int = 0
    hq_reads: int = 0
    max_min_anchor: int = 0
    max_mapq: int = 0


@dataclass(frozen=True, slots=True)
class JunctionRecord:
    seqid: str
    start: int
    end: int
    strand: str
    primary_reads: int
    unique_reads: int
    hq_reads: int
    max_min_anchor: int
    max_mapq: int


def cigar_junctions(
    reference_start: int, cigartuples: Iterable[tuple[int, int]]
) -> list[tuple[int, int, int, int]]:
    """Return 1-based inclusive introns and adjacent aligned anchors."""

    cigar = list(cigartuples)
    anchors: list[int] = []
    current = 0
    for operation, length in cigar:
        if operation == 3:
            anchors.append(current)
            current = 0
        elif operation in MATCH_CIGAR_OPS:
            current += length
    anchors.append(current)

    result: list[tuple[int, int, int, int]] = []
    position = reference_start
    junction_index = 0
    for operation, length in cigar:
        if operation == 3:
            result.append(
                (position + 1, position + length, anchors[junction_index], anchors[junction_index + 1])
            )
            junction_index += 1
        if operation in REFERENCE_CIGAR_OPS:
            position += length
    return result


def transcript_strand(alignment: Any, protocol: str) -> tuple[str, str]:
    """Prefer XS; otherwise infer transcript strand from the declared protocol."""

    if alignment.has_tag("XS"):
        strand = alignment.get_tag("XS")
        if strand in {"+", "-"}:
            return strand, "XS"
    if protocol == "unstranded":
        return ".", "unknown"

    aligned = "-" if alignment.is_reverse else "+"
    opposite = "+" if aligned == "-" else "-"
    if protocol == "R":
        return opposite, "inferred"
    if protocol == "F":
        return aligned, "inferred"
    if protocol == "RF":
        if alignment.is_read2:
            return aligned, "inferred"
        return opposite, "inferred"
    if protocol == "FR":
        if alignment.is_read2:
            return opposite, "inferred"
        return aligned, "inferred"
    raise ValueError(f"unsupported RNA strandedness protocol: {protocol}")


def primary_usable(alignment: Any) -> bool:
    return not (
        alignment.is_unmapped
        or alignment.is_secondary
        or alignment.is_supplementary
        or alignment.is_qcfail
        or alignment.is_duplicate
    )


def count_bam_junctions(
    bam_path: Path,
    *,
    protocol: str,
    min_mapq: int,
    min_anchor: int,
    io_threads: int,
) -> tuple[dict[tuple[str, int, int, str], JunctionCount], dict[str, int], list[str]]:
    import pysam

    counts: dict[tuple[str, int, int, str], JunctionCount] = defaultdict(JunctionCount)
    stats: dict[str, int] = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), "rb", threads=io_threads) as bam:
        references = list(bam.references)
        for alignment in bam.fetch():
            stats["alignments_seen"] += 1
            if not primary_usable(alignment):
                stats["alignments_skipped_by_flag"] += 1
                continue
            if not alignment.cigartuples or not any(
                operation == 3 for operation, _ in alignment.cigartuples
            ):
                continue
            stats["primary_spliced_alignments"] += 1
            strand, source = transcript_strand(alignment, protocol)
            seqid = bam.get_reference_name(alignment.reference_id)
            nh = alignment.get_tag("NH") if alignment.has_tag("NH") else None
            if nh is None:
                stats["spliced_alignments_missing_nh"] += 1
            unique = nh == 1
            for start, end, left_anchor, right_anchor in cigar_junctions(
                alignment.reference_start, alignment.cigartuples
            ):
                stats["primary_junction_observations"] += 1
                stats[f"strand_source_{source}_observations"] += 1
                key = (seqid, start, end, strand)
                count = counts[key]
                count.primary_reads += 1
                count.max_min_anchor = max(count.max_min_anchor, min(left_anchor, right_anchor))
                count.max_mapq = max(count.max_mapq, alignment.mapping_quality)
                if unique:
                    count.unique_reads += 1
                    stats["unique_junction_observations"] += 1
                if (
                    strand in {"+", "-"}
                    and unique
                    and alignment.mapping_quality >= min_mapq
                    and left_anchor >= min_anchor
                    and right_anchor >= min_anchor
                ):
                    count.hq_reads += 1
                    stats["hq_junction_observations"] += 1
    stats["distinct_junctions"] = len(counts)
    stats["distinct_hq_junctions"] = sum(value.hq_reads > 0 for value in counts.values())
    return counts, stats, references


def junction_rows(
    counts: dict[tuple[str, int, int, str], JunctionCount], references: Sequence[str]
) -> Iterator[tuple[Any, ...]]:
    rank = {seqid: index for index, seqid in enumerate(references)}
    for seqid, start, end, strand in sorted(
        counts,
        key=lambda key: (rank.get(key[0], len(rank)), key[1], key[2], key[3]),
    ):
        count = counts[(seqid, start, end, strand)]
        yield (
            seqid,
            start,
            end,
            strand,
            count.primary_reads,
            count.unique_reads,
            count.hq_reads,
            count.max_min_anchor,
            count.max_mapq,
        )


def recurrence_metrics(
    hq_counts: Sequence[int], sample_conditions: Sequence[str]
) -> dict[str, int]:
    if len(hq_counts) != len(sample_conditions):
        raise ValueError("HQ count and condition vectors differ in length")
    by_condition: dict[str, list[int]] = defaultdict(list)
    for count, condition in zip(hq_counts, sample_conditions, strict=True):
        by_condition[condition].append(count)

    def conditions_any(threshold: int) -> int:
        return sum(any(count >= threshold for count in counts) for counts in by_condition.values())

    def conditions_replicates(read_threshold: int, replicate_threshold: int) -> int:
        return sum(
            sum(count >= read_threshold for count in counts) >= replicate_threshold
            for counts in by_condition.values()
        )

    return {
        "sample_count_hq_ge1": sum(count >= 1 for count in hq_counts),
        "sample_count_hq_ge3": sum(count >= 3 for count in hq_counts),
        "sample_count_hq_ge5": sum(count >= 5 for count in hq_counts),
        "condition_count_hq_ge1": conditions_any(1),
        "condition_count_hq_ge3": conditions_any(3),
        "condition_count_ge2rep_hq_ge1": conditions_replicates(1, 2),
        "condition_count_ge2rep_hq_ge3": conditions_replicates(3, 2),
        "condition_count_ge3rep_hq_ge1": conditions_replicates(1, 3),
        "condition_count_ge3rep_hq_ge3": conditions_replicates(3, 3),
    }


def write_bgzf_indexed(
    output: Path, columns: Sequence[str], rows: Iterable[Sequence[Any]]
) -> None:
    import pysam

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.partial.{os.getpid()}")
    temporary_index = Path(f"{temporary}.csi")
    try:
        with pysam.BGZFile(str(temporary), "wb") as handle:
            handle.write(("#" + "\t".join(columns) + "\n").encode())
            for row in rows:
                handle.write(("\t".join(map(str, row)) + "\n").encode())
        pysam.tabix_index(
            str(temporary),
            seq_col=0,
            start_col=1,
            end_col=2,
            meta_char="#",
            zerobased=False,
            csi=True,
            force=True,
        )
        temporary.replace(output)
        temporary_index.replace(Path(f"{output}.csi"))
    finally:
        temporary.unlink(missing_ok=True)
        temporary_index.unlink(missing_ok=True)


def write_metric_table(path: Path, rows: Iterable[tuple[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(rows)
    temporary.replace(path)


def parse_junction_file(
    path: Path, sequence_rank: dict[str, int]
) -> Iterator[tuple[tuple[int, int, int, str], JunctionRecord]]:
    previous: tuple[int, int, int, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        if handle.readline().rstrip("\n") != "#" + "\t".join(JUNCTION_COLUMNS):
            raise ValueError(f"unexpected junction header: {path}")
        for number, raw in enumerate(handle, 2):
            values = raw.rstrip("\n").split("\t")
            if len(values) != len(JUNCTION_COLUMNS):
                raise ValueError(f"junction column mismatch: {path}:{number}")
            if values[0] not in sequence_rank:
                raise ValueError(f"junction seqid absent from genome: {values[0]}")
            record = JunctionRecord(
                values[0], *(int(value) for value in values[1:3]), values[3],
                *(int(value) for value in values[4:])
            )
            key = (sequence_rank[record.seqid], record.start, record.end, record.strand)
            if previous is not None and key <= previous:
                raise ValueError(f"junction file is not strictly sorted: {path}:{number}")
            previous = key
            yield key, record


def aggregate_rows(
    sample_paths: Sequence[Path],
    sample_conditions: Sequence[str],
    sequence_rank: dict[str, int],
) -> Iterator[tuple[Any, ...]]:
    streams = [iter(parse_junction_file(path, sequence_rank)) for path in sample_paths]
    heap: list[tuple[tuple[int, int, int, str], int, JunctionRecord]] = []
    for index, stream in enumerate(streams):
        try:
            key, record = next(stream)
        except StopIteration:
            continue
        heapq.heappush(heap, (key, index, record))

    while heap:
        key = heap[0][0]
        observed: dict[int, JunctionRecord] = {}
        while heap and heap[0][0] == key:
            _, index, record = heapq.heappop(heap)
            if index in observed:
                raise ValueError(f"duplicate junction in sample file: {record}")
            observed[index] = record
            try:
                next_key, next_record = next(streams[index])
            except StopIteration:
                continue
            heapq.heappush(heap, (next_key, index, next_record))
        representative = next(iter(observed.values()))
        hq = [observed[index].hq_reads if index in observed else 0 for index in range(len(sample_paths))]
        recurrence = recurrence_metrics(hq, sample_conditions)
        records = list(observed.values())
        yield (
            representative.seqid,
            representative.start,
            representative.end,
            representative.strand,
            sum(record.primary_reads for record in records),
            sum(record.unique_reads for record in records),
            sum(hq),
            max(hq),
            *(recurrence[column] for column in AGGREGATE_COLUMNS[8:17]),
            max(record.max_min_anchor for record in records),
            max(record.max_mapq for record in records),
            ",".join(map(str, hq)),
        )
