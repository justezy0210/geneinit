"""Predictor-independent RNA, protein, repeat, and concordance metrics."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
import csv
import gzip
import json
import math
from pathlib import Path
import re
from typing import Iterable, Iterator, Mapping, Sequence

from geneinit.contracts.streaming import DiskBackedUniqueIds
from geneinit.evidence.junctions import AGGREGATE_COLUMNS
from geneinit.evidence.repeats import REPEAT_SCOPE, is_classified_te, repeat_classification
from geneinit.evidence.metrics import (
    Interval,
    IntervalIndex,
    clip_intervals,
    interval_bp,
    interval_metrics,
    intron_metrics,
    introns_from_blocks,
    merge_intervals,
)
from geneinit.models.candidate import CdsPart, parse_attributes, reverse_complement


DIRECT_GROUP = re.compile(r"^(.+)_(MSTRG\.[0-9]+)$")
NONHC_GROUP = re.compile(r"^(.+)_([0-9]+_g)$")

EVIDENCE_COLUMNS = (
    "candidate_id", "structure_id", "seqid", "strand", "start", "end",
    "predictors", "predictor_count", "cds_part_count", "intron_count",
    "translated_protein_aa", "orf_status", "common_qc_pass",
    "splice_motifs", "major_splice_count", "minor_splice_count",
    "noncanonical_splice_count", "ambiguous_splice_count",
    "rna_best_transcript_id", "rna_exon_overlap_bp", "rna_model_exon_precision",
    "rna_transcript_exon_recall", "rna_exon_f1", "rna_intron_matches",
    "rna_intron_precision", "rna_intron_recall", "rna_intron_f1",
    "rna_hq_read_counts", "rna_sample_hq_ge1_counts",
    "rna_condition_ge2rep_hq_ge1_counts", "rna_condition_ge2rep_hq_ge3_counts",
    "rna_exact_supported_intron_count", "rna_recurrent_supported_intron_count",
    "rna_strong_supported_intron_count", "rna_tier",
    "protein_chain_id", "protein_chain_source", "protein_query_id",
    "protein_target_id", "protein_tied_best_chain_count",
    "protein_tied_best_pairs_json", "protein_chain_complete",
    "protein_cds_overlap_bp",
    "protein_model_precision", "protein_chain_recall", "protein_cds_f1",
    "protein_phase_matches", "protein_phase_comparable_parts",
    "protein_phase_fraction", "protein_intron_matches", "protein_intron_precision",
    "protein_intron_recall", "protein_intron_f1", "protein_start_match",
    "protein_stop_match", "protein_tier", "protein_retained_target_count",
    "protein_seed_target_count", "protein_target_limit", "protein_target_linked",
    "protein_multiplicity_status", "predictor_other_candidate_id",
    "predictor_cds_f1", "predictor_intron_f1", "predictor_exact_cds_chain",
    "predictor_exact_intron_chain", "predictor_exact_phase_chain", "predictor_tier",
    "repeat_scope", "all_repeat_overlap_bp", "all_repeat_cds_fraction",
    "all_repeat_class_bp_json",
    "repeat_overlap_bp", "repeat_cds_fraction", "repeat_top_class",
    "repeat_class_bp_json", "introns_json", "cds_signature",
)


class EvidenceError(RuntimeError):
    """Raised when evidence inputs cannot be reconciled safely."""


@dataclass(frozen=True, slots=True)
class CandidateRecord:
    candidate_id: str
    structure_id: str
    seqid: str
    strand: str
    cds: tuple[CdsPart, ...]
    predictors: tuple[str, ...]
    translated_cds_length: int
    orf_status: str

    @property
    def intervals(self) -> tuple[Interval, ...]:
        return tuple((part.start, part.end) for part in self.cds)

    @property
    def start(self) -> int:
        return self.cds[0].start

    @property
    def end(self) -> int:
        return self.cds[-1].end

    @property
    def introns(self) -> tuple[Interval, ...]:
        return introns_from_blocks(self.intervals)

    @property
    def cds_bp(self) -> int:
        return interval_bp(self.intervals)

    @property
    def protein_aa(self) -> int:
        amino_acids = self.translated_cds_length // 3
        return max(0, amino_acids - int(self.orf_status == "complete"))


@dataclass(frozen=True, slots=True)
class RnaTranscript:
    transcript_id: str
    seqid: str
    strand: str
    exons: tuple[Interval, ...]

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def introns(self) -> tuple[Interval, ...]:
        return introns_from_blocks(self.exons)


@dataclass(frozen=True, slots=True)
class ProteinChain:
    source: str
    seqid: str
    strand: str
    query_id: str
    target_id: str
    cds_parts: tuple[CdsPart, ...]
    introns: tuple[Interval, ...]
    starts: tuple[Interval, ...]
    stops: tuple[Interval, ...]

    @property
    def chain_id(self) -> str:
        return f"{self.source}:{self.query_id}_{self.target_id}"

    @property
    def intervals(self) -> tuple[Interval, ...]:
        return tuple((part.start, part.end) for part in self.cds_parts)

    @property
    def start(self) -> int:
        values = self.intervals + self.introns + self.starts + self.stops
        return min(start for start, _ in values)

    @property
    def end(self) -> int:
        values = self.intervals + self.introns + self.starts + self.stops
        return max(end for _, end in values)

    @property
    def complete(self) -> bool:
        return bool(self.starts and self.stops)


@dataclass(frozen=True, slots=True)
class RepeatFeature:
    seqid: str
    start: int
    end: int
    classification: str


@dataclass(frozen=True, slots=True)
class JunctionAggregate:
    total_hq_reads: int
    sample_count_hq_ge1: int
    condition_count_ge2rep_hq_ge1: int
    condition_count_ge2rep_hq_ge3: int


@dataclass(slots=True)
class _ProteinBuilder:
    source: str
    seqid: str
    strand: str
    query_id: str
    target_id: str
    cds_parts: list[CdsPart] = field(default_factory=list)
    introns: list[Interval] = field(default_factory=list)
    starts: list[Interval] = field(default_factory=list)
    stops: list[Interval] = field(default_factory=list)


def parse_cds_signature(text: str) -> tuple[CdsPart, ...]:
    parts: list[CdsPart] = []
    for token in text.split(","):
        match = re.fullmatch(r"([0-9]+)-([0-9]+):([012])", token)
        if match is None:
            raise EvidenceError(f"invalid CDS signature: {text!r}")
        parts.append(CdsPart(*map(int, match.groups())))
    if not parts:
        raise EvidenceError("candidate CDS signature is empty")
    return tuple(sorted(parts))


def iter_candidates(path: Path) -> Iterator[CandidateRecord]:
    """Stream candidates in the table's canonical coordinate order."""

    observed = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id", "structure_id", "seqid", "strand", "predictors",
            "translated_cds_length", "orf_status", "cds_signature",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise EvidenceError(f"candidate table columns missing: {sorted(missing)}")
        previous_key: tuple[object, ...] | None = None
        with DiskBackedUniqueIds() as unique_ids:
            for number, row in enumerate(reader, 2):
                candidate_id = row["candidate_id"]
                if not candidate_id or not unique_ids.add(candidate_id):
                    raise EvidenceError(
                        f"duplicate/empty candidate ID at {path}:{number}"
                    )
                predictors = tuple(
                    sorted(filter(None, row["predictors"].split(",")))
                )
                if not predictors:
                    raise EvidenceError(
                        f"candidate lacks predictor lineage at {path}:{number}"
                    )
                cds = parse_cds_signature(row["cds_signature"])
                key = (
                    row["seqid"], cds[0].start, cds[-1].end,
                    row["strand"], candidate_id,
                )
                if previous_key is not None and key < previous_key:
                    raise EvidenceError(
                        "candidate table is not in canonical coordinate order "
                        f"at {path}:{number}"
                    )
                previous_key = key
                observed += 1
                yield CandidateRecord(
                    candidate_id=candidate_id,
                    structure_id=row["structure_id"],
                    seqid=row["seqid"],
                    strand=row["strand"],
                    cds=cds,
                    predictors=predictors,
                    translated_cds_length=int(row["translated_cds_length"]),
                    orf_status=row["orf_status"],
                )
    if not observed:
        raise EvidenceError(f"candidate table contains no models: {path}")


def load_candidates(path: Path) -> list[CandidateRecord]:
    """Compatibility loader for callers that explicitly need all candidates."""

    return list(iter_candidates(path))


def iter_candidate_blocks(path: Path) -> Iterator[list[CandidateRecord]]:
    """Yield exact independent span-overlap blocks without a size threshold.

    Candidate tables are coordinate ordered. Once the next model starts after
    the maximum end of the current block, no overlap-based evidence or
    selection edge can cross that boundary. Opposite strands are deliberately
    kept in the same conservative span block so input order is preserved.
    """

    block: list[CandidateRecord] = []
    block_seqid = ""
    block_end = 0
    for candidate in iter_candidates(path):
        if block and (
            candidate.seqid != block_seqid or candidate.start > block_end
        ):
            yield block
            block = []
        if not block:
            block_seqid = candidate.seqid
            block_end = candidate.end
        else:
            block_end = max(block_end, candidate.end)
        block.append(candidate)
    if block:
        yield block


def load_stringtie_transcripts(path: Path) -> list[RnaTranscript]:
    builders: dict[str, tuple[str, str, list[Interval]]] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise EvidenceError(f"{path}:{number} does not have 9 columns")
            if fields[2] != "exon":
                continue
            attrs, _ = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id") or attrs.get("Parent")
            if not transcript_id:
                raise EvidenceError(f"StringTie exon lacks transcript_id: {path}:{number}")
            seqid, strand = fields[0], fields[6]
            if strand == ".":
                # A strandless merged transcript cannot support a strand-aware
                # coding model. StringTie may emit isolated single-exon rows of
                # this kind even for stranded libraries, so ignore them rather
                # than turning them into evidence for both strands.
                continue
            if strand not in {"+", "-"}:
                raise EvidenceError(f"StringTie exon has unknown strand: {path}:{number}")
            interval = (int(fields[3]), int(fields[4]))
            if transcript_id not in builders:
                builders[transcript_id] = (seqid, strand, [interval])
            else:
                prior_seqid, prior_strand, exons = builders[transcript_id]
                if (prior_seqid, prior_strand) != (seqid, strand):
                    raise EvidenceError(f"StringTie transcript ID is reused: {transcript_id}")
                exons.append(interval)
    return sorted(
        (
            RnaTranscript(transcript_id, seqid, strand, merge_intervals(exons))
            for transcript_id, (seqid, strand, exons) in builders.items()
        ),
        key=lambda item: (item.seqid, item.start, item.end, item.strand, item.transcript_id),
    )


def load_trace(path: Path) -> dict[str, tuple[str, int, int]]:
    trace: dict[str, tuple[str, int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise EvidenceError(f"invalid non-HC trace row: {path}:{number}")
            label, seqid, start_text, end_text = fields
            if label in trace:
                raise EvidenceError(f"duplicate non-HC trace label: {label}")
            trace[label] = (seqid, int(start_text), int(end_text))
    if not trace:
        raise EvidenceError(f"non-HC trace contains no regions: {path}")
    return trace


def _parse_protein_group(group: str, source: str) -> tuple[str, str]:
    pattern = DIRECT_GROUP if source == "direct" else NONHC_GROUP
    match = pattern.fullmatch(group)
    if match is None:
        raise EvidenceError(f"invalid {source} ProtHint grp: {group}")
    return match.group(1), match.group(2)


def load_protein_chains(
    path: Path,
    source: str,
    *,
    trace: Mapping[str, tuple[str, int, int]] | None = None,
) -> list[ProteinChain]:
    if source not in {"direct", "nonhc"}:
        raise ValueError(f"unsupported protein-chain source: {source}")
    if source == "nonhc" and trace is None:
        raise EvidenceError("non-HC protein chains require a trace")
    builders: dict[tuple[str, str, str], _ProteinBuilder] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise EvidenceError(f"{path}:{number} does not have 9 columns")
            raw_seqid, feature, strand = fields[0], fields[2].lower(), fields[6]
            if feature not in {"cdspart", "intron", "start", "stop"}:
                continue
            attrs, _ = parse_attributes(fields[8])
            group = attrs.get("grp")
            if not group:
                raise EvidenceError(f"ProtHint row lacks grp: {path}:{number}")
            query_id, target_id = _parse_protein_group(group, source)
            start, end = int(fields[3]), int(fields[4])
            seqid = raw_seqid
            if source == "nonhc":
                if raw_seqid not in trace:
                    raise EvidenceError(f"non-HC region absent from trace: {raw_seqid}")
                seqid, trace_start, trace_end = trace[raw_seqid]
                if start < 1 or end < start or end > trace_end - trace_start + 1:
                    raise EvidenceError(f"non-HC feature exceeds trace: {path}:{number}")
                start += trace_start - 1
                end += trace_start - 1
            key = (raw_seqid, group, strand)
            builder = builders.setdefault(
                key, _ProteinBuilder(source, seqid, strand, query_id, target_id)
            )
            if builder.seqid != seqid:
                raise EvidenceError(f"protein group crosses genome sequences: {group}")
            if feature == "cdspart":
                phase = int(fields[7]) if fields[7] in {"0", "1", "2"} else -1
                builder.cds_parts.append(CdsPart(start, end, phase))
            elif feature == "intron":
                builder.introns.append((start, end))
            elif feature == "start":
                builder.starts.append((start, end))
            else:
                builder.stops.append((start, end))
    chains: list[ProteinChain] = []
    for builder in builders.values():
        if not (builder.cds_parts or builder.introns or builder.starts or builder.stops):
            continue
        chains.append(
            ProteinChain(
                builder.source, builder.seqid, builder.strand,
                builder.query_id, builder.target_id,
                tuple(sorted(set(builder.cds_parts))),
                tuple(sorted(set(builder.introns))),
                tuple(sorted(set(builder.starts))),
                tuple(sorted(set(builder.stops))),
            )
        )
    if not chains:
        raise EvidenceError(f"no retained protein chains: {path}")
    return sorted(chains, key=lambda item: (item.seqid, item.start, item.end, item.strand, item.chain_id))


def load_diamond_targets(
    direct: Path, nonhc: Path
) -> dict[str, frozenset[str]]:
    targets: dict[str, set[str]] = defaultdict(set)
    for path, namespace in ((direct, "direct"), (nonhc, "nonhc")):
        with path.open(encoding="utf-8") as handle:
            for number, raw in enumerate(handle, 1):
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 2 or not all(fields):
                    raise EvidenceError(f"invalid DIAMOND target row: {path}:{number}")
                target, query = fields
                targets[query].add(f"{namespace}:{target}")
    if not targets:
        raise EvidenceError("DIAMOND target files contain no source-target pairs")
    return {query: frozenset(values) for query, values in targets.items()}


def retained_protein_targets(
    chains: Sequence[ProteinChain], raw_targets: Mapping[str, frozenset[str]]
) -> dict[str, frozenset[str]]:
    retained: dict[str, set[str]] = defaultdict(set)
    for chain in chains:
        target = f"{chain.source}:{chain.target_id}"
        if target not in raw_targets.get(chain.query_id, frozenset()):
            raise EvidenceError(
                f"retained protein chain lacks raw DIAMOND pair: {chain.query_id}->{target}"
            )
        retained[chain.query_id].add(target)
    return {query: frozenset(values) for query, values in retained.items()}


def load_repeats(path: Path) -> list[RepeatFeature]:
    repeats: list[RepeatFeature] = []
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise EvidenceError(f"{path}:{number} does not have 9 columns")
            attrs, _ = parse_attributes(fields[8])
            classification = repeat_classification(attrs, fields[2])
            repeats.append(
                RepeatFeature(fields[0], int(fields[3]), int(fields[4]), classification)
            )
    return sorted(repeats, key=lambda item: (item.seqid, item.start, item.end, item.classification))


def load_junction_aggregate(path: Path) -> dict[tuple[str, int, int, str], JunctionAggregate]:
    records: dict[tuple[str, int, int, str], JunctionAggregate] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n")
        expected = "#" + "\t".join(AGGREGATE_COLUMNS)
        if header != expected:
            raise EvidenceError(f"unexpected junction aggregate header: {path}")
        for number, raw in enumerate(handle, 2):
            values = raw.rstrip("\n").split("\t")
            if len(values) != len(AGGREGATE_COLUMNS):
                raise EvidenceError(f"junction column mismatch: {path}:{number}")
            row = dict(zip(AGGREGATE_COLUMNS, values, strict=True))
            key = (row["seqid"], int(row["intron_start"]), int(row["intron_end"]), row["strand"])
            if key in records:
                raise EvidenceError(f"duplicate junction aggregate key: {path}:{number}")
            records[key] = JunctionAggregate(
                total_hq_reads=int(row["total_hq_reads"]),
                sample_count_hq_ge1=int(row["sample_count_hq_ge1"]),
                condition_count_ge2rep_hq_ge1=int(row["condition_count_ge2rep_hq_ge1"]),
                condition_count_ge2rep_hq_ge3=int(row["condition_count_ge2rep_hq_ge3"]),
            )
    return records


def terminal_intervals(candidate: CandidateRecord) -> tuple[tuple[Interval, ...], tuple[Interval, ...]]:
    ordered = sorted(candidate.cds, key=lambda part: part.start, reverse=candidate.strand == "-")
    first, last = ordered[0], ordered[-1]
    if candidate.strand == "+":
        start = (first.start + first.phase, first.start + first.phase + 2)
        stop = (last.end - 2, last.end)
    else:
        start = (first.end - first.phase - 2, first.end - first.phase)
        stop = (last.start, last.start + 2)
    starts = (start,) if start[0] >= first.start and start[1] <= first.end else ()
    stops = (stop,) if stop[0] >= last.start and stop[1] <= last.end else ()
    return starts, stops


def protein_phase_metrics(
    candidate: CandidateRecord, chain: ProteinChain
) -> tuple[int, int, float | None]:
    ordered = sorted(candidate.cds, key=lambda part: part.start, reverse=candidate.strand == "-")
    first_phase = ordered[0].phase
    offsets: list[tuple[CdsPart, int]] = []
    cumulative = 0
    for part in ordered:
        offsets.append((part, cumulative))
        cumulative += part.length
    matches = comparable = 0
    for protein_part in chain.cds_parts:
        if protein_part.phase not in {0, 1, 2}:
            continue
        for candidate_part, offset in offsets:
            if protein_part.start < candidate_part.start or protein_part.end > candidate_part.end:
                continue
            within = (
                protein_part.start - candidate_part.start
                if candidate.strand == "+"
                else candidate_part.end - protein_part.end
            )
            expected = (first_phase - offset - within) % 3
            comparable += 1
            matches += int(protein_part.phase == expected)
            break
    return matches, comparable, matches / comparable if comparable else None


def _none_protein_metrics() -> dict[str, object]:
    return {
        "protein_chain_id": "NA", "protein_chain_source": "NA",
        "protein_query_id": "NA", "protein_target_id": "NA",
        "protein_tied_best_chain_count": 0,
        "protein_tied_best_pairs_json": "[]",
        "protein_chain_complete": "NA", "protein_cds_overlap_bp": 0,
        "protein_model_precision": None, "protein_chain_recall": None,
        "protein_cds_f1": None, "protein_phase_matches": 0,
        "protein_phase_comparable_parts": 0, "protein_phase_fraction": None,
        "protein_intron_matches": 0, "protein_intron_precision": None,
        "protein_intron_recall": None, "protein_intron_f1": None,
        "protein_start_match": 0, "protein_stop_match": 0,
    }


def best_protein_metrics(
    candidate: CandidateRecord, index: IntervalIndex[ProteinChain]
) -> dict[str, object]:
    best_rank: tuple[object, ...] | None = None
    tied: list[tuple[str, ProteinChain, dict[str, object]]] = []
    candidate_starts, candidate_stops = terminal_intervals(candidate)
    for chain in index.query(candidate.seqid, candidate.start, candidate.end, candidate.strand):
        overlap, precision, recall, f1 = interval_metrics(candidate.intervals, chain.intervals)
        matches, intron_precision, intron_recall, intron_f1 = intron_metrics(
            candidate.introns, chain.introns
        )
        phase_matches, phase_total, phase_fraction = protein_phase_metrics(candidate, chain)
        start_match = int(bool(set(candidate_starts) & set(chain.starts)))
        stop_match = int(bool(set(candidate_stops) & set(chain.stops)))
        if overlap == 0 and matches == 0 and start_match == 0 and stop_match == 0:
            continue
        rank = (
            int(matches > 0), matches, intron_f1 if intron_f1 is not None else -1.0,
            recall if recall is not None else -1.0, start_match + stop_match,
            precision if precision is not None else -1.0,
            phase_fraction if phase_fraction is not None else -1.0,
        )
        metrics: dict[str, object] = {
            "protein_chain_id": chain.chain_id, "protein_chain_source": chain.source,
            "protein_query_id": chain.query_id, "protein_target_id": f"{chain.source}:{chain.target_id}",
            "protein_chain_complete": int(chain.complete), "protein_cds_overlap_bp": overlap,
            "protein_model_precision": precision, "protein_chain_recall": recall,
            "protein_cds_f1": f1, "protein_phase_matches": phase_matches,
            "protein_phase_comparable_parts": phase_total,
            "protein_phase_fraction": phase_fraction, "protein_intron_matches": matches,
            "protein_intron_precision": intron_precision,
            "protein_intron_recall": intron_recall, "protein_intron_f1": intron_f1,
            "protein_start_match": start_match, "protein_stop_match": stop_match,
        }
        if best_rank is None or rank > best_rank:
            best_rank = rank
            tied = [(chain.chain_id, chain, metrics)]
        elif rank == best_rank:
            tied.append((chain.chain_id, chain, metrics))
    if not tied:
        return _none_protein_metrics()
    _stable_id, _chain, representative = min(tied, key=lambda item: item[0])
    pairs = sorted(
        {
            (chain.query_id, f"{chain.source}:{chain.target_id}")
            for _stable, chain, _metrics in tied
        }
    )
    representative["protein_tied_best_chain_count"] = len(tied)
    representative["protein_tied_best_pairs_json"] = json.dumps(
        pairs, separators=(",", ":")
    )
    return representative


def protein_tier(metrics: Mapping[str, object], intron_count: int) -> int:
    if metrics["protein_chain_id"] == "NA":
        return 0
    values = [metrics[name] for name in ("protein_model_precision", "protein_chain_recall", "protein_cds_f1")]
    strong = all(isinstance(value, (float, int)) and float(value) >= 0.8 for value in values)
    if intron_count:
        value = metrics["protein_intron_precision"]
        strong = strong and isinstance(value, (float, int)) and float(value) >= 0.8
    terminals = int(metrics["protein_start_match"]) + int(metrics["protein_stop_match"])
    if strong and (metrics["protein_chain_source"] == "direct" or (int(metrics["protein_chain_complete"]) == 1 and terminals == 2)):
        return 4
    if strong:
        return 3
    f1 = metrics["protein_cds_f1"]
    if int(metrics["protein_intron_matches"]) > 0 or (isinstance(f1, (float, int)) and float(f1) >= 0.5):
        return 2
    return 1


def best_rna_metrics(
    candidate: CandidateRecord,
    index: IntervalIndex[RnaTranscript],
    junctions: Mapping[tuple[str, int, int, str], JunctionAggregate],
) -> dict[str, object]:
    best: tuple[tuple[object, ...], str, dict[str, object]] | None = None
    for transcript in index.query(candidate.seqid, candidate.start, candidate.end, candidate.strand):
        relevant_exons = clip_intervals(transcript.exons, candidate.start, candidate.end)
        overlap, precision, recall, f1 = interval_metrics(candidate.intervals, relevant_exons)
        relevant_introns = tuple(
            intron for intron in transcript.introns
            if intron[0] >= candidate.start and intron[1] <= candidate.end
        )
        matches, intron_precision, intron_recall, intron_f1 = intron_metrics(
            candidate.introns, relevant_introns
        )
        if not overlap and not matches:
            continue
        rank = (
            int(matches > 0), matches, intron_f1 if intron_f1 is not None else -1.0,
            f1 if f1 is not None else -1.0, overlap,
        )
        metrics: dict[str, object] = {
            "rna_best_transcript_id": transcript.transcript_id,
            "rna_exon_overlap_bp": overlap, "rna_model_exon_precision": precision,
            "rna_transcript_exon_recall": recall, "rna_exon_f1": f1,
            "rna_intron_matches": matches, "rna_intron_precision": intron_precision,
            "rna_intron_recall": intron_recall, "rna_intron_f1": intron_f1,
        }
        if best is None or rank > best[0] or (rank == best[0] and transcript.transcript_id < best[1]):
            best = rank, transcript.transcript_id, metrics
    structure = best[2] if best else {
        "rna_best_transcript_id": "NA", "rna_exon_overlap_bp": 0,
        "rna_model_exon_precision": None, "rna_transcript_exon_recall": None,
        "rna_exon_f1": None, "rna_intron_matches": 0,
        "rna_intron_precision": None, "rna_intron_recall": None, "rna_intron_f1": None,
    }
    observed = [junctions.get((candidate.seqid, start, end, candidate.strand)) for start, end in candidate.introns]
    hq = [record.total_hq_reads if record else 0 for record in observed]
    samples = [record.sample_count_hq_ge1 if record else 0 for record in observed]
    rep1 = [record.condition_count_ge2rep_hq_ge1 if record else 0 for record in observed]
    rep3 = [record.condition_count_ge2rep_hq_ge3 if record else 0 for record in observed]
    exact = sum(value >= 1 for value in hq)
    recurrent = sum(value >= 2 and sample >= 2 for value, sample in zip(hq, samples, strict=True))
    strong = sum(value >= 3 and replicate >= 1 for value, replicate in zip(hq, rep1, strict=True))
    if not candidate.introns:
        transcript_present = structure["rna_best_transcript_id"] != "NA"
        precision = structure["rna_model_exon_precision"]
        recall = structure["rna_transcript_exon_recall"]
        tier = 2 if transcript_present and isinstance(precision, (float, int)) and isinstance(recall, (float, int)) and min(float(precision), float(recall)) >= 0.8 else int(transcript_present)
    elif all(value >= 1 for value in rep3):
        tier = 4
    elif all(value >= 3 for value in hq) or strong / len(candidate.introns) >= 0.5:
        tier = 3
    elif exact or int(structure["rna_intron_matches"]):
        tier = 2
    else:
        tier = int(structure["rna_best_transcript_id"] != "NA")
    return {
        **structure,
        "rna_hq_read_counts": json.dumps(hq, separators=(",", ":")),
        "rna_sample_hq_ge1_counts": json.dumps(samples, separators=(",", ":")),
        "rna_condition_ge2rep_hq_ge1_counts": json.dumps(rep1, separators=(",", ":")),
        "rna_condition_ge2rep_hq_ge3_counts": json.dumps(rep3, separators=(",", ":")),
        "rna_exact_supported_intron_count": exact,
        "rna_recurrent_supported_intron_count": recurrent,
        "rna_strong_supported_intron_count": strong,
        "rna_tier": tier,
    }


def predictor_metrics(
    candidate: CandidateRecord, index: IntervalIndex[CandidateRecord]
) -> dict[str, object]:
    if len(candidate.predictors) > 1:
        return {
            "predictor_other_candidate_id": candidate.candidate_id,
            "predictor_cds_f1": 1.0, "predictor_intron_f1": 1.0 if candidate.introns else None,
            "predictor_exact_cds_chain": 1, "predictor_exact_intron_chain": int(bool(candidate.introns)),
            "predictor_exact_phase_chain": 1, "predictor_tier": 3,
        }
    best: tuple[tuple[object, ...], str, dict[str, object]] | None = None
    for other in index.query(candidate.seqid, candidate.start, candidate.end, candidate.strand):
        if other.candidate_id == candidate.candidate_id:
            continue
        # A normalized partner can contain both the current predictor and an
        # independent predictor because exact structures were collapsed before
        # concordance scoring.  It remains valid evidence whenever it contributes
        # at least one predictor lineage absent from the current candidate.
        if not (set(other.predictors) - set(candidate.predictors)):
            continue
        _overlap, _precision, _recall, cds_f1 = interval_metrics(candidate.intervals, other.intervals)
        _matches, _ip, _ir, intron_f1 = intron_metrics(candidate.introns, other.introns)
        exact_cds = int(candidate.intervals == other.intervals)
        exact_introns = int(bool(candidate.introns or other.introns) and candidate.introns == other.introns)
        exact_phase = int(candidate.cds == other.cds)
        rank = (exact_phase, exact_cds, exact_introns, intron_f1 or -1.0, cds_f1 or -1.0)
        tier = 3 if exact_phase else 2 if exact_cds or (exact_introns and (cds_f1 or 0) >= 0.8) else 1 if max(cds_f1 or 0, intron_f1 or 0) >= 0.8 else 0
        metrics = {
            "predictor_other_candidate_id": other.candidate_id,
            "predictor_cds_f1": cds_f1, "predictor_intron_f1": intron_f1,
            "predictor_exact_cds_chain": exact_cds,
            "predictor_exact_intron_chain": exact_introns,
            "predictor_exact_phase_chain": exact_phase, "predictor_tier": tier,
        }
        if best is None or rank > best[0] or (rank == best[0] and other.candidate_id < best[1]):
            best = rank, other.candidate_id, metrics
    return best[2] if best else {
        "predictor_other_candidate_id": "NA", "predictor_cds_f1": None,
        "predictor_intron_f1": None, "predictor_exact_cds_chain": 0,
        "predictor_exact_intron_chain": 0, "predictor_exact_phase_chain": 0,
        "predictor_tier": 0,
    }


def repeat_metrics(
    candidate: CandidateRecord, index: IntervalIndex[RepeatFeature]
) -> dict[str, object]:
    all_overlaps: list[Interval] = []
    classified_overlaps: list[Interval] = []
    by_class: dict[str, list[Interval]] = defaultdict(list)
    for repeat in index.query(candidate.seqid, candidate.start, candidate.end):
        for left, right in candidate.intervals:
            start, end = max(left, repeat.start), min(right, repeat.end)
            if start <= end:
                all_overlaps.append((start, end))
                if is_classified_te(repeat.classification):
                    classified_overlaps.append((start, end))
                by_class[repeat.classification].append((start, end))
    all_total = interval_bp(all_overlaps)
    total = interval_bp(classified_overlaps)
    all_class_bp = {name: interval_bp(intervals) for name, intervals in by_class.items()}
    class_bp = {name: bp for name, bp in all_class_bp.items() if is_classified_te(name)}
    top = min(class_bp, key=lambda name: (-class_bp[name], name)) if class_bp else "NA"
    return {
        "repeat_scope": REPEAT_SCOPE,
        "all_repeat_overlap_bp": all_total,
        "all_repeat_cds_fraction": all_total / candidate.cds_bp,
        "all_repeat_class_bp_json": json.dumps(all_class_bp, sort_keys=True, separators=(",", ":")),
        "repeat_overlap_bp": total,
        "repeat_cds_fraction": total / candidate.cds_bp,
        "repeat_top_class": top,
        "repeat_class_bp_json": json.dumps(class_bp, sort_keys=True, separators=(",", ":")),
    }


def splice_metrics(candidate: CandidateRecord, genome: Mapping[str, str]) -> dict[str, object]:
    sequence = genome[candidate.seqid]
    motifs: list[str] = []
    for start, end in candidate.introns:
        if end - start + 1 < 2:
            motif = "SHORT"
        elif candidate.strand == "+":
            motif = sequence[start - 1 : start + 1].upper() + "-" + sequence[end - 2 : end].upper()
        else:
            donor = reverse_complement(sequence[end - 2 : end]).upper()
            acceptor = reverse_complement(sequence[start - 1 : start + 1]).upper()
            motif = donor + "-" + acceptor
        motifs.append(motif)
    major = sum(motif == "GT-AG" for motif in motifs)
    minor = sum(motif in {"GC-AG", "AT-AC"} for motif in motifs)
    ambiguous = sum(not re.fullmatch(r"[ACGT]{2}-[ACGT]{2}", motif) for motif in motifs)
    noncanonical = len(motifs) - major - minor - ambiguous
    return {
        "splice_motifs": json.dumps(motifs, separators=(",", ":")),
        "major_splice_count": major, "minor_splice_count": minor,
        "noncanonical_splice_count": noncanonical,
        "ambiguous_splice_count": ambiguous,
    }


def protein_multiplicity(
    metrics: Mapping[str, object],
    retained: Mapping[str, frozenset[str]],
    seeds: Mapping[str, frozenset[str]],
    expected_copies: int,
) -> dict[str, object]:
    query = str(metrics["protein_query_id"])
    limit = max(1, math.floor(expected_copies * 1.0))
    if query == "NA":
        return {
            "protein_retained_target_count": "NA", "protein_seed_target_count": "NA",
            "protein_target_limit": limit, "protein_target_linked": 0,
            "protein_multiplicity_status": "NO_PROTEIN_CHAIN",
        }
    try:
        pairs = [
            (str(item[0]), str(item[1]))
            for item in json.loads(str(metrics["protein_tied_best_pairs_json"]))
        ]
    except (KeyError, TypeError, ValueError, IndexError, json.JSONDecodeError) as error:
        raise EvidenceError("tied-best protein pairs must be a JSON pair array") from error
    if not pairs:
        raise EvidenceError("protein chain has no tied-best query-target pair")
    statuses: list[str] = []
    retained_counts: list[int] = []
    seed_counts: list[int] = []
    linked_values: list[bool] = []
    for tied_query, target in pairs:
        retained_targets = retained.get(tied_query, frozenset())
        seed_targets = seeds.get(tied_query, frozenset())
        linked = target in retained_targets and target in seed_targets
        retained_counts.append(len(retained_targets))
        seed_counts.append(len(seed_targets))
        linked_values.append(linked)
        if not linked:
            statuses.append("UNRESOLVED_PROTEIN_TARGET")
        elif len(retained_targets) > limit:
            statuses.append("AMBIGUOUS_PROTEIN_MULTI_HIT")
        elif len(seed_targets) > limit:
            statuses.append("UNIQUE_RETAINED_PROTEIN_HIT_SEED_MULTI_HIT_AUDIT")
        else:
            statuses.append("UNIQUE_PROTEIN_HIT")
    priority = {
        "UNIQUE_PROTEIN_HIT": 0,
        "UNIQUE_RETAINED_PROTEIN_HIT_SEED_MULTI_HIT_AUDIT": 1,
        "AMBIGUOUS_PROTEIN_MULTI_HIT": 2,
        "UNRESOLVED_PROTEIN_TARGET": 3,
    }
    status = max(statuses, key=lambda item: priority[item])
    return {
        "protein_retained_target_count": max(retained_counts),
        "protein_seed_target_count": max(seed_counts),
        "protein_target_limit": limit,
        "protein_target_linked": int(all(linked_values)),
        "protein_multiplicity_status": status,
    }


def _score_candidate(
    candidate: CandidateRecord,
    candidate_index: IntervalIndex[CandidateRecord],
    transcript_index: IntervalIndex[RnaTranscript],
    chain_index: IntervalIndex[ProteinChain],
    repeat_index: IntervalIndex[RepeatFeature],
    junctions: Mapping[tuple[str, int, int, str], JunctionAggregate],
    genome: Mapping[str, str],
    retained_targets: Mapping[str, frozenset[str]],
    seed_targets: Mapping[str, frozenset[str]],
    expected_protein_copies: int,
) -> dict[str, object]:
    protein = best_protein_metrics(candidate, chain_index)
    row: dict[str, object] = {
        "candidate_id": candidate.candidate_id,
        "structure_id": candidate.structure_id,
        "seqid": candidate.seqid,
        "strand": candidate.strand,
        "start": candidate.start,
        "end": candidate.end,
        "predictors": ",".join(candidate.predictors),
        "predictor_count": len(candidate.predictors),
        "cds_part_count": len(candidate.cds),
        "intron_count": len(candidate.introns),
        "translated_protein_aa": candidate.protein_aa,
        "orf_status": candidate.orf_status,
        "common_qc_pass": int(candidate.orf_status == "complete"),
    }
    row.update(splice_metrics(candidate, genome))
    row.update(best_rna_metrics(candidate, transcript_index, junctions))
    row.update(protein)
    row["protein_tier"] = protein_tier(protein, len(candidate.introns))
    row.update(
        protein_multiplicity(
            protein, retained_targets, seed_targets, expected_protein_copies
        )
    )
    row.update(predictor_metrics(candidate, candidate_index))
    row.update(repeat_metrics(candidate, repeat_index))
    row["introns_json"] = json.dumps(candidate.introns, separators=(",", ":"))
    row["cds_signature"] = ",".join(
        f"{part.start}-{part.end}:{part.phase}" for part in candidate.cds
    )
    return row


def score_candidate_blocks(
    candidate_blocks: Iterable[Sequence[CandidateRecord]],
    transcripts: Sequence[RnaTranscript],
    chains: Sequence[ProteinChain],
    repeats: Sequence[RepeatFeature],
    junctions: Mapping[tuple[str, int, int, str], JunctionAggregate],
    genome: Mapping[str, str],
    retained_targets: Mapping[str, frozenset[str]],
    seed_targets: Mapping[str, frozenset[str]],
    expected_protein_copies: int,
) -> Iterator[dict[str, object]]:
    """Score independent candidate blocks with one global evidence index."""

    transcript_index = IntervalIndex(
        transcripts, seqid=lambda x: x.seqid, start=lambda x: x.start,
        end=lambda x: x.end, strand=lambda x: x.strand,
    )
    chain_index = IntervalIndex(
        chains, seqid=lambda x: x.seqid, start=lambda x: x.start,
        end=lambda x: x.end, strand=lambda x: x.strand,
    )
    repeat_index = IntervalIndex(
        repeats, seqid=lambda x: x.seqid, start=lambda x: x.start,
        end=lambda x: x.end,
    )
    for candidates in candidate_blocks:
        candidate_index = IntervalIndex(
            candidates, seqid=lambda x: x.seqid, start=lambda x: x.start,
            end=lambda x: x.end, strand=lambda x: x.strand,
        )
        for candidate in candidates:
            yield _score_candidate(
                candidate, candidate_index, transcript_index, chain_index,
                repeat_index, junctions, genome, retained_targets, seed_targets,
                expected_protein_copies,
            )


def score_candidates(
    candidates: Sequence[CandidateRecord],
    transcripts: Sequence[RnaTranscript],
    chains: Sequence[ProteinChain],
    repeats: Sequence[RepeatFeature],
    junctions: Mapping[tuple[str, int, int, str], JunctionAggregate],
    genome: Mapping[str, str],
    retained_targets: Mapping[str, frozenset[str]],
    seed_targets: Mapping[str, frozenset[str]],
    expected_protein_copies: int,
) -> list[dict[str, object]]:
    return list(
        score_candidate_blocks(
            (candidates,), transcripts, chains, repeats, junctions, genome,
            retained_targets, seed_targets, expected_protein_copies,
        )
    )


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format(value, ".10g")
    return str(value)


def write_evidence_table(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=EVIDENCE_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: format_value(row[column]) for column in EVIDENCE_COLUMNS})
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
