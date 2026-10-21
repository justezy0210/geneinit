"""Deterministic provisional gene grouping from selected coding models."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO

from geneinit.utr.core import (
    EMPTY,
    Feature,
    TranscriptModel,
    UtrError,
    _write_table,
    atomic_text,
    format_gff3_attributes,
    load_raw_gff3,
    sha256_file,
)


PAIR_COLUMNS = (
    "selection_locus_id", "left_candidate_id", "right_candidate_id",
    "left_group_id", "right_group_id", "same_provisional_gene",
    "span_overlap_bp", "cds_overlap_bp", "shorter_cds_overlap_fraction",
    "same_coding_frame", "same_frame_cds_overlap_ge3",
    "shared_cds_intron_count", "shared_transcript_intron_count",
    "structural_reasons", "active_fusion_boundary", "action",
)
GROUP_COLUMNS = (
    "provisional_gene_id", "temporary_transcript_id", "candidate_id",
    "canonical_transcript", "canonical_rank", "seqid", "strand",
    "gene_start", "gene_end", "transcript_start", "transcript_end",
    "selection_locus_id", "grouping_basis", "routes", "active_fusion",
    "utr_refine_evidence", "predictors", "source_gene_ids",
    "source_transcript_ids",
)


@dataclass(frozen=True, slots=True)
class PairAssessment:
    reasons: tuple[str, ...]
    active_fusion_boundary: bool
    should_group: bool
    span_overlap_bp: int
    cds_overlap_bp: int
    shorter_cds_overlap_fraction: float
    same_coding_frame: bool
    same_frame_cds_overlap_ge3: bool
    shared_cds_intron_count: int
    shared_transcript_intron_count: int


class DisjointSet:
    def __init__(self, keys: Iterable[str]) -> None:
        self.parent = {key: key for key in keys}

    def find(self, key: str) -> str:
        if self.parent[key] != key:
            self.parent[key] = self.find(self.parent[key])
        return self.parent[key]

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _intervals(features: Iterable[Feature]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((feature.start, feature.end) for feature in features))


def _introns(intervals: Sequence[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (left_end + 1, right_start - 1)
        for (_, left_end), (right_start, _) in zip(intervals, intervals[1:])
        if left_end < right_start - 1
    )


def _interval_bp(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in intervals)


def _overlap_bp(
    left: Sequence[tuple[int, int]], right: Sequence[tuple[int, int]]
) -> int:
    left_index = right_index = shared = 0
    while left_index < len(left) and right_index < len(right):
        left_start, left_end = left[left_index]
        right_start, right_end = right[right_index]
        shared += max(0, min(left_end, right_end) - max(left_start, right_start) + 1)
        if left_end < right_end:
            left_index += 1
        else:
            right_index += 1
    return shared


def _coding_anchor(model: TranscriptModel) -> int:
    cds = sorted(model.cds, key=lambda feature: (feature.start, feature.end))
    first = cds[0] if model.strand == "+" else cds[-1]
    phase = int(first.phase)
    return first.start + phase if model.strand == "+" else first.end - phase


def _mstrg_gene(value: str) -> str:
    parts = value.split(".")
    if len(parts) >= 2 and parts[0] == "MSTRG" and parts[1].isdigit():
        return ".".join(parts[:2])
    return ""


def _lineage_sets(
    candidate_table: Path, source_models: Path
) -> dict[str, frozenset[tuple[str, str]]]:
    with candidate_table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_id", "structure_id"}
        if required - set(reader.fieldnames or ()):
            raise UtrError("candidate table lacks candidate_id/structure_id")
        structure_to_candidate = {row["structure_id"]: row["candidate_id"] for row in reader}
    lineages: defaultdict[str, set[tuple[str, str]]] = defaultdict(set)
    with source_models.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"structure_id", "predictor", "source_gene_id"}
        if required - set(reader.fieldnames or ()):
            raise UtrError("source model table lacks structural lineage columns")
        for row in reader:
            candidate_id = structure_to_candidate.get(row["structure_id"])
            gene_id = row["source_gene_id"]
            if candidate_id and gene_id not in EMPTY:
                lineages[candidate_id].add((row["predictor"], gene_id))
    return {key: frozenset(values) for key, values in lineages.items()}


def assess_pair(
    left: TranscriptModel,
    right: TranscriptModel,
    left_metadata: Mapping[str, str],
    right_metadata: Mapping[str, str],
    lineages: Mapping[str, frozenset[tuple[str, str]]],
) -> PairAssessment:
    if (left.seqid, left.strand) != (right.seqid, right.strand):
        raise UtrError("gene grouping pair crosses reference sequence or strand")
    span_overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    left_cds, right_cds = _intervals(left.cds), _intervals(right.cds)
    cds_overlap = _overlap_bp(left_cds, right_cds)
    shorter = min(_interval_bp(left_cds), _interval_bp(right_cds))
    reciprocal = cds_overlap / shorter if shorter else 0.0
    left_cds_introns, right_cds_introns = set(_introns(left_cds)), set(_introns(right_cds))
    left_tx_introns = set(_introns(_intervals(left.exons)))
    right_tx_introns = set(_introns(_intervals(right.exons)))
    shared_cds = left_cds_introns & right_cds_introns
    shared_tx = left_tx_introns & right_tx_introns
    same_frame = (_coding_anchor(left) - _coding_anchor(right)) % 3 == 0
    same_frame_overlap = same_frame and cds_overlap >= 3
    reasons: set[str] = set()
    # Predictor/RNA labels may connect alternatives only when the transcript
    # spans actually overlap. This prevents one erroneous fusion label from
    # joining two non-overlapping genes.
    if span_overlap and lineages.get(left.candidate_id, frozenset()) & lineages.get(
        right.candidate_id, frozenset()
    ):
        reasons.add("SAME_PREDICTOR_GENE")
    left_mstrg = _mstrg_gene(left.utr_refine_evidence)
    right_mstrg = _mstrg_gene(right.utr_refine_evidence)
    if span_overlap and left_mstrg and left_mstrg == right_mstrg:
        reasons.add("SAME_STRINGTIE_GENE")
    if left_cds_introns and left_cds_introns == right_cds_introns:
        reasons.add("SAME_CDS_INTRON_CHAIN")
    if shared_cds:
        reasons.add("SHARED_CDS_INTRON")
    if shared_tx:
        reasons.add("SHARED_TRANSCRIPT_INTRON")
    if len(left_cds) == len(right_cds) == 1 and reciprocal >= 0.8 and same_frame:
        reasons.add("SINGLE_CDS_EXON_SAME_FRAME_RMO80")
    boundary = left_metadata["active_fusion"] == "1" or right_metadata["active_fusion"] == "1"
    return PairAssessment(
        reasons=tuple(sorted(reasons)),
        active_fusion_boundary=boundary,
        should_group=bool(reasons) and not boundary,
        span_overlap_bp=span_overlap,
        cds_overlap_bp=cds_overlap,
        shorter_cds_overlap_fraction=reciprocal,
        same_coding_frame=same_frame,
        same_frame_cds_overlap_ge3=same_frame_overlap,
        shared_cds_intron_count=len(shared_cds),
        shared_transcript_intron_count=len(shared_tx),
    )


def _canonical_key(model: TranscriptModel, metadata: Mapping[str, str]) -> tuple[object, ...]:
    routes = metadata["routes"]
    rna = "RNA" in routes
    protein = "PROTEIN" in routes
    cds_bp = _interval_bp(_intervals(model.cds))
    return (-int(rna), -int(protein), -cds_bp, model.candidate_id)


def _stable_group_id(candidates: Iterable[str]) -> str:
    payload = json.dumps(sorted(candidates), separators=(",", ":"))
    return f"GGRP-{hashlib.sha256(payload.encode()).hexdigest()[:24]}"


def _write_gff_row(
    handle: TextIO,
    seqid: str,
    source: str,
    feature_type: str,
    start: int,
    end: int,
    strand: str,
    phase: str,
    attributes: Iterable[tuple[str, str]],
) -> None:
    handle.write(
        f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{strand}\t{phase}\t"
        f"{format_gff3_attributes(attributes)}\n"
    )


def group_coding_models(
    coding_gff3: Path,
    coding_mapping: Path,
    candidate_table: Path,
    source_models: Path,
    reference_sequence_ids: Sequence[str],
    output_gff3: Path,
    output_mapping: Path,
    pair_audit: Path,
    summary_path: Path,
    complete_path: Path,
) -> None:
    order, by_transcript = load_raw_gff3(coding_gff3)
    by_candidate = {model.candidate_id: model for model in by_transcript.values()}
    if len(by_candidate) != len(by_transcript):
        raise UtrError("coding model GFF3 contains duplicate candidate IDs")
    with coding_mapping.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id", "temporary_transcript_id", "selection_locus_id", "routes",
            "active_fusion", "predictors", "source_gene_ids", "source_transcript_ids",
            "utr_refine_evidence",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise UtrError(f"coding mapping lacks columns: {sorted(missing)}")
        metadata = {row["candidate_id"]: row for row in reader}
    if set(metadata) != set(by_candidate):
        raise UtrError("coding GFF3 and mapping candidate sets differ")
    lineages = _lineage_sets(candidate_table, source_models)
    loci: defaultdict[str, list[str]] = defaultdict(list)
    for candidate_id in by_candidate:
        loci[metadata[candidate_id]["selection_locus_id"]].append(candidate_id)

    disjoint = DisjointSet(by_candidate)
    assessments: dict[tuple[str, str], PairAssessment] = {}
    direct_reasons: defaultdict[str, set[str]] = defaultdict(set)
    for locus_id in sorted(loci):
        candidates = sorted(loci[locus_id])
        for left_id, right_id in combinations(candidates, 2):
            assessment = assess_pair(
                by_candidate[left_id], by_candidate[right_id], metadata[left_id],
                metadata[right_id], lineages,
            )
            assessments[(left_id, right_id)] = assessment
            if assessment.should_group:
                disjoint.union(left_id, right_id)
                direct_reasons[left_id].update(assessment.reasons)
                direct_reasons[right_id].update(assessment.reasons)

    components: defaultdict[str, list[str]] = defaultdict(list)
    for candidate_id in by_candidate:
        components[disjoint.find(candidate_id)].append(candidate_id)
    sequence_rank = {seqid: index for index, seqid in enumerate(reference_sequence_ids)}
    groups = sorted(
        components.values(),
        key=lambda values: (
            sequence_rank[by_candidate[values[0]].seqid],
            min(by_candidate[value].start for value in values),
            max(by_candidate[value].end for value in values),
            by_candidate[values[0]].strand,
            min(values),
        ),
    )
    assignment: dict[str, tuple[str, int]] = {}
    group_rows: list[dict[str, object]] = []
    with atomic_text(output_gff3) as handle:
        handle.write("##gff-version 3\n")
        for values in groups:
            group_id = _stable_group_id(values)
            ranked = sorted(values, key=lambda value: _canonical_key(by_candidate[value], metadata[value]))
            models = [by_candidate[value] for value in ranked]
            gene_start = min(model.start for model in models)
            gene_end = max(model.end for model in models)
            seqid, strand = models[0].seqid, models[0].strand
            basis = sorted({reason for value in values for reason in direct_reasons[value]}) or ["SINGLETON"]
            _write_gff_row(
                handle, seqid, "GeneInIt", "gene", gene_start, gene_end, strand, ".",
                (("ID", group_id), ("transcript_count", str(len(models))),
                 ("canonical_transcript", models[0].transcript_id),
                 ("grouping_basis", ",".join(basis))),
            )
            for rank, model in enumerate(models, 1):
                assignment[model.candidate_id] = (group_id, rank)
                attrs = [
                    ("ID", model.transcript_id), ("Parent", group_id),
                    ("candidate_id", model.candidate_id), ("canonical", "true" if rank == 1 else "false"),
                    ("selection_locus_id", metadata[model.candidate_id]["selection_locus_id"]),
                ]
                if model.utr_refine_evidence not in EMPTY:
                    attrs.append(("utr_refine_evidence", model.utr_refine_evidence))
                _write_gff_row(
                    handle, model.seqid, model.source, "mRNA", model.start, model.end,
                    model.strand, ".", attrs,
                )
                for feature_type, features in (
                    ("exon", model.exons), ("CDS", model.cds),
                    ("five_prime_UTR", model.five_prime_utrs),
                    ("three_prime_UTR", model.three_prime_utrs),
                ):
                    for index, feature in enumerate(sorted(features, key=lambda item: (item.start, item.end, item.phase)), 1):
                        _write_gff_row(
                            handle, model.seqid, feature.source, feature_type,
                            feature.start, feature.end, model.strand,
                            feature.phase if feature_type == "CDS" else ".",
                            (("ID", f"{feature_type}-{model.transcript_id}-{index:03d}"),
                             ("Parent", model.transcript_id)),
                        )
                source = metadata[model.candidate_id]
                group_rows.append(
                    {
                        "provisional_gene_id": group_id,
                        "temporary_transcript_id": model.transcript_id,
                        "candidate_id": model.candidate_id,
                        "canonical_transcript": int(rank == 1),
                        "canonical_rank": rank,
                        "seqid": model.seqid,
                        "strand": model.strand,
                        "gene_start": gene_start,
                        "gene_end": gene_end,
                        "transcript_start": model.start,
                        "transcript_end": model.end,
                        "selection_locus_id": source["selection_locus_id"],
                        "grouping_basis": ",".join(sorted(direct_reasons[model.candidate_id])) or "SINGLETON",
                        "routes": source["routes"],
                        "active_fusion": source["active_fusion"],
                        "utr_refine_evidence": source["utr_refine_evidence"],
                        "predictors": source["predictors"],
                        "source_gene_ids": source["source_gene_ids"],
                        "source_transcript_ids": source["source_transcript_ids"],
                    }
                )
    _write_table(output_mapping, GROUP_COLUMNS, group_rows)

    pair_rows: list[dict[str, object]] = []
    for (left_id, right_id), assessment in sorted(assessments.items()):
        left_group, right_group = assignment[left_id][0], assignment[right_id][0]
        same_gene = left_group == right_group
        if assessment.active_fusion_boundary:
            action = "KEEP_SEPARATE_ACTIVE_FUSION"
        elif assessment.should_group:
            action = "GROUP_DIRECT"
        elif same_gene:
            action = "GROUP_TRANSITIVE"
        elif assessment.same_frame_cds_overlap_ge3:
            action = "POSTQC_SAME_FRAME_OVERLAP"
        else:
            action = "KEEP_SEPARATE_NO_LINK"
        pair_rows.append(
            {
                "selection_locus_id": metadata[left_id]["selection_locus_id"],
                "left_candidate_id": left_id, "right_candidate_id": right_id,
                "left_group_id": left_group, "right_group_id": right_group,
                "same_provisional_gene": int(same_gene),
                "span_overlap_bp": assessment.span_overlap_bp,
                "cds_overlap_bp": assessment.cds_overlap_bp,
                "shorter_cds_overlap_fraction": f"{assessment.shorter_cds_overlap_fraction:.8f}",
                "same_coding_frame": int(assessment.same_coding_frame),
                "same_frame_cds_overlap_ge3": int(assessment.same_frame_cds_overlap_ge3),
                "shared_cds_intron_count": assessment.shared_cds_intron_count,
                "shared_transcript_intron_count": assessment.shared_transcript_intron_count,
                "structural_reasons": ",".join(assessment.reasons) or ".",
                "active_fusion_boundary": int(assessment.active_fusion_boundary),
                "action": action,
            }
        )
    _write_table(pair_audit, PAIR_COLUMNS, pair_rows)
    action_counts = Counter(row["action"] for row in pair_rows)
    summary_rows = [
        ("status", "PASS_MIGRATION_PENDING"),
        ("input_transcript_count", len(by_candidate)),
        ("provisional_gene_count", len(groups)),
        ("multi_isoform_gene_count", sum(len(values) > 1 for values in groups)),
        ("pair_count", len(pair_rows)),
        ("postqc_same_frame_overlap_pair_count", action_counts["POSTQC_SAME_FRAME_OVERLAP"]),
    ]
    summary_rows.extend((f"action::{key}", value) for key, value in sorted(action_counts.items()))
    _write_table(summary_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in summary_rows))
    completion = (
        ("schema_version", 1), ("status", "PASS_MIGRATION_PENDING"),
        ("input_transcript_count", len(by_candidate)), ("provisional_gene_count", len(groups)),
        ("grouped_gff3_sha256", sha256_file(output_gff3)),
        ("group_mapping_sha256", sha256_file(output_mapping)),
        ("pair_audit_sha256", sha256_file(pair_audit)),
        ("summary_sha256", sha256_file(summary_path)),
    )
    _write_table(complete_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in completion))
