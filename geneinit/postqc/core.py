"""Transcript-first structural correction without comparative-genomics input.

Gene containers are treated as provisional. Coding-compatible transcript nodes
are rebuilt into connected components, supported bridge candidates are
recovered, and unsupported articulation transcripts are retired. The final
fixed-point gate requires a one-to-one correspondence between coding components
and gene containers.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
import hashlib
from itertools import zip_longest
import json
import math
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence, TextIO

from geneinit.contracts.streaming import DiskBackedUniqueIds
from geneinit.evidence.te_classification import TeClassifier, confirmed_only, model_identity, protein_sha256
from geneinit.evidence.repeats import validate_repeat_scope
from geneinit.evidence.components import JunctionAggregate, load_junction_aggregate
from geneinit.evidence.metrics import IntervalIndex
from geneinit.models.candidate import CdsPart, audit_orf, read_fasta
from geneinit.utr.core import (
    EMPTY,
    Feature,
    TranscriptModel,
    _write_table,
    atomic_text,
    format_gff3_attributes,
    load_raw_gff3,
    parse_cds_signature,
    sha256_file,
)


EDGE_COLUMNS = (
    "left_transcript_id", "right_transcript_id", "left_candidate_id",
    "right_candidate_id", "left_input_gene_id", "right_input_gene_id",
    "seqid", "strand", "span_overlap_bp", "cds_overlap_bp",
    "same_frame_cds_overlap_bp", "shared_cds_intron_count", "coding_link",
    "active_fusion_boundary", "present_initial", "present_final",
    "introduced_by_promotion", "crosses_input_gene_boundary", "disposition",
)
TRANSCRIPT_COLUMNS = (
    "transcript_id", "candidate_id", "source_status", "input_gene_id",
    "postqc_gene_id", "canonical_transcript", "canonical_rank", "seqid",
    "strand", "start", "end", "selection_locus_id", "routes",
    "external_support", "active_fusion", "operation", "reason",
)
GENE_COLUMNS = (
    "postqc_gene_id", "seqid", "strand", "start", "end", "transcript_count",
    "input_gene_count", "input_gene_ids", "promoted_transcript_count",
    "operation",
)
PROMOTION_COLUMNS = (
    "candidate_id", "transcript_id", "seqid", "strand", "start", "end",
    "initial_component_count", "hit_transcript_count", "hit_transcripts",
    "rna_tier", "protein_tier", "status",
)
BRIDGE_REVIEW_COLUMNS = (
    "candidate_id", "transcript_id", "seqid", "strand", "start", "end",
    "intron_start", "intron_end", "intron_bp", "left_transcript_id",
    "right_transcript_id", "total_hq_reads", "sample_count_hq_ge1",
    "condition_count_ge2rep_hq_ge1", "condition_count_ge2rep_hq_ge3",
    "automatic_max_intron_bp", "review_max_intron_bp", "threshold_class",
    "status", "required_evidence",
)
RETIREMENT_COLUMNS = (
    "transcript_id", "candidate_id", "input_gene_id", "seqid", "strand",
    "neighbor_component_count", "protected_neighbor_component_count",
    "status", "reason",
)
ITERATION_COLUMNS = (
    "iteration", "state_before_sha256", "state_after_sha256",
    "supported_bridge_promotions", "strong_rna_bridge_promotions",
    "bridge_review_candidates", "connector_retirements", "changed",
    "converged",
)


class PostQcError(RuntimeError):
    """Raised when post-QC correction cannot satisfy structural invariants."""


@dataclass(slots=True)
class Node:
    model: TranscriptModel
    input_gene_id: str
    source_status: str
    selection_locus_id: str
    routes: str
    external_support: bool
    rna_tier: int
    protein_tier: int
    active_fusion: bool = False
    retired: bool = False
    postqc_gene_id: str = ""
    canonical_rank: int = 0
    operation: str = ""
    reason: str = ""

    @property
    def transcript_id(self) -> str:
        return self.model.transcript_id

    @property
    def candidate_id(self) -> str:
        return self.model.candidate_id

    @property
    def start(self) -> int:
        return self.model.start

    @property
    def end(self) -> int:
        return self.model.end


@dataclass(frozen=True, slots=True)
class CodingEdge:
    left: str
    right: str
    span_overlap_bp: int
    cds_overlap_bp: int
    same_frame_cds_overlap_bp: int
    shared_cds_intron_count: int
    active_fusion_boundary: bool = False

    @property
    def coding_link(self) -> bool:
        return not self.active_fusion_boundary and (
            self.same_frame_cds_overlap_bp >= 3
            or self.shared_cds_intron_count > 0
        )

    def linked(self, minimum_same_frame_bp: int) -> bool:
        return not self.active_fusion_boundary and (
            self.same_frame_cds_overlap_bp >= minimum_same_frame_bp
            or self.shared_cds_intron_count > 0
        )


@dataclass(frozen=True, slots=True)
class RnaBridgeThresholds:
    requested_mode: str
    effective_mode: str
    training_intron_count: int
    automatic_quantile: float
    automatic_quantile_bp: int
    review_quantile: float
    review_quantile_bp: int
    automatic_max_intron_bp: int
    review_max_intron_bp: int


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


def _cds_introns(model: TranscriptModel) -> tuple[tuple[int, int], ...]:
    intervals = _intervals(model.cds)
    return tuple(
        (left_end + 1, right_start - 1)
        for (_, left_end), (right_start, _) in zip(intervals, intervals[1:])
        if left_end < right_start - 1
    )


def _exon_introns(model: TranscriptModel) -> tuple[tuple[int, int], ...]:
    intervals = _intervals(model.exons)
    return tuple(
        (left_end + 1, right_start - 1)
        for (_, left_end), (right_start, _) in zip(intervals, intervals[1:])
        if left_end < right_start - 1
    )


def _nearest_rank(values: Sequence[int], quantile: float) -> int:
    if not values:
        return 0
    if not 0 < quantile <= 1:
        raise PostQcError("intron length quantiles must be in (0, 1]")
    ordered = sorted(values)
    return ordered[max(0, math.ceil(quantile * len(ordered)) - 1)]


def derive_rna_bridge_thresholds(
    trusted_intron_lengths: Sequence[int],
    length_policy: Mapping[str, object],
) -> RnaBridgeThresholds:
    mode = str(length_policy.get("mode", ""))
    lengths = [int(value) for value in trusted_intron_lengths]
    if any(value < 1 for value in lengths):
        raise PostQcError("trusted intron lengths must be positive")
    if mode == "fixed":
        automatic = int(length_policy["fixed_automatic_max_intron_bp"])
        review = int(length_policy["fixed_review_max_intron_bp"])
        if automatic < 1 or review < automatic:
            raise PostQcError("fixed RNA bridge intron thresholds are invalid")
        return RnaBridgeThresholds(
            requested_mode=mode, effective_mode=mode,
            training_intron_count=len(lengths), automatic_quantile=0,
            automatic_quantile_bp=0, review_quantile=0,
            review_quantile_bp=0, automatic_max_intron_bp=automatic,
            review_max_intron_bp=review,
        )
    if mode != "adaptive":
        raise PostQcError("RNA bridge intron length policy mode must be adaptive or fixed")

    minimum_count = int(length_policy["minimum_training_intron_count"])
    automatic_quantile = float(length_policy["automatic_quantile"])
    review_quantile = float(length_policy["review_quantile"])
    automatic_floor = int(length_policy["automatic_floor_bp"])
    automatic_ceiling = int(length_policy["automatic_ceiling_bp"])
    review_multiplier = float(length_policy["review_multiplier"])
    review_ceiling = int(length_policy["review_ceiling_bp"])
    fallback_automatic = int(length_policy["fallback_automatic_max_intron_bp"])
    fallback_review = int(length_policy["fallback_review_max_intron_bp"])
    if minimum_count < 1:
        raise PostQcError("minimum adaptive intron training count must be positive")
    if not 0 < automatic_quantile <= review_quantile <= 1:
        raise PostQcError("adaptive intron quantiles must satisfy 0 < automatic <= review <= 1")
    if not 1 <= automatic_floor <= automatic_ceiling <= review_ceiling:
        raise PostQcError("adaptive RNA bridge intron ceilings/floors are invalid")
    if review_multiplier < 1:
        raise PostQcError("adaptive RNA bridge review multiplier must be at least one")
    if fallback_automatic < 1 or fallback_review < fallback_automatic:
        raise PostQcError("adaptive RNA bridge fallback thresholds are invalid")

    automatic_quantile_bp = _nearest_rank(lengths, automatic_quantile)
    review_quantile_bp = _nearest_rank(lengths, review_quantile)
    if len(lengths) < minimum_count:
        return RnaBridgeThresholds(
            requested_mode=mode, effective_mode="adaptive_fallback",
            training_intron_count=len(lengths),
            automatic_quantile=automatic_quantile,
            automatic_quantile_bp=automatic_quantile_bp,
            review_quantile=review_quantile,
            review_quantile_bp=review_quantile_bp,
            automatic_max_intron_bp=fallback_automatic,
            review_max_intron_bp=fallback_review,
        )
    automatic = min(
        automatic_ceiling, max(automatic_floor, automatic_quantile_bp)
    )
    review = min(
        review_ceiling,
        max(automatic, math.ceil(automatic * review_multiplier), review_quantile_bp),
    )
    return RnaBridgeThresholds(
        requested_mode=mode, effective_mode=mode,
        training_intron_count=len(lengths),
        automatic_quantile=automatic_quantile,
        automatic_quantile_bp=automatic_quantile_bp,
        review_quantile=review_quantile,
        review_quantile_bp=review_quantile_bp,
        automatic_max_intron_bp=automatic,
        review_max_intron_bp=review,
    )


def trusted_intron_lengths(
    nodes: Mapping[str, Node],
    evidence_by_candidate: Mapping[str, Mapping[str, str]],
    junctions: Mapping[tuple[str, int, int, str], JunctionAggregate],
    minimum_hq_reads: int,
    minimum_replicated_conditions: int,
) -> list[int]:
    trusted: set[tuple[str, int, int, str]] = set()
    for node in nodes.values():
        if node.retired or node.source_status != "RELEASED":
            continue
        evidence = evidence_by_candidate.get(node.candidate_id)
        if (
            evidence is None
            or evidence.get("orf_status") != "complete"
            or int(evidence.get("common_qc_pass") or 0) != 1
        ):
            continue
        exon_introns = set(_exon_introns(node.model))
        for intron_start, intron_end in _cds_introns(node.model):
            if (intron_start, intron_end) not in exon_introns:
                continue
            key = (node.model.seqid, intron_start, intron_end, node.model.strand)
            aggregate = junctions.get(key)
            if aggregate is None or (
                aggregate.total_hq_reads < minimum_hq_reads
                or aggregate.condition_count_ge2rep_hq_ge1
                < minimum_replicated_conditions
            ):
                continue
            trusted.add(key)
    return sorted(end - start + 1 for _, start, end, _ in trusted)


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


def same_frame_overlap_bp(left: TranscriptModel, right: TranscriptModel) -> int:
    """Count overlapping CDS bases whose local GFF3 coding frame agrees."""

    if (left.seqid, left.strand) != (right.seqid, right.strand):
        return 0
    shared = 0
    for left_part in left.cds:
        left_anchor = (
            left_part.start + int(left_part.phase)
            if left.strand == "+"
            else left_part.end - int(left_part.phase)
        )
        for right_part in right.cds:
            start = max(left_part.start, right_part.start)
            end = min(left_part.end, right_part.end)
            if start > end:
                continue
            right_anchor = (
                right_part.start + int(right_part.phase)
                if right.strand == "+"
                else right_part.end - int(right_part.phase)
            )
            if (left_anchor - right_anchor) % 3 == 0:
                shared += end - start + 1
    return shared


def assess_edge(left: Node, right: Node) -> CodingEdge:
    if (left.model.seqid, left.model.strand) != (right.model.seqid, right.model.strand):
        raise PostQcError("post-QC edge crosses reference sequence or strand")
    return CodingEdge(
        left=left.transcript_id,
        right=right.transcript_id,
        span_overlap_bp=max(0, min(left.end, right.end) - max(left.start, right.start) + 1),
        cds_overlap_bp=_overlap_bp(_intervals(left.model.cds), _intervals(right.model.cds)),
        same_frame_cds_overlap_bp=same_frame_overlap_bp(left.model, right.model),
        shared_cds_intron_count=len(set(_cds_introns(left.model)) & set(_cds_introns(right.model))),
        active_fusion_boundary=left.active_fusion or right.active_fusion,
    )


def build_edges(nodes: Mapping[str, Node]) -> list[CodingEdge]:
    by_axis: defaultdict[tuple[str, str], list[Node]] = defaultdict(list)
    for node in nodes.values():
        if not node.retired:
            by_axis[(node.model.seqid, node.model.strand)].append(node)
    edges: list[CodingEdge] = []
    for axis in sorted(by_axis):
        ordered = sorted(by_axis[axis], key=lambda item: (item.start, item.end, item.transcript_id))
        for index, left in enumerate(ordered):
            for right_index in range(index + 1, len(ordered)):
                right = ordered[right_index]
                if right.start > left.end:
                    break
                edge = assess_edge(left, right)
                if edge.span_overlap_bp:
                    edges.append(edge)
    return edges


def coding_components(
    nodes: Mapping[str, Node], edges: Sequence[CodingEdge], minimum_same_frame_bp: int = 3
) -> tuple[dict[str, list[str]], dict[str, str]]:
    active = sorted(key for key, node in nodes.items() if not node.retired)
    dsu = DisjointSet(active)
    for edge in edges:
        if edge.linked(minimum_same_frame_bp):
            dsu.union(edge.left, edge.right)
    components: defaultdict[str, list[str]] = defaultdict(list)
    for key in active:
        components[dsu.find(key)].append(key)
    ordered = {
        root: sorted(values) for root, values in sorted(
            components.items(), key=lambda item: min(item[1])
        )
    }
    transcript_to_root = {
        transcript_id: root for root, values in ordered.items() for transcript_id in values
    }
    return ordered, transcript_to_root


def _read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise PostQcError(f"{path} lacks columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def _candidate_model(row: Mapping[str, str]) -> TranscriptModel:
    candidate_id = row["candidate_id"]
    parts = parse_cds_signature(row["cds_signature"])
    transcript_id = f"GIPROM-{candidate_id.removeprefix('GCAN-')}"
    cds = [
        Feature(row["seqid"], "GeneInIt-postQC", "CDS", part.start, part.end, row["strand"], str(part.phase), {})
        for part in parts
    ]
    exons = [
        Feature(row["seqid"], "GeneInIt-postQC", "exon", part.start, part.end, row["strand"], ".", {})
        for part in parts
    ]
    return TranscriptModel(
        gene_id=".", transcript_id=transcript_id, candidate_id=candidate_id,
        seqid=row["seqid"], strand=row["strand"], source="GeneInIt-postQC",
        exons=exons, cds=cds,
    )


def _signature(model: TranscriptModel) -> tuple[object, ...]:
    return (
        model.seqid, model.strand,
        tuple(sorted((feature.start, feature.end, feature.phase) for feature in model.cds)),
    )


def _eligible_bridge(
    candidate: Mapping[str, str], evidence: Mapping[str, str], audit: Mapping[str, str],
    bridge_policy: Mapping[str, object],
) -> bool:
    exclude_active_fusion = bool(bridge_policy["exclude_active_fusion"])
    if confirmed_only(bridge_policy):
        if audit.get("te_policy_id") != bridge_policy["te_policy_id"]:
            raise PostQcError("confirmed-only bridge requires a current TE selection audit")
        te_allowed = audit.get("candidate_pool") == "GENERAL"
    else:
        te_allowed = float(evidence["repeat_cds_fraction"] or 0) < float(bridge_policy["maximum_repeat_cds_fraction"])
    return (
        audit["selected"] == "0"
        and audit["hard_qc_pass"] == "1"
        and (not exclude_active_fusion or audit["active_fusion"] == "0")
        and evidence["orf_status"] == "complete"
        and evidence["common_qc_pass"] == "1"
        and te_allowed
        and (
            int(evidence["rna_tier"] or 0) >= int(bridge_policy["minimum_rna_tier"])
            or int(evidence["protein_tier"] or 0)
            >= int(bridge_policy["minimum_protein_tier"])
        )
        and candidate["candidate_id"] == evidence["candidate_id"] == audit["candidate_id"]
    )


def _iter_table_rows(
    path: Path, required: set[str]
) -> Iterable[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise PostQcError(f"{path} lacks columns: {sorted(missing)}")
        for row in reader:
            yield {
                key: (value or "").strip() for key, value in row.items()
            }


def load_postqc_candidate_inputs(
    candidate_table: Path,
    component_evidence: Path,
    selection_audit: Path,
    selected_candidate_ids: set[str],
    bridge_policy: Mapping[str, object],
) -> tuple[
    list[dict[str, str]],
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
]:
    """Lockstep-stream large tables and retain only post-QC working rows."""

    candidate_rows = _iter_table_rows(
        candidate_table,
        {"candidate_id", "seqid", "strand", "start", "end", "cds_signature"},
    )
    evidence_rows = _iter_table_rows(
        component_evidence,
        {
            "candidate_id", "orf_status", "common_qc_pass", "rna_tier",
            "protein_tier", "repeat_cds_fraction", "repeat_scope",
        },
    )
    audit_rows = _iter_table_rows(
        selection_audit,
        {
            "candidate_id", "locus_id", "selected", "hard_qc_pass",
            "active_fusion", "routes",
        },
    )
    retained_candidates: list[dict[str, str]] = []
    evidence_by_candidate: dict[str, dict[str, str]] = {}
    audit_by_candidate: dict[str, dict[str, str]] = {}
    sentinel = object()
    previous_key: tuple[object, ...] | None = None
    with DiskBackedUniqueIds() as unique_ids:
        for number, triplet in enumerate(
            zip_longest(
                candidate_rows, evidence_rows, audit_rows, fillvalue=sentinel
            ),
            2,
        ):
            candidate, evidence, audit = triplet
            if sentinel in triplet:
                raise PostQcError(
                    "candidate, component-evidence, and selection-audit row "
                    "counts differ"
                )
            assert isinstance(candidate, dict)
            assert isinstance(evidence, dict)
            try:
                validate_repeat_scope(evidence)
            except ValueError as error:
                raise PostQcError(str(error)) from error
            assert isinstance(audit, dict)
            candidate_id = candidate["candidate_id"]
            if (
                candidate_id != evidence["candidate_id"]
                or candidate_id != audit["candidate_id"]
            ):
                raise PostQcError(
                    "candidate/evidence/audit order differs at logical row "
                    f"{number}: {candidate_id}, {evidence['candidate_id']}, "
                    f"{audit['candidate_id']}"
                )
            key = (
                candidate["seqid"], int(candidate["start"]),
                int(candidate["end"]), candidate["strand"], candidate_id,
            )
            if previous_key is not None and key < previous_key:
                raise PostQcError(
                    "candidate tables are not in canonical coordinate order"
                )
            if not candidate_id or not unique_ids.add(candidate_id):
                raise PostQcError(
                    f"duplicate/empty candidate ID: {candidate_id!r}"
                )
            previous_key = key
            selected = candidate_id in selected_candidate_ids
            eligible = _eligible_bridge(candidate, evidence, audit, bridge_policy)
            if selected or eligible:
                evidence_by_candidate[candidate_id] = evidence
                audit_by_candidate[candidate_id] = audit
            if eligible:
                retained_candidates.append(candidate)
    missing_selected = selected_candidate_ids - set(evidence_by_candidate)
    if missing_selected:
        raise PostQcError(
            "selected candidates are absent from lockstep evidence inputs: "
            f"{sorted(missing_selected)[:10]}"
        )
    return retained_candidates, evidence_by_candidate, audit_by_candidate


def recover_supported_bridges(
    nodes: dict[str, Node],
    candidate_rows: Sequence[Mapping[str, str]],
    evidence_by_candidate: Mapping[str, Mapping[str, str]],
    audit_by_candidate: Mapping[str, Mapping[str, str]],
    minimum_same_frame_bp: int = 3,
    bridge_policy: Mapping[str, object] | None = None,
) -> list[dict[str, object]]:
    if bridge_policy is None:
        bridge_policy = {
            "minimum_rna_tier": 2,
            "minimum_protein_tier": 2,
            "maximum_repeat_cds_fraction": 0.8,
            "exclude_active_fusion": True,
        }
    initial_edges = build_edges(nodes)
    _components, roots = coding_components(nodes, initial_edges, minimum_same_frame_bp)
    initial_transcript_ids = frozenset(
        transcript_id
        for transcript_id in roots
        if nodes[transcript_id].source_status == "RELEASED"
    )
    released_index = IntervalIndex(
        (nodes[transcript_id] for transcript_id in initial_transcript_ids),
        seqid=lambda node: node.model.seqid,
        start=lambda node: node.start,
        end=lambda node: node.end,
        strand=lambda node: node.model.strand,
    )
    existing_signatures = {_signature(node.model) for node in nodes.values()}
    promotions: list[dict[str, object]] = []
    for candidate in candidate_rows:
        candidate_id = candidate["candidate_id"]
        evidence = evidence_by_candidate.get(candidate_id)
        audit = audit_by_candidate.get(candidate_id)
        if evidence is None or audit is None or not _eligible_bridge(
            candidate, evidence, audit, bridge_policy
        ):
            continue
        model = _candidate_model(candidate)
        if _signature(model) in existing_signatures:
            continue
        probe = Node(
            model=model, input_gene_id=".", source_status="PROMOTED_SUPPORTED_BRIDGE",
            selection_locus_id=audit["locus_id"], routes="POSTQC_RNA_OR_PROTEIN_BRIDGE",
            external_support=True, rna_tier=int(evidence["rna_tier"] or 0),
            protein_tier=int(evidence["protein_tier"] or 0),
        )
        # A recovery candidate must bridge two components from the annotation
        # that entered this stage. Earlier promotions are possible isoforms,
        # not new evidence anchors for later candidates.
        hits, hit_roots = _linked_component_hits(
            probe, released_index, roots, minimum_same_frame_bp
        )
        if len(hit_roots) < 2:
            continue
        nodes[model.transcript_id] = probe
        existing_signatures.add(_signature(model))
        promotions.append(
            {
                "candidate_id": candidate_id,
                "transcript_id": model.transcript_id,
                "seqid": model.seqid, "strand": model.strand,
                "start": model.start, "end": model.end,
                "initial_component_count": len(hit_roots),
                "hit_transcript_count": len(hits),
                "hit_transcripts": ",".join(sorted(hits)),
                "rna_tier": probe.rna_tier, "protein_tier": probe.protein_tier,
                "status": "PROMOTED_SUPPORTED_BRIDGE",
            }
        )
    return promotions


def _linked_component_hits(
    probe: Node,
    released_index: IntervalIndex[Node],
    roots: Mapping[str, str],
    minimum_same_frame_bp: int,
) -> tuple[tuple[str, ...], frozenset[str]]:
    """Return released transcripts and component roots truly linked to a probe."""

    hits: list[str] = []
    for node in released_index.query(
        probe.model.seqid, probe.start, probe.end, probe.model.strand
    ):
        transcript_id = node.transcript_id
        if node.retired:
            continue
        if assess_edge(probe, node).linked(minimum_same_frame_bp):
            hits.append(transcript_id)
    ordered_hits = tuple(sorted(hits))
    return ordered_hits, frozenset(roots[transcript_id] for transcript_id in ordered_hits)


def _phased_parts(
    genomic_intervals: Sequence[tuple[int, int]], strand: str
) -> tuple[CdsPart, ...]:
    ordered = sorted(genomic_intervals, reverse=strand == "-")
    cumulative = 0
    parts: list[CdsPart] = []
    for start, end in ordered:
        phase = (3 - cumulative % 3) % 3
        parts.append(CdsPart(start, end, phase))
        cumulative += end - start + 1
    return tuple(sorted(parts))


def recover_strong_rna_bridges(
    nodes: dict[str, Node],
    junctions: Mapping[tuple[str, int, int, str], JunctionAggregate],
    genome: Mapping[str, str],
    minimum_same_frame_bp: int = 3,
    minimum_hq_reads: int = 3,
    minimum_replicated_conditions: int = 1,
    maximum_intron_bp: int = 5000,
    review_maximum_intron_bp: int = 20000,
    review_rows: list[dict[str, object]] | None = None,
    classify_generated=None,
) -> list[dict[str, object]]:
    """Reconstruct complete-ORF bridges at exact, replicated splice boundaries."""

    current_edges = build_edges(nodes)
    _components, roots = coding_components(nodes, current_edges, minimum_same_frame_bp)
    exon_starts: defaultdict[tuple[str, str, int], set[str]] = defaultdict(set)
    exon_ends: defaultdict[tuple[str, str, int], set[str]] = defaultdict(set)
    cds_starts: defaultdict[tuple[str, str, int], set[str]] = defaultdict(set)
    cds_ends: defaultdict[tuple[str, str, int], set[str]] = defaultdict(set)
    active_fusion_nodes = tuple(
        node
        for node in nodes.values()
        if (
            not node.retired
            and node.source_status == "RELEASED"
            and node.active_fusion
        )
    )
    active_fusion_index = IntervalIndex(
        active_fusion_nodes,
        seqid=lambda node: node.model.seqid,
        start=lambda node: node.start,
        end=lambda node: node.end,
        strand=lambda node: node.model.strand,
    )
    for transcript_id, node in nodes.items():
        # Synthetic/promoted transcripts are products of this QC stage, not
        # independent evidence anchors for another synthetic bridge. An active
        # fusion is a hard structural boundary, not a bridge anchor.
        if (
            node.retired
            or node.source_status != "RELEASED"
            or node.active_fusion
        ):
            continue
        model = node.model
        for feature in model.exons:
            exon_starts[(model.seqid, model.strand, feature.start)].add(transcript_id)
            exon_ends[(model.seqid, model.strand, feature.end)].add(transcript_id)
        for feature in model.cds:
            cds_starts[(model.seqid, model.strand, feature.start)].add(transcript_id)
            cds_ends[(model.seqid, model.strand, feature.end)].add(transcript_id)

    existing_signatures = {_signature(node.model) for node in nodes.values()}
    if maximum_intron_bp < 1 or review_maximum_intron_bp < maximum_intron_bp:
        raise PostQcError("RNA bridge automatic/review intron thresholds are invalid")
    initial_transcript_ids = frozenset(
        transcript_id
        for transcript_id in roots
        if nodes[transcript_id].source_status == "RELEASED"
    )
    released_index = IntervalIndex(
        (nodes[transcript_id] for transcript_id in initial_transcript_ids),
        seqid=lambda node: node.model.seqid,
        start=lambda node: node.start,
        end=lambda node: node.end,
        strand=lambda node: node.model.strand,
    )
    proposals: dict[
        tuple[object, ...],
        tuple[
            Node, str, str, int, int, JunctionAggregate,
            tuple[str, ...], frozenset[str],
        ],
    ] = {}
    for (seqid, intron_start, intron_end, strand), aggregate in sorted(junctions.items()):
        if (
            aggregate.total_hq_reads < minimum_hq_reads
            or aggregate.condition_count_ge2rep_hq_ge1 < minimum_replicated_conditions
        ):
            continue
        left_end, right_start = intron_start - 1, intron_end + 1
        left_ids = exon_ends[(seqid, strand, left_end)] & cds_ends[(seqid, strand, left_end)]
        right_ids = exon_starts[(seqid, strand, right_start)] & cds_starts[(seqid, strand, right_start)]
        for left_id in sorted(left_ids):
            for right_id in sorted(right_ids):
                if roots[left_id] == roots[right_id]:
                    continue
                left_parts = sorted(
                    (feature.start, feature.end)
                    for feature in nodes[left_id].model.cds
                    if feature.end <= left_end
                )
                right_parts = sorted(
                    (feature.start, feature.end)
                    for feature in nodes[right_id].model.cds
                    if feature.start >= right_start
                )
                intervals = left_parts + right_parts
                if not left_parts or not right_parts:
                    continue
                ordered_genomic = sorted(intervals)
                if any(right[0] <= left[1] for left, right in zip(ordered_genomic, ordered_genomic[1:])):
                    continue
                parts = _phased_parts(ordered_genomic, strand)
                orf = audit_orf(seqid, strand, parts, genome)
                if orf.status != "complete":
                    continue
                payload = json.dumps(
                    {
                        "seqid": seqid, "strand": strand,
                        "parts": [(part.start, part.end, part.phase) for part in parts],
                    },
                    sort_keys=True, separators=(",", ":"),
                )
                digest = hashlib.sha256(payload.encode()).hexdigest()[:24]
                candidate_id = f"GRNA-{digest}"
                transcript_id = f"GIRNA-{digest}"
                model = TranscriptModel(
                    gene_id=".", transcript_id=transcript_id, candidate_id=candidate_id,
                    seqid=seqid, strand=strand, source="GeneInIt-RNA-bridge",
                    exons=[
                        Feature(seqid, "GeneInIt-RNA-bridge", "exon", part.start, part.end, strand, ".", {})
                        for part in parts
                    ],
                    cds=[
                        Feature(seqid, "GeneInIt-RNA-bridge", "CDS", part.start, part.end, strand, str(part.phase), {})
                        for part in parts
                    ],
                )
                # Preserve active-fusion boundaries even when both exact RNA
                # anchors are non-fusion models on opposite sides. A synthetic
                # bridge may not overlap or span any active-fusion model.
                if next(
                    active_fusion_index.query(
                        model.seqid, model.start, model.end, model.strand
                    ),
                    None,
                ) is not None:
                    continue
                probe = Node(
                    model=model, input_gene_id=".",
                    source_status="PROMOTED_STRONG_RNA_BRIDGE",
                    selection_locus_id=".",
                    routes="POSTQC_STRONG_RNA_BRIDGE_COMPLETE_ORF",
                    external_support=True, rna_tier=3, protein_tier=0,
                )
                hits, hit_roots = _linked_component_hits(
                    probe, released_index, roots, minimum_same_frame_bp
                )
                if len(hit_roots) < 2:
                    continue
                signature = _signature(model)
                if signature in existing_signatures:
                    continue
                proposal = (
                    probe, left_id, right_id, intron_start, intron_end,
                    aggregate, hits, hit_roots,
                )
                previous = proposals.get(signature)
                if previous is None or (left_id, right_id, transcript_id) < (
                    previous[1], previous[2], previous[0].transcript_id
                ):
                    proposals[signature] = proposal

    rows: list[dict[str, object]] = []
    for signature, (
        probe, left_id, right_id, intron_start, intron_end, aggregate,
        hits, hit_roots,
    ) in sorted(
        proposals.items(), key=lambda item: item[1][0].transcript_id
    ):
        model = probe.model
        intron_bp = intron_end - intron_start + 1
        if intron_bp > maximum_intron_bp:
            if review_rows is not None:
                within_review = intron_bp <= review_maximum_intron_bp
                review_rows.append(
                    {
                        "candidate_id": model.candidate_id,
                        "transcript_id": model.transcript_id,
                        "seqid": model.seqid, "strand": model.strand,
                        "start": model.start, "end": model.end,
                        "intron_start": intron_start, "intron_end": intron_end,
                        "intron_bp": intron_bp,
                        "left_transcript_id": left_id,
                        "right_transcript_id": right_id,
                        "total_hq_reads": aggregate.total_hq_reads,
                        "sample_count_hq_ge1": aggregate.sample_count_hq_ge1,
                        "condition_count_ge2rep_hq_ge1": (
                            aggregate.condition_count_ge2rep_hq_ge1
                        ),
                        "condition_count_ge2rep_hq_ge3": (
                            aggregate.condition_count_ge2rep_hq_ge3
                        ),
                        "automatic_max_intron_bp": maximum_intron_bp,
                        "review_max_intron_bp": review_maximum_intron_bp,
                        "threshold_class": (
                            "ADDITIONAL_EVIDENCE" if within_review else "EXTREME"
                        ),
                        "status": (
                            "HELD_LONG_INTRON_ADDITIONAL_EVIDENCE_REQUIRED"
                            if within_review
                            else "KEPT_SEPARATE_EXTREME_INTRON"
                        ),
                        "required_evidence": (
                            "SPANNING_PROTEIN_OR_FULL_LENGTH_TRANSCRIPT"
                            if within_review
                            else "FULL_LENGTH_TRANSCRIPT_OR_PROTEIN_PLUS_RNA"
                        ),
                    }
                )
            continue
        if classify_generated is not None and not classify_generated(
            model, nodes[left_id].candidate_id, nodes[right_id].candidate_id
        ):
            continue
        nodes[model.transcript_id] = probe
        existing_signatures.add(signature)
        rows.append(
            {
                "candidate_id": model.candidate_id,
                "transcript_id": model.transcript_id,
                "seqid": model.seqid, "strand": model.strand,
                "start": model.start, "end": model.end,
                "initial_component_count": len(hit_roots),
                "hit_transcript_count": len(hits),
                "hit_transcripts": ",".join(hits),
                "rna_tier": 3, "protein_tier": 0,
                "status": (
                    "PROMOTED_STRONG_RNA_BRIDGE_COMPLETE_ORF:"
                    f"{intron_start}-{intron_end}"
                ),
            }
        )
    return rows


def _adjacency(
    nodes: Mapping[str, Node], edges: Sequence[CodingEdge], minimum_same_frame_bp: int = 3
) -> dict[str, set[str]]:
    adjacency = {key: set() for key, node in nodes.items() if not node.retired}
    for edge in edges:
        if edge.linked(minimum_same_frame_bp):
            adjacency[edge.left].add(edge.right)
            adjacency[edge.right].add(edge.left)
    return adjacency


def _connector_candidates(
    adjacency: Mapping[str, set[str]],
    nodes: Mapping[str, Node],
    minimum_protected_components: int,
) -> dict[str, tuple[int, int, str]]:
    """Find support-separating articulation nodes in O(V+E) time."""

    discovery: dict[str, int] = {}
    low: dict[str, int] = {}
    parent: dict[str, str | None] = {}
    children: defaultdict[str, list[str]] = defaultdict(list)
    subtree_size: dict[str, int] = {}
    subtree_support: dict[str, int] = {}
    result: dict[str, tuple[int, int, str]] = {}
    clock = 0

    for root in sorted(adjacency):
        if root in discovery or nodes[root].retired:
            continue
        component: list[str] = []
        parent[root] = None
        clock += 1
        discovery[root] = low[root] = clock
        stack: list[tuple[str, Iterable[str]]] = [
            (root, iter(sorted(adjacency[root])))
        ]
        while stack:
            current, neighbors = stack[-1]
            try:
                neighbor = next(neighbors)
            except StopIteration:
                stack.pop()
                component.append(current)
                subtree_size[current] = 1 + sum(
                    subtree_size[child] for child in children[current]
                )
                subtree_support[current] = int(nodes[current].external_support) + sum(
                    subtree_support[child] for child in children[current]
                )
                parent_id = parent[current]
                if parent_id is not None:
                    low[parent_id] = min(low[parent_id], low[current])
                continue
            if nodes[neighbor].retired:
                continue
            if neighbor not in discovery:
                parent[neighbor] = current
                children[current].append(neighbor)
                clock += 1
                discovery[neighbor] = low[neighbor] = clock
                stack.append((neighbor, iter(sorted(adjacency[neighbor]))))
            elif neighbor != parent[current]:
                low[current] = min(low[current], discovery[neighbor])

        component_size = len(component)
        component_support = sum(nodes[value].external_support for value in component)
        component_key = min(component)
        for transcript_id in component:
            node = nodes[transcript_id]
            if (
                node.retired
                or node.external_support
                or node.active_fusion
                or len(adjacency[transcript_id]) < 2
            ):
                continue
            cut_children = [
                child for child in children[transcript_id]
                if parent[transcript_id] is None
                or low[child] >= discovery[transcript_id]
            ]
            piece_sizes = [subtree_size[child] for child in cut_children]
            piece_supports = [
                subtree_support[child] for child in cut_children
            ]
            if parent[transcript_id] is not None:
                residual_size = component_size - 1 - sum(piece_sizes)
                if residual_size:
                    piece_sizes.append(residual_size)
                    piece_supports.append(
                        component_support
                        - int(node.external_support)
                        - sum(subtree_support[child] for child in cut_children)
                    )
            neighbor_components = len(piece_sizes)
            protected_components = sum(value > 0 for value in piece_supports)
            if protected_components >= minimum_protected_components:
                result[transcript_id] = (
                    neighbor_components, protected_components, component_key
                )
    return result


def retire_unsupported_connectors(
    nodes: dict[str, Node], minimum_same_frame_bp: int = 3,
    minimum_protected_components: int = 2,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    edges = build_edges(nodes)
    adjacency = _adjacency(nodes, edges, minimum_same_frame_bp)
    candidates = _connector_candidates(
        adjacency, nodes, minimum_protected_components
    )
    chosen_by_component: dict[str, str] = {}
    for transcript_id in sorted(candidates):
        component_key = candidates[transcript_id][2]
        chosen_by_component.setdefault(component_key, transcript_id)
    for transcript_id in sorted(chosen_by_component.values()):
        neighbor_count, protected_count, _component_key = candidates[transcript_id]
        node = nodes[transcript_id]
        node.retired = True
        node.operation = "RETIRE_UNSUPPORTED_FUSION_CONNECTOR"
        node.reason = "WEAK_ARTICULATION_BETWEEN_SUPPORTED_CODING_COMPONENTS"
        rows.append(
            {
                "transcript_id": transcript_id, "candidate_id": node.candidate_id,
                "input_gene_id": node.input_gene_id, "seqid": node.model.seqid,
                "strand": node.model.strand,
                "neighbor_component_count": neighbor_count,
                "protected_neighbor_component_count": protected_count,
                "status": "RETIRED_UNSUPPORTED_FUSION_CONNECTOR",
                "reason": node.reason,
            }
        )
    return rows


def _node_state_sha256(nodes: Mapping[str, Node]) -> str:
    """Hash only the mutable structural state used by convergence decisions."""

    payload = [
        {
            "transcript_id": transcript_id,
            "candidate_id": node.candidate_id,
            "source_status": node.source_status,
            "retired": node.retired,
            "seqid": node.model.seqid,
            "strand": node.model.strand,
            "cds": [
                (feature.start, feature.end, feature.phase)
                for feature in sorted(
                    node.model.cds,
                    key=lambda item: (item.start, item.end, item.phase),
                )
            ],
        }
        for transcript_id, node in sorted(nodes.items())
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _deduplicate_rows(
    rows: Iterable[Mapping[str, object]],
) -> list[dict[str, object]]:
    observed: set[str] = set()
    unique: list[dict[str, object]] = []
    for row in rows:
        key = json.dumps(dict(row), sort_keys=True, separators=(",", ":"))
        if key in observed:
            continue
        observed.add(key)
        unique.append(dict(row))
    return unique


def iterate_structural_qc(
    nodes: dict[str, Node],
    maximum_iterations: int,
    apply_iteration: Callable[
        [],
        tuple[
            list[dict[str, object]],
            list[dict[str, object]],
            list[dict[str, object]],
            list[dict[str, object]],
        ],
    ],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    int,
]:
    """Apply structural correction until a no-change iteration proves convergence."""

    if maximum_iterations < 1:
        raise PostQcError("post-QC maximum_iterations must be positive")
    promotions: list[dict[str, object]] = []
    bridge_reviews: list[dict[str, object]] = []
    retirements: list[dict[str, object]] = []
    iteration_rows: list[dict[str, object]] = []
    for iteration in range(1, maximum_iterations + 1):
        before = _node_state_sha256(nodes)
        supported, strong_rna, reviews, retired = apply_iteration()
        after = _node_state_sha256(nodes)
        changed = before != after
        promotions.extend(supported)
        promotions.extend(strong_rna)
        bridge_reviews.extend(reviews)
        retirements.extend(retired)
        iteration_rows.append(
            {
                "iteration": iteration,
                "state_before_sha256": before,
                "state_after_sha256": after,
                "supported_bridge_promotions": len(supported),
                "strong_rna_bridge_promotions": len(strong_rna),
                "bridge_review_candidates": len(reviews),
                "connector_retirements": len(retired),
                "changed": int(changed),
                "converged": int(not changed),
            }
        )
        if not changed:
            return (
                _deduplicate_rows(promotions),
                _deduplicate_rows(bridge_reviews),
                _deduplicate_rows(retirements),
                iteration_rows,
                iteration,
            )
    raise PostQcError(
        "post-QC did not converge within "
        f"maximum_iterations={maximum_iterations}; last state still changed"
    )


def _stable_gene_id(nodes: Iterable[Node]) -> str:
    payload = json.dumps(sorted(node.candidate_id for node in nodes), separators=(",", ":"))
    return f"GPQC-{hashlib.sha256(payload.encode()).hexdigest()[:24]}"


def _cds_bp(model: TranscriptModel) -> int:
    return sum(feature.end - feature.start + 1 for feature in model.cds)


def _canonical_key(node: Node) -> tuple[object, ...]:
    # Canonical rank is intentionally based on coding evidence. UTRs are
    # reconstructed only after post-QC and therefore cannot participate here.
    return (
        -int(node.external_support),
        -int(node.source_status == "RELEASED"),
        -max(node.rna_tier, node.protein_tier),
        -_cds_bp(node.model),
        node.candidate_id,
    )


def edge_audit_rows(
    nodes: Mapping[str, Node],
    initial_edges: Sequence[CodingEdge],
    final_edges: Sequence[CodingEdge],
    minimum_same_frame_bp: int,
) -> list[dict[str, object]]:
    """Audit the union of input and final edges, including promoted nodes."""

    initial_by_pair = {(edge.left, edge.right): edge for edge in initial_edges}
    final_by_pair = {(edge.left, edge.right): edge for edge in final_edges}
    rows: list[dict[str, object]] = []
    for pair in sorted(set(initial_by_pair) | set(final_by_pair)):
        edge = (
            final_by_pair[pair] if pair in final_by_pair else initial_by_pair[pair]
        )
        left, right = nodes[edge.left], nodes[edge.right]
        present_initial = pair in initial_by_pair
        present_final = pair in final_by_pair
        if edge.active_fusion_boundary:
            disposition = "KEEP_SEPARATE_ACTIVE_FUSION"
        elif edge.linked(minimum_same_frame_bp):
            disposition = (
                "RETIRED_CONNECTOR" if left.retired or right.retired
                else "SAME_POSTQC_GENE"
                if left.postqc_gene_id == right.postqc_gene_id
                else "ERROR_UNRESOLVED_CODING_LINK"
            )
        else:
            disposition = "KEEP_SEPARATE_FRAME_OR_STRUCTURE"
        rows.append(
            {
                "left_transcript_id": edge.left,
                "right_transcript_id": edge.right,
                "left_candidate_id": left.candidate_id,
                "right_candidate_id": right.candidate_id,
                "left_input_gene_id": left.input_gene_id,
                "right_input_gene_id": right.input_gene_id,
                "seqid": left.model.seqid,
                "strand": left.model.strand,
                "span_overlap_bp": edge.span_overlap_bp,
                "cds_overlap_bp": edge.cds_overlap_bp,
                "same_frame_cds_overlap_bp": edge.same_frame_cds_overlap_bp,
                "shared_cds_intron_count": edge.shared_cds_intron_count,
                "coding_link": int(edge.linked(minimum_same_frame_bp)),
                "active_fusion_boundary": int(edge.active_fusion_boundary),
                "present_initial": int(present_initial),
                "present_final": int(present_final),
                "introduced_by_promotion": int(
                    not present_initial
                    and (
                        left.source_status.startswith("PROMOTED_")
                        or right.source_status.startswith("PROMOTED_")
                    )
                ),
                "crosses_input_gene_boundary": int(
                    left.input_gene_id != right.input_gene_id
                ),
                "disposition": disposition,
            }
        )
    return rows


def _component_operation(
    component_nodes: Sequence[Node], old_gene_component_counts: Mapping[str, int]
) -> str:
    old_genes = {node.input_gene_id for node in component_nodes if node.input_gene_id not in EMPTY}
    merged = len(old_genes) > 1
    split = any(old_gene_component_counts[gene_id] > 1 for gene_id in old_genes)
    promoted = any(node.source_status.startswith("PROMOTED_") for node in component_nodes)
    if merged and split:
        return "MERGE_AND_SPLIT_TRANSCRIPT_COMPONENTS"
    if merged:
        return "CONSOLIDATE_GENE_IDS_AS_ISOFORMS"
    if split:
        return "SPLIT_DISCONNECTED_CODING_COMPONENTS"
    if promoted:
        return "ADD_SUPPORTED_BRIDGE_TRANSCRIPT"
    return "UNCHANGED"


def _write_gff_row(
    handle: TextIO, model: TranscriptModel, feature_type: str, start: int, end: int,
    phase: str, attributes: Iterable[tuple[str, str]], source: str,
) -> None:
    handle.write(
        f"{model.seqid}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{model.strand}\t{phase}\t"
        f"{format_gff3_attributes(attributes)}\n"
    )


def _write_corrected_gff3(
    path: Path,
    ordered_components: Sequence[Sequence[str]],
    nodes: Mapping[str, Node],
    operations: Mapping[str, str],
) -> None:
    with atomic_text(path) as handle:
        handle.write("##gff-version 3\n")
        for transcript_ids in ordered_components:
            component_nodes = [nodes[value] for value in transcript_ids]
            gene_id = component_nodes[0].postqc_gene_id
            ranked = sorted(component_nodes, key=_canonical_key)
            gene_start = min(node.start for node in component_nodes)
            gene_end = max(node.end for node in component_nodes)
            _write_gff_row(
                handle, ranked[0].model, "gene", gene_start, gene_end, ".",
                (("ID", gene_id), ("transcript_count", str(len(ranked))),
                 ("canonical_transcript", ranked[0].transcript_id),
                 ("postqc_operation", operations[gene_id])),
                "GeneInIt-postQC",
            )
            for rank, node in enumerate(ranked, 1):
                model = node.model
                attrs = [
                    ("ID", model.transcript_id), ("Parent", gene_id),
                    ("candidate_id", model.candidate_id),
                    ("canonical", "true" if rank == 1 else "false"),
                    ("source_status", node.source_status),
                ]
                if model.utr_refine_evidence not in EMPTY:
                    attrs.append(("utr_refine_evidence", model.utr_refine_evidence))
                _write_gff_row(
                    handle, model, "mRNA", model.start, model.end, ".", attrs,
                    model.source,
                )
                for feature_type, features in (
                    ("exon", model.exons), ("CDS", model.cds),
                    ("five_prime_UTR", model.five_prime_utrs),
                    ("three_prime_UTR", model.three_prime_utrs),
                ):
                    for index, feature in enumerate(sorted(features, key=lambda item: (item.start, item.end, item.phase)), 1):
                        _write_gff_row(
                            handle, model, feature_type, feature.start, feature.end,
                            feature.phase if feature_type == "CDS" else ".",
                            (("ID", f"{feature_type}-{model.transcript_id}-{index:03d}"),
                             ("Parent", model.transcript_id)), feature.source,
                        )


def run_postqc(
    grouped_gff3: Path,
    grouping_mapping: Path,
    candidate_table: Path,
    component_evidence: Path,
    selection_audit: Path,
    junction_path: Path,
    genome_path: Path,
    policy_path: Path,
    reference_sequence_ids: Sequence[str],
    output_gff3: Path,
    transcript_decisions: Path,
    gene_decisions: Path,
    edge_audit: Path,
    promotions_path: Path,
    bridge_reviews_path: Path,
    retirements_path: Path,
    summary_path: Path,
    complete_path: Path,
    iterations_path: Path | None = None,
    te_evidence_path: Path | None = None,
    assembly_id: str | None = None,
    te_classification_path: Path | None = None,
) -> None:
    with policy_path.open(encoding="utf-8") as handle:
        policy_document = json.load(handle)
    postqc_schema = policy_document.get("postqc_schema")
    if postqc_schema not in {"geneinit_postqc_policy_v1", "geneinit_postqc_policy_v2"}:
        raise PostQcError("policy does not contain a supported GeneInIt post-QC schema")
    postqc_policy = policy_document["postqc"]
    maximum_iterations = int(postqc_policy["maximum_iterations"])
    if maximum_iterations < 1:
        raise PostQcError("post-QC maximum_iterations must be positive")
    minimum_same_frame_bp = int(postqc_policy["minimum_same_frame_overlap_bp"])
    if minimum_same_frame_bp < 1:
        raise PostQcError("post-QC minimum same-frame overlap must be positive")
    bridge_policy = postqc_policy["bridge_recovery"]
    rna_bridge_policy = postqc_policy["rna_bridge_reconstruction"]
    if postqc_schema == "geneinit_postqc_policy_v1":
        fixed_automatic = int(rna_bridge_policy["maximum_intron_bp"])
        length_policy: Mapping[str, object] = {
            "mode": "fixed",
            "fixed_automatic_max_intron_bp": fixed_automatic,
            "fixed_review_max_intron_bp": max(fixed_automatic, fixed_automatic * 4),
        }
    else:
        length_policy = rna_bridge_policy["intron_length_policy"]
    minimum_protected_components = int(
        postqc_policy["unsupported_connector"]["minimum_protected_neighbor_components"]
    )
    if minimum_protected_components < 2:
        raise PostQcError("connector retirement requires at least two protected components")
    order, models_by_transcript = load_raw_gff3(grouped_gff3)
    grouping_rows = _read_table(
        grouping_mapping,
        {
            "candidate_id", "temporary_transcript_id", "provisional_gene_id",
            "selection_locus_id", "routes", "active_fusion",
        },
    )
    grouping_by_candidate = {row["candidate_id"]: row for row in grouping_rows}
    if set(grouping_by_candidate) != {model.candidate_id for model in models_by_transcript.values()}:
        raise PostQcError("grouped GFF3 and grouping mapping candidate sets differ")
    selected_candidate_ids = {
        model.candidate_id for model in models_by_transcript.values()
    }
    (
        candidate_rows,
        evidence_by_candidate,
        audit_by_candidate,
    ) = load_postqc_candidate_inputs(
        candidate_table,
        component_evidence,
        selection_audit,
        selected_candidate_ids,
        bridge_policy,
    )

    nodes: dict[str, Node] = {}
    for transcript_id in order:
        model = models_by_transcript[transcript_id]
        mapping = grouping_by_candidate[model.candidate_id]
        evidence = evidence_by_candidate[model.candidate_id]
        selection = audit_by_candidate[model.candidate_id]
        if selection["selected"] != "1":
            raise PostQcError(
                f"grouped transcript was not selected: {model.candidate_id}"
            )
        if selection["active_fusion"] != mapping["active_fusion"]:
            raise PostQcError(
                f"active-fusion lineage mismatch: {model.candidate_id}"
            )
        routes = mapping["routes"]
        rna_tier = int(evidence["rna_tier"] or 0)
        protein_tier = int(evidence["protein_tier"] or 0)
        external = "RNA" in routes or "PROTEIN" in routes or rna_tier >= 2 or protein_tier >= 2
        nodes[transcript_id] = Node(
            model=model,
            input_gene_id=model.gene_id,
            source_status="RELEASED",
            selection_locus_id=mapping["selection_locus_id"],
            routes=routes,
            external_support=external,
            rna_tier=rna_tier,
            protein_tier=protein_tier,
            active_fusion=mapping["active_fusion"] == "1",
        )
    if len(nodes) != len(models_by_transcript):
        raise PostQcError("duplicate transcript IDs in grouped annotation")

    junctions = load_junction_aggregate(junction_path)
    minimum_hq_reads = int(rna_bridge_policy["minimum_hq_reads"])
    minimum_replicated_conditions = int(
        rna_bridge_policy["minimum_replicated_conditions"]
    )
    training_lengths = trusted_intron_lengths(
        nodes, evidence_by_candidate, junctions, minimum_hq_reads,
        minimum_replicated_conditions,
    )
    thresholds = derive_rna_bridge_thresholds(training_lengths, length_policy)

    initial_edges = build_edges(nodes)
    genome = read_fasta(genome_path)
    classifier = TeClassifier(te_evidence_path, assembly_id) if confirmed_only(bridge_policy) else None
    generated_classifications = {}

    def classify_generated(model, left_candidate, right_candidate):
        identity = model_identity(model)
        parts = tuple(CdsPart(f.start, f.end, int(f.phase)) for f in sorted(model.cds, key=lambda f: f.start))
        sequence = audit_orf(model.seqid, model.strand, parts, genome).protein_sequence
        result = classifier.classify(identity, protein_sha256(sequence))
        generated_classifications[model.candidate_id] = {
            **identity, **result, "origin_stage": "POSTQC",
            "allowed_route": "POSTQC_STRONG_RNA_BRIDGE_COMPLETE_ORF",
            "parent_candidate_ids": [left_candidate, right_candidate],
        }
        return result["candidate_pool"] == "GENERAL"

    def apply_iteration() -> tuple[
        list[dict[str, object]],
        list[dict[str, object]],
        list[dict[str, object]],
        list[dict[str, object]],
    ]:
        supported = recover_supported_bridges(
            nodes, candidate_rows, evidence_by_candidate, audit_by_candidate,
            minimum_same_frame_bp, bridge_policy,
        )
        reviews: list[dict[str, object]] = []
        strong_rna = recover_strong_rna_bridges(
            nodes, junctions, genome, minimum_same_frame_bp,
            minimum_hq_reads, minimum_replicated_conditions,
            thresholds.automatic_max_intron_bp,
            thresholds.review_max_intron_bp,
            reviews,
            classify_generated if classifier is not None else None,
        )
        retired = retire_unsupported_connectors(
            nodes, minimum_same_frame_bp, minimum_protected_components
        )
        return supported, strong_rna, reviews, retired

    (
        promotions,
        bridge_reviews,
        retirements,
        iteration_rows,
        converged_iteration,
    ) = iterate_structural_qc(
        nodes, maximum_iterations, apply_iteration
    )
    if te_classification_path is not None:
        with atomic_text(te_classification_path) as handle:
            json.dump({"schema": "geneinit_postqc_te_v2", "candidates": generated_classifications}, handle, indent=2)
            handle.write("\n")
    final_edges = build_edges(nodes)
    final_components, final_roots = coding_components(
        nodes, final_edges, minimum_same_frame_bp
    )

    old_gene_roots: defaultdict[str, set[str]] = defaultdict(set)
    for transcript_id, node in nodes.items():
        if not node.retired and node.input_gene_id not in EMPTY:
            old_gene_roots[node.input_gene_id].add(final_roots[transcript_id])
    old_gene_component_counts = {key: len(values) for key, values in old_gene_roots.items()}
    sequence_rank = {seqid: index for index, seqid in enumerate(reference_sequence_ids)}
    ordered_components = sorted(
        final_components.values(),
        key=lambda values: (
            sequence_rank[nodes[values[0]].model.seqid],
            min(nodes[value].start for value in values),
            max(nodes[value].end for value in values),
            nodes[values[0]].model.strand,
            min(values),
        ),
    )
    operations: dict[str, str] = {}
    gene_rows: list[dict[str, object]] = []
    for transcript_ids in ordered_components:
        component_nodes = [nodes[value] for value in transcript_ids]
        gene_id = _stable_gene_id(component_nodes)
        operation = _component_operation(component_nodes, old_gene_component_counts)
        operations[gene_id] = operation
        ranked = sorted(component_nodes, key=_canonical_key)
        for rank, node in enumerate(ranked, 1):
            node.postqc_gene_id = gene_id
            node.canonical_rank = rank
            if not node.operation:
                node.operation = operation
                node.reason = "TRANSCRIPT_FIRST_CODING_COMPONENT"
        old_genes = sorted({node.input_gene_id for node in component_nodes if node.input_gene_id not in EMPTY})
        gene_rows.append(
            {
                "postqc_gene_id": gene_id,
                "seqid": component_nodes[0].model.seqid,
                "strand": component_nodes[0].model.strand,
                "start": min(node.start for node in component_nodes),
                "end": max(node.end for node in component_nodes),
                "transcript_count": len(component_nodes),
                "input_gene_count": len(old_genes),
                "input_gene_ids": ",".join(old_genes) or ".",
                "promoted_transcript_count": sum(node.source_status.startswith("PROMOTED_") for node in component_nodes),
                "operation": operation,
            }
        )

    _write_corrected_gff3(output_gff3, ordered_components, nodes, operations)
    transcript_rows = [
        {
            "transcript_id": node.transcript_id,
            "candidate_id": node.candidate_id,
            "source_status": node.source_status,
            "input_gene_id": node.input_gene_id,
            "postqc_gene_id": node.postqc_gene_id if not node.retired else ".",
            "canonical_transcript": int(not node.retired and node.canonical_rank == 1),
            "canonical_rank": node.canonical_rank if not node.retired else 0,
            "seqid": node.model.seqid, "strand": node.model.strand,
            "start": node.start, "end": node.end,
            "selection_locus_id": node.selection_locus_id,
            "routes": node.routes,
            "external_support": int(node.external_support),
            "active_fusion": int(node.active_fusion),
            "operation": node.operation,
            "reason": node.reason,
        }
        for node in sorted(nodes.values(), key=lambda item: (item.model.seqid, item.start, item.end, item.transcript_id))
    ]
    _write_table(transcript_decisions, TRANSCRIPT_COLUMNS, transcript_rows)
    _write_table(gene_decisions, GENE_COLUMNS, gene_rows)
    _write_table(promotions_path, PROMOTION_COLUMNS, promotions)
    _write_table(bridge_reviews_path, BRIDGE_REVIEW_COLUMNS, bridge_reviews)
    _write_table(retirements_path, RETIREMENT_COLUMNS, retirements)
    if iterations_path is not None:
        _write_table(iterations_path, ITERATION_COLUMNS, iteration_rows)

    final_edge_by_pair = {(edge.left, edge.right): edge for edge in final_edges}
    edge_rows = edge_audit_rows(
        nodes, initial_edges, final_edges, minimum_same_frame_bp
    )
    _write_table(edge_audit, EDGE_COLUMNS, edge_rows)

    # Fixed-point validation: each coding component is exactly one gene and
    # each gene is exactly one coding component after correction.
    unresolved_links = 0
    for edge in final_edge_by_pair.values():
        if edge.linked(minimum_same_frame_bp) and nodes[edge.left].postqc_gene_id != nodes[edge.right].postqc_gene_id:
            unresolved_links += 1
    roots_per_gene: defaultdict[str, set[str]] = defaultdict(set)
    genes_per_root: defaultdict[str, set[str]] = defaultdict(set)
    for transcript_id, root in final_roots.items():
        gene_id = nodes[transcript_id].postqc_gene_id
        roots_per_gene[gene_id].add(root)
        genes_per_root[root].add(gene_id)
    split_genes = sum(len(values) != 1 for values in roots_per_gene.values())
    merged_roots = sum(len(values) != 1 for values in genes_per_root.values())
    if unresolved_links or split_genes or merged_roots:
        raise PostQcError(
            f"post-QC fixed point failed: links={unresolved_links}, split_genes={split_genes}, merged_roots={merged_roots}"
        )
    operation_counts = Counter(row["operation"] for row in gene_rows)
    summary_rows = [
        ("status", "PASS_MIGRATION_PENDING"),
        ("input_gene_count", len({node.input_gene_id for node in nodes.values() if node.input_gene_id not in EMPTY and node.source_status == "RELEASED"})),
        ("input_transcript_count", sum(node.source_status == "RELEASED" for node in nodes.values())),
        ("promoted_transcript_count", len(promotions)),
        ("retired_transcript_count", len(retirements)),
        ("output_gene_count", len(gene_rows)),
        ("output_transcript_count", sum(not node.retired for node in nodes.values())),
        ("initial_overlap_pair_count", len(initial_edges)),
        ("minimum_same_frame_overlap_bp", minimum_same_frame_bp),
        ("maximum_iterations", maximum_iterations),
        ("iterations_executed", len(iteration_rows)),
        ("converged_iteration", converged_iteration),
        ("rna_bridge_threshold_mode_requested", thresholds.requested_mode),
        ("rna_bridge_threshold_mode_effective", thresholds.effective_mode),
        ("rna_bridge_training_intron_count", thresholds.training_intron_count),
        ("rna_bridge_automatic_quantile", thresholds.automatic_quantile),
        ("rna_bridge_automatic_quantile_bp", thresholds.automatic_quantile_bp),
        ("rna_bridge_review_quantile", thresholds.review_quantile),
        ("rna_bridge_review_quantile_bp", thresholds.review_quantile_bp),
        ("maximum_synthetic_rna_bridge_intron_bp", thresholds.automatic_max_intron_bp),
        ("rna_bridge_review_maximum_intron_bp", thresholds.review_max_intron_bp),
        ("rna_bridge_review_candidate_count", len(bridge_reviews)),
        ("rna_bridge_additional_evidence_candidate_count", sum(row["threshold_class"] == "ADDITIONAL_EVIDENCE" for row in bridge_reviews)),
        ("rna_bridge_extreme_candidate_count", sum(row["threshold_class"] == "EXTREME" for row in bridge_reviews)),
        ("final_coding_link_count", sum(edge.linked(minimum_same_frame_bp) for edge in final_edges)),
        ("unresolved_coding_link_count", unresolved_links),
        ("multi_component_gene_count", split_genes),
        ("multi_gene_component_count", merged_roots),
    ]
    summary_rows.extend((f"operation::{key}", value) for key, value in sorted(operation_counts.items()))
    _write_table(summary_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in summary_rows))
    completion: list[tuple[str, object]] = [
        ("schema_version", 1), ("status", "PASS_MIGRATION_PENDING"),
        ("fixed_point_status", "PASS"),
        ("maximum_iterations", maximum_iterations),
        ("iterations_executed", len(iteration_rows)),
        ("converged_iteration", converged_iteration),
        ("policy_id", policy_document["policy_id"]),
        ("policy_sha256", sha256_file(policy_path)),
        ("rna_bridge_threshold_mode_effective", thresholds.effective_mode),
        ("rna_bridge_training_intron_count", thresholds.training_intron_count),
        ("maximum_synthetic_rna_bridge_intron_bp", thresholds.automatic_max_intron_bp),
        ("rna_bridge_review_maximum_intron_bp", thresholds.review_max_intron_bp),
        ("rna_bridge_review_candidate_count", len(bridge_reviews)),
        ("output_gene_count", len(gene_rows)),
        ("output_transcript_count", sum(not node.retired for node in nodes.values())),
        ("promoted_transcript_count", len(promotions)),
        ("retired_transcript_count", len(retirements)),
        ("corrected_gff3_sha256", sha256_file(output_gff3)),
        ("transcript_decisions_sha256", sha256_file(transcript_decisions)),
        ("gene_decisions_sha256", sha256_file(gene_decisions)),
        ("edge_audit_sha256", sha256_file(edge_audit)),
        ("promotions_sha256", sha256_file(promotions_path)),
        ("bridge_reviews_sha256", sha256_file(bridge_reviews_path)),
        ("retirements_sha256", sha256_file(retirements_path)),
        ("summary_sha256", sha256_file(summary_path)),
    ]
    if te_classification_path is not None:
        completion.append(("te_classification_sha256", sha256_file(te_classification_path)))
    if iterations_path is not None:
        completion.append(("iterations_sha256", sha256_file(iterations_path)))
    _write_table(complete_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in completion))
