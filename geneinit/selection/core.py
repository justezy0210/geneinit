"""Clean implementation of the comparison-free v2.2H evidence core.

This module deliberately keeps component metrics independent.  It applies
explicit boolean routes, repeat gates, and a deterministic conflict graph; it
never computes a weighted evidence score.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence

from geneinit.contracts.streaming import DiskBackedUniqueIds
from geneinit.evidence.metrics import interval_bp, overlap_bp
from geneinit.evidence.rna_cds import (
    RNA_CDS_BASE_COLUMNS, RNA_CDS_COLUMNS, RNA_CDS_SAMPLE_COLUMNS, RNA_CDS_STATUSES,
)
from geneinit.evidence.repeats import validate_repeat_scope
from geneinit.evidence.te_classification import TE_COLUMNS, confirmed_only, repeat_neutral, partition_pools
from geneinit.models.candidate import CdsPart


AUDIT_COLUMNS = (
    "candidate_id", "structure_id", "seqid", "strand", "start", "end",
    "locus_id", "policy_id", "hard_qc_pass", "rna_calibration_mode",
    "strong_recurrent_rna", "recurrent_rna", "strong_unique_protein",
    "strong_unique_full_protein", "moderate_unique_protein", "predictor_exact",
    "repeat_scope", "repeat_cds_fraction", "all_repeat_cds_fraction",
    "strong_intron_count", "intron_count", "repeat_override_reason",
    "rna_recurrent_intron_count", "rna_recurrent_intron_fraction",
    "rna_nonrecurrent_introns_json",
    "repeat_risk", "repeat_override", "balanced_base_selected",
    "base_supplement", "eligible", "selected", "active_fusion",
    "rna_cds_conflict_penalty", "rna_cds_complete_support", "selection_rank_json",
    "fusion_center_precedence",
    "routes", "decision", "blockers", "reasons",
)

AUDIT_COLUMNS += TE_COLUMNS + RNA_CDS_COLUMNS + RNA_CDS_SAMPLE_COLUMNS

EDGE_COLUMNS = (
    "locus_id", "seqid", "strand", "left_candidate_id", "right_candidate_id",
    "cds_overlap_bp", "shorter_cds_overlap_fraction", "shared_intron_count",
    "exact_intron_chain", "hard_overlap", "relation",
)

COMPONENT_REQUIRED_COLUMNS = {
    "candidate_id", "structure_id", "seqid", "strand", "start", "end",
    "cds_signature", "introns_json", "intron_count", "orf_status",
    "common_qc_pass", "translated_protein_aa", "rna_tier",
    "rna_hq_read_counts", "rna_sample_hq_ge1_counts",
    "rna_condition_ge2rep_hq_ge1_counts", "protein_tier",
    "protein_multiplicity_status", "protein_target_linked",
    "predictor_tier", "predictors", "repeat_scope", "repeat_cds_fraction", "splice_motifs",
    "noncanonical_splice_count", "ambiguous_splice_count",
}


class SelectionError(RuntimeError):
    """Raised when the selection graph or evidence table is inconsistent."""


def _integer(row: Mapping[str, str], key: str) -> int:
    try:
        return int(row[key])
    except (KeyError, ValueError) as error:
        raise SelectionError(f"{key} must be an integer") from error


def _numeric(row: Mapping[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "NA")
    if value in {"", ".", "NA"}:
        return default
    try:
        return float(value)
    except ValueError as error:
        raise SelectionError(f"{key} must be numeric or NA") from error


def _json_ints(row: Mapping[str, str], key: str) -> tuple[int, ...]:
    try:
        values = json.loads(row[key])
        result = tuple(int(value) for value in values)
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise SelectionError(f"{key} must be a JSON integer array") from error
    return result


def _parse_cds(text: str) -> tuple[CdsPart, ...]:
    parts: list[CdsPart] = []
    for token in text.split(","):
        try:
            coordinates, phase = token.rsplit(":", 1)
            start, end = coordinates.split("-", 1)
            parts.append(CdsPart(int(start), int(end), int(phase)))
        except (ValueError, TypeError) as error:
            raise SelectionError(f"invalid CDS signature: {text!r}") from error
    if not parts:
        raise SelectionError("empty CDS signature")
    return tuple(sorted(parts))


@dataclass(slots=True)
class Record:
    row: dict[str, str]
    cds: tuple[CdsPart, ...]
    locus_id: str = ""
    hard_qc: bool = False
    rna_mode: str = ""
    strong_rna: bool = False
    recurrent_rna: bool = False
    strong_unique_protein: bool = False
    strong_unique_full_protein: bool = False
    moderate_unique_protein: bool = False
    predictor_exact: bool = False
    repeat_risk: str = "LOW"
    repeat_override: bool = False
    repeat_override_reason: str = "NOT_EVALUATED"
    routes: list[str] = field(default_factory=list)
    eligible: bool = False
    selected: bool = False
    rna_cds_ranking: bool = False
    rna_cds_competition: bool = False
    rna_cds_conflict: bool = False
    whole_model_ranking: str = ""
    ranking_protein_fit: bool = False
    ranking_multisample_rna: bool = False
    ranking_directional_protein_fit: bool = False
    active_fusion: bool = False
    fusion_center_precedence: bool = False
    fusion_witnesses: set[str] = field(default_factory=set)
    fusion_spanning_introns: set[tuple[int, int]] = field(default_factory=set)
    balanced_base_selected: bool = False
    base_supplement: bool = False
    decision: str = "REJECT_EVIDENCE"
    blockers: list[str] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        try:
            validate_repeat_scope(self.row)
        except ValueError as error:
            raise SelectionError(str(error)) from error

    @property
    def key(self) -> str:
        return self.row["candidate_id"]

    @property
    def seqid(self) -> str:
        return self.row["seqid"]

    @property
    def strand(self) -> str:
        return self.row["strand"]

    @property
    def start(self) -> int:
        return _integer(self.row, "start")

    @property
    def end(self) -> int:
        return _integer(self.row, "end")

    @property
    def intervals(self) -> tuple[tuple[int, int], ...]:
        return tuple((part.start, part.end) for part in self.cds)

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        try:
            values = json.loads(self.row["introns_json"])
            return tuple((int(value[0]), int(value[1])) for value in values)
        except (KeyError, TypeError, ValueError, IndexError, json.JSONDecodeError) as error:
            raise SelectionError(f"invalid introns_json for {self.key}") from error

    @property
    def recurrent_introns(self) -> frozenset[tuple[int, int]]:
        hq = _json_ints(self.row, "rna_hq_read_counts")
        samples = _json_ints(self.row, "rna_sample_hq_ge1_counts")
        if len(hq) != len(self.introns) or len(samples) != len(self.introns):
            raise SelectionError(f"RNA evidence vector length mismatch for {self.key}")
        return frozenset(
            intron for intron, count, sample_count in zip(self.introns, hq, samples, strict=True)
            if count >= 2 and sample_count >= 2
        )

    @property
    def strong_introns(self) -> frozenset[tuple[int, int]]:
        hq = _json_ints(self.row, "rna_hq_read_counts")
        replicates = _json_ints(self.row, "rna_condition_ge2rep_hq_ge1_counts")
        if len(hq) != len(self.introns) or len(replicates) != len(self.introns):
            raise SelectionError(f"RNA strong vector length mismatch for {self.key}")
        return frozenset(
            intron for intron, count, replicate in zip(self.introns, hq, replicates, strict=True)
            if count >= 3 and replicate >= 1
        )


@dataclass(frozen=True, slots=True)
class Edge:
    locus_id: str
    seqid: str
    strand: str
    left: str
    right: str
    cds_overlap_bp: int
    shorter_fraction: float
    shared_intron_count: int
    exact_intron_chain: bool
    hard_overlap: bool
    relation: str


def load_component_rows(path: Path) -> list[Record]:
    records: list[Record] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = COMPONENT_REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise SelectionError(f"component table columns missing: {sorted(missing)}")
        with DiskBackedUniqueIds() as unique_ids:
            for number, raw in enumerate(reader, 2):
                row = {key: (value or "").strip() for key, value in raw.items()}
                key = row["candidate_id"]
                if not key or not unique_ids.add(key):
                    raise SelectionError(
                        f"duplicate/empty candidate at {path}:{number}"
                    )
                records.append(
                    Record(row=row, cds=_parse_cds(row["cds_signature"]))
                )
    if not records:
        raise SelectionError(f"component evidence table is empty: {path}")
    return records


def iter_component_blocks(path: Path) -> Iterator[list[Record]]:
    """Stream exact independent span-overlap blocks in coordinate order."""

    observed = 0
    block: list[Record] = []
    block_seqid = ""
    block_end = 0
    previous_key: tuple[object, ...] | None = None
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = COMPONENT_REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise SelectionError(
                f"component table columns missing: {sorted(missing)}"
            )
        with DiskBackedUniqueIds() as unique_ids:
            for number, raw in enumerate(reader, 2):
                row = {key: (value or "").strip() for key, value in raw.items()}
                candidate_id = row["candidate_id"]
                if not candidate_id or not unique_ids.add(candidate_id):
                    raise SelectionError(
                        f"duplicate/empty candidate at {path}:{number}"
                    )
                record = Record(row=row, cds=_parse_cds(row["cds_signature"]))
                key = (
                    record.seqid, record.start, record.end,
                    record.strand, record.key,
                )
                if previous_key is not None and key < previous_key:
                    raise SelectionError(
                        "component table is not in canonical coordinate order at "
                        f"{path}:{number}"
                    )
                previous_key = key
                if block and (
                    record.seqid != block_seqid or record.start > block_end
                ):
                    yield block
                    block = []
                if not block:
                    block_seqid = record.seqid
                    block_end = record.end
                else:
                    block_end = max(block_end, record.end)
                block.append(record)
                observed += 1
    if block:
        yield block
    if not observed:
        raise SelectionError(f"component evidence table is empty: {path}")


def _stable_locus_id(assembly_id: str, members: Sequence[Record]) -> str:
    payload = {
        "schema": 1,
        "assembly_id": assembly_id,
        "seqid": members[0].seqid,
        "strand": members[0].strand,
        "members": sorted(record.key for record in members),
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return f"GLOC-{digest[:24]}"


def build_loci_and_edges(
    assembly_id: str, records: Sequence[Record]
) -> tuple[dict[str, list[Record]], list[Edge]]:
    by_axis: dict[tuple[str, str], list[Record]] = defaultdict(list)
    for record in records:
        by_axis[(record.seqid, record.strand)].append(record)
    preliminary: list[tuple[Record, Record, int, float, int, bool, bool, str]] = []
    for axis in sorted(by_axis):
        ordered = sorted(by_axis[axis], key=lambda item: (item.start, item.end, item.key))
        for index, left in enumerate(ordered):
            for right_index in range(index + 1, len(ordered)):
                right = ordered[right_index]
                if right.start > left.end:
                    break
                shared = len(set(left.introns) & set(right.introns))
                overlap = overlap_bp(left.intervals, right.intervals)
                if not shared and overlap == 0:
                    continue
                shorter = min(interval_bp(left.intervals), interval_bp(right.intervals))
                fraction = overlap / shorter if shorter else 0.0
                exact = bool(left.introns) and left.introns == right.introns
                hard = not shared and fraction >= 0.5
                relation = "SAME_CHAIN" if exact else "ALT_SPLICE" if shared else "HARD_OVERLAP" if hard else "SOFT_OVERLAP"
                preliminary.append((left, right, overlap, fraction, shared, exact, hard, relation))

    # A selection locus is a connected component of actual CDS/shared-intron
    # competition edges, not merely a transcript-span overlap component. This
    # prevents long UTR/intron spans with no coding conflict from joining loci.
    parent = {record.key: record.key for record in records}

    def find(key: str) -> str:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    for left, right, *_ in preliminary:
        union(left.key, right.key)
    grouped: dict[str, list[Record]] = defaultdict(list)
    for record in records:
        grouped[find(record.key)].append(record)
    loci: dict[str, list[Record]] = {}
    locus_by_key: dict[str, str] = {}
    for members in grouped.values():
        locus_id = _stable_locus_id(assembly_id, members)
        if locus_id in loci:
            raise SelectionError(f"locus identifier collision: {locus_id}")
        loci[locus_id] = members
        for record in members:
            record.locus_id = locus_id
            locus_by_key[record.key] = locus_id
    edges = [
        Edge(
            locus_by_key[left.key], left.seqid, left.strand, left.key, right.key,
            overlap, fraction, shared, exact, hard, relation,
        )
        for left, right, overlap, fraction, shared, exact, hard, relation in preliminary
    ]
    return loci, edges


def _recurrent_support_pass(
    supported: int, total: int, fraction: float, maximum_missing: int | None,
) -> bool:
    return (
        total > 0
        and supported / total >= fraction
        and (maximum_missing is None or total - supported <= maximum_missing)
    )


def _protein_core_pass(record: Record, minimum: float) -> bool:
    keys = (
        "protein_model_precision", "protein_chain_recall", "protein_cds_f1",
        "protein_phase_fraction",
    )
    if any(_numeric(record.row, key, -1.0) < minimum for key in keys):
        return False
    if record.introns and _numeric(record.row, "protein_intron_precision", -1.0) < minimum:
        return False
    return True


def calibrate_rna_mode(
    sample_count: int, replicates_per_condition: int, policy: Mapping[str, Any],
    *, strong_intron_count: int = 0,
) -> str:
    # Aggregate evidence records conditions with >=2 supporting replicates.
    # The minimum replicate count across unrelated conditions cannot veto it.
    if int(policy["rna"]["min_replicates_strong"]) != 2:
        raise SelectionError("junction evidence supports min_replicates_strong=2 only")
    if sample_count >= 2 and strong_intron_count > 0:
        return "REPLICATED_STRONG"
    if sample_count >= int(policy["rna"]["min_samples_recurrent"]):
        return "MULTISAMPLE_RECURRENT"
    return "DISABLED_WEAK_DESIGN"


def _noncanonical_rna_rescue(record: Record) -> bool:
    if _integer(record.row, "noncanonical_splice_count") == 0:
        return True
    motifs = tuple(json.loads(record.row["splice_motifs"]))
    hq = _json_ints(record.row, "rna_hq_read_counts")
    replicates = _json_ints(record.row, "rna_condition_ge2rep_hq_ge1_counts")
    samples = _json_ints(record.row, "rna_sample_hq_ge1_counts")
    if not (
        len(motifs) == len(record.introns)
        and len(hq) == len(record.introns)
        and len(replicates) == len(record.introns)
        and len(samples) == len(record.introns)
    ):
        return False
    full_rna = (
        hq[0] >= 2 and samples[0] >= 2
        if len(record.introns) == 1
        else bool(record.introns) and all(count >= 1 for count in hq)
    )
    if not full_rna:
        return False
    return all(
        motif in {"GT-AG", "GC-AG", "AT-AC"}
        or (hq[index] >= 3 and replicates[index] >= 1)
        for index, motif in enumerate(motifs)
    )


def rna_cds_conflict_penalty(row: Mapping[str, str], policy: Mapping[str, Any]) -> bool:
    """New policies require fresh boundary evidence; legacy tables cannot bypass it."""
    if not policy.get("rna_cds_boundary_ranking", False):
        return False
    competition = policy.get("rna_cds_competition_ranking", False)
    missing = set(RNA_CDS_COLUMNS if competition else RNA_CDS_BASE_COLUMNS) - row.keys()
    versions = {"2"} if competition else {"1", "2"}
    if missing or str(row.get("rna_cds_boundary_version")) not in versions:
        raise SelectionError("RNA CDS boundary evidence missing/stale; recompute components")
    status = row["rna_cds_status"]
    if status not in RNA_CDS_STATUSES:
        raise SelectionError(f"invalid RNA CDS status: {status}")
    if competition:
        try:
            conflicts = json.loads(row["rna_cds_unopposed_conflicts_json"])
            fraction = float(row["rna_cds_unopposed_conflict_fraction"])
        except (ValueError, TypeError) as exc:
            raise SelectionError("invalid RNA CDS competition evidence") from exc
        if (not isinstance(conflicts, list) or not 0 <= fraction <= 1
                or bool(conflicts) != (fraction > 0)):
            raise SelectionError("inconsistent RNA CDS competition evidence")
        return bool(conflicts)
    return status == "CONFLICT"


def whole_model_ranking_mode(policy: Mapping[str, Any]) -> str:
    """Require an explicit experiment and fresh RNA competition observations."""
    if "whole_model_ranking" not in policy:
        return ""
    mode = policy["whole_model_ranking"]
    if mode not in (
        "minimal", "coherent", "protein_first", "corroborated", "corroborated_directional",
        "corroborated_guarded", "corroborated_multisample",
    ):
        raise SelectionError("unknown whole_model_ranking experiment")
    if not (
        policy.get("rna_cds_boundary_ranking")
        and policy.get("rna_cds_competition_ranking")
    ):
        raise SelectionError("whole_model_ranking requires RNA CDS competition evidence")
    if mode == "corroborated_multisample":
        minimum = policy.get("single_exon_min_full_cds_sample_sources")
        if type(minimum) is not int or minimum < 2:
            raise SelectionError("multisample ranking requires explicit sample-GTF minimum >=2")
    return mode


def whole_model_rna_support(record: Record) -> tuple[bool, bool]:
    """Return witness existence, never assembly counts or replication claims."""
    try:
        full = json.loads(record.row["rna_cds_support_transcripts_json"])
        chain = json.loads(record.row["rna_cds_chain_support_transcripts_json"])
        strong_chain = _integer(record.row, "rna_cds_self_strong_chain")
    except (KeyError, TypeError, ValueError) as error:
        raise SelectionError("invalid whole-model RNA witnesses; recompute components") from error
    if (
        not isinstance(full, list) or not isinstance(chain, list)
        or not all(isinstance(item, str) and item for item in (*full, *chain))
        or strong_chain not in (0, 1)
        or bool(full) != (record.row["rna_cds_status"] in {"SUPPORTED", "MIXED"})
        or (strong_chain and (not chain or not record.introns
                              or len(record.strong_introns) != len(record.introns)))
    ):
        raise SelectionError("inconsistent whole-model RNA witnesses; recompute components")
    return bool(full), bool(strong_chain)


def multisample_rna_support_for_ranking(record: Record, policy: Mapping[str, Any]) -> bool:
    """Validate fresh full-CDS file provenance, without inferring replication."""
    if record.whole_model_ranking != "corroborated_multisample":
        return False
    try:
        version = record.row["rna_cds_sample_provenance_version"]
        sources = json.loads(record.row["rna_cds_support_sample_sources_json"])
        count = _integer(record.row, "rna_cds_support_sample_count")
    except (KeyError, TypeError, ValueError, SelectionError) as error:
        raise SelectionError("RNA CDS sample provenance missing/invalid; recompute components") from error
    if (
        str(version) != "1" or not isinstance(sources, list)
        or not all(isinstance(source, str) and source
                   and Path(source).is_absolute() and str(Path(source)) == source
                   and ".." not in Path(source).parts for source in sources)
        or sources != sorted(set(sources)) or count != len(sources)
        or (sources and not whole_model_rna_support(record)[0])
    ):
        raise SelectionError("RNA CDS sample provenance stale/inconsistent; recompute components")
    return count >= policy["single_exon_min_full_cds_sample_sources"]


def protein_fit_for_ranking(record: Record, policy: Mapping[str, Any]) -> bool:
    """Mapped-chain agreement for ranking only; never grant an admission route.

    Copy-number uniqueness, exact termini and chain completeness are separate
    observations. Single-exon concordance can require matching protein termini
    or, in the multisample experiment, enough distinct full-CDS sample-GTF sources.
    Core agreement includes phase and is not whole-reference sequence coverage.
    """
    return (
        record.whole_model_ranking in {
            "corroborated", "corroborated_directional", "corroborated_guarded",
            "corroborated_multisample",
        }
        and _integer(record.row, "protein_tier") >= int(policy["protein"]["strong_min_tier"])
        and _protein_core_pass(record, float(policy["protein"]["min_core_metric"]))
        and (
            record.whole_model_ranking not in {"corroborated_guarded", "corroborated_multisample"}
            or bool(record.introns)
            or (_integer(record.row, "protein_start_match") + _integer(record.row, "protein_stop_match")
                >= int(policy["protein"]["required_terminal_matches"]))
            or record.ranking_multisample_rna
        )
    )


def directional_protein_fit_for_ranking(record: Record, policy: Mapping[str, Any]) -> bool:
    """Candidate-side mapped-chain fit, not whole-reference coverage or admission.

    Only the directional experiment can combine this fit with full-CDS RNA,
    a self-strong coding chain and strong recurrent RNA. Reference-side recall
    and F1 remain required by the existing corroborated path.
    """
    if record.whole_model_ranking != "corroborated_directional":
        return False
    keys = ["protein_model_precision", "protein_phase_fraction"]
    if record.introns:
        keys.append("protein_intron_precision")
    return (
        _integer(record.row, "protein_tier") >= int(policy["protein"]["moderate_min_tier"])
        and all(_numeric(record.row, key, -1.0) >= float(policy["protein"]["min_core_metric"])
                for key in keys)
    )


def _dual_evidence_rank(record: Record) -> bool:
    legacy = {"STRONG_RECURRENT_RNA", "STRONG_UNIQUE_FULL_PROTEIN"} <= set(record.routes)
    if not record.whole_model_ranking or record.whole_model_ranking == "protein_first":
        return legacy
    full, strong_chain = whole_model_rna_support(record)
    if record.whole_model_ranking in {
        "corroborated", "corroborated_directional", "corroborated_guarded",
        "corroborated_multisample",
    }:
        return (
            legacy or (record.ranking_protein_fit and (full or strong_chain))
            or (record.ranking_directional_protein_fit
                and full and strong_chain and record.strong_rna)
        )
    supported = record.strong_unique_full_protein and (
        full or (record.strong_rna and strong_chain)
    )
    return legacy or supported if record.whole_model_ranking == "minimal" else supported


def add_whole_model_ranking_reasons(record: Record) -> None:
    if not record.whole_model_ranking or not record.eligible:
        return
    full, strong_chain = whole_model_rna_support(record)
    record.reasons.append("WHOLE_MODEL_RANKING_" + record.whole_model_ranking.upper())
    if record.strong_rna and not (full or strong_chain):
        record.reasons.append("RNA50_ADMISSION_WITHOUT_WHOLE_RNA_WITNESS")
    if record.whole_model_ranking in {
        "corroborated", "corroborated_directional", "corroborated_guarded",
        "corroborated_multisample",
    }:
        if record.whole_model_ranking == "corroborated_guarded" and not record.introns:
            record.reasons.append("SINGLE_EXON_CONCORDANCE_TERMINAL_GUARD")
        if record.whole_model_ranking == "corroborated_multisample" and not record.introns:
            record.reasons.append("SINGLE_EXON_CONCORDANCE_TERMINAL_OR_SAMPLE_GTF_GUARD")
        if record.ranking_protein_fit and (full or strong_chain):
            record.reasons.extend((
                "RNA_PROTEIN_CDS_CONCORDANCE",
                "PROTEIN_FIT_RANK_ONLY_ADMISSION_UNCHANGED",
                "DUAL_RANK_FULL_CDS_RNA_WITNESS" if full
                else "DUAL_RANK_STRONG_CODING_CHAIN_WITNESS",
            ))
            if record.row["protein_multiplicity_status"] == "AMBIGUOUS_PROTEIN_MULTI_HIT":
                record.reasons.append("MULTI_HIT_FIT_USED_FOR_RANK_ONLY")
            if not record.introns and full:
                record.reasons.append("SINGLE_EXON_RNA_PROTEIN_DUAL_RANK")
                if record.ranking_multisample_rna:
                    record.reasons.append("MULTIPLE_SAMPLE_GTFS_SUPPORT_FULL_CDS")
        elif record.ranking_directional_protein_fit and full and strong_chain and record.strong_rna:
            record.reasons.extend((
                "RNA_PROTEIN_DIRECTIONAL_CONCORDANCE",
                "PROTEIN_FIT_RANK_ONLY_ADMISSION_UNCHANGED",
                "DUAL_RANK_FULL_CDS_RNA_AND_STRONG_CODING_CHAIN",
                "CANDIDATE_PRECISION_RANK_ONLY_REFERENCE_RECALL_NOT_REQUIRED",
            ))
            if record.row["protein_multiplicity_status"] == "AMBIGUOUS_PROTEIN_MULTI_HIT":
                record.reasons.append("MULTI_HIT_FIT_USED_FOR_RANK_ONLY")
        return
    if (
        record.whole_model_ranking != "protein_first"
        and record.strong_unique_full_protein and (full or strong_chain)
        and _dual_evidence_rank(record)
    ):
        record.reasons.append(
            "DUAL_RANK_FULL_CDS_RNA_WITNESS" if full else "DUAL_RANK_STRONG_CODING_CHAIN_WITNESS"
        )
        if not record.introns and full:
            record.reasons.append("SINGLE_EXON_RNA_PROTEIN_DUAL_RANK")
    elif (
        record.whole_model_ranking == "coherent"
        and record.strong_rna and record.strong_unique_full_protein
    ):
        record.reasons.append("RNA_PROTEIN_ADMISSION_WITHOUT_WHOLE_RNA_DUAL_RANK")


def fusion_whole_model_support(record: Record) -> bool:
    """Restrict automatic precedence, not the candidate's admission."""
    if not record.whole_model_ranking:
        return True
    full, strong_chain = whole_model_rna_support(record)
    return not record.rna_cds_conflict and (full or strong_chain)


def evaluate_record(
    record: Record,
    policy: Mapping[str, Any],
    *,
    sample_count: int,
    replicates_per_condition: int,
    additional_routes: Callable[[Record], Sequence[str]] | None = None,
) -> None:
    record.whole_model_ranking = whole_model_ranking_mode(policy)
    record.ranking_multisample_rna = multisample_rna_support_for_ranking(record, policy)
    record.ranking_protein_fit = protein_fit_for_ranking(record, policy)
    record.ranking_directional_protein_fit = directional_protein_fit_for_ranking(record, policy)
    record.rna_cds_ranking = bool(policy.get("rna_cds_boundary_ranking", False))
    record.rna_cds_competition = bool(policy.get("rna_cds_competition_ranking", False))
    record.rna_cds_conflict = rna_cds_conflict_penalty(record.row, policy)
    if record.whole_model_ranking:
        whole_model_rna_support(record)
    if confirmed_only(policy):
        record.row["te_policy_id"] = policy["te_policy_id"]
    hard = policy["hard_qc"]
    intron_count = len(record.introns)
    protein_aa = _integer(record.row, "translated_protein_aa")
    record.hard_qc = (
        record.row["orf_status"] == "complete"
        and _integer(record.row, "common_qc_pass") == int(hard["common_qc_pass"])
        and _integer(record.row, "ambiguous_splice_count") == 0
        and _noncanonical_rna_rescue(record)
        and protein_aa >= int(hard["min_protein_aa"])
        and (intron_count > 0 or protein_aa >= int(hard["single_exon_min_protein_aa"]))
    )
    strong_count = len(record.strong_introns)
    record.rna_mode = calibrate_rna_mode(
        sample_count, replicates_per_condition, policy, strong_intron_count=strong_count
    )
    recurrent_count = len(record.recurrent_introns)
    rna_policy = policy["rna"]
    maximum_missing = rna_policy["max_nonrecurrent_introns"]
    if maximum_missing is not None:
        maximum_missing = int(maximum_missing)
    record.strong_rna = (
        record.rna_mode == "REPLICATED_STRONG"
        and _integer(record.row, "rna_tier") >= int(rna_policy["strong_min_tier"])
        and _recurrent_support_pass(strong_count, intron_count, float(rna_policy["min_recurrent_fraction"]), maximum_missing)
    )
    record.recurrent_rna = (
        record.rna_mode != "DISABLED_WEAK_DESIGN"
        and _integer(record.row, "rna_tier") >= int(rna_policy["moderate_min_tier"])
        and _recurrent_support_pass(recurrent_count, intron_count, float(rna_policy["min_recurrent_fraction"]), maximum_missing)
    )

    uniqueness = record.row["protein_multiplicity_status"] in {
        "UNIQUE_PROTEIN_HIT", "UNIQUE_RETAINED_PROTEIN_HIT_SEED_MULTI_HIT_AUDIT"
    }
    protein_policy = policy["protein"]
    protein_tier = _integer(record.row, "protein_tier")
    record.moderate_unique_protein = uniqueness and protein_tier >= int(protein_policy["moderate_min_tier"])
    record.strong_unique_protein = (
        uniqueness
        and protein_tier >= int(protein_policy["strong_min_tier"])
        and _protein_core_pass(record, float(protein_policy["min_core_metric"]))
    )
    terminal_count = _integer(record.row, "protein_start_match") + _integer(record.row, "protein_stop_match")
    record.strong_unique_full_protein = (
        record.strong_unique_protein
        and terminal_count >= int(protein_policy["required_terminal_matches"])
    )
    record.predictor_exact = (
        _integer(record.row, "predictor_tier") >= int(policy["predictor"]["exact_min_tier"])
        and (
            not bool(policy["predictor"]["require_both_predictors"])
            or {"augustus", "genemark_etp"} <= set(record.row["predictors"].split(","))
        )
    )
    repeat_fraction = _numeric(record.row, "repeat_cds_fraction")
    record.repeat_risk = (
        "SEVERE" if repeat_fraction >= float(policy["repeat"]["severe_cds_fraction"])
        else "HIGH" if repeat_fraction >= float(policy["repeat"]["high_cds_fraction"])
        else "LOW"
    )
    if not record.hard_qc:
        record.reasons.append("HARD_QC_FAIL")
        return
    if intron_count == 0:
        if record.strong_unique_full_protein:
            record.routes.append("STRONG_UNIQUE_FULL_PROTEIN")
    else:
        if record.strong_rna:
            record.routes.append("STRONG_RECURRENT_RNA")
        if record.recurrent_rna:
            record.routes.append(f"RECURRENT_RNA_{100 * float(rna_policy['min_recurrent_fraction']):g}")
        if record.strong_unique_full_protein:
            record.routes.append("STRONG_UNIQUE_FULL_PROTEIN")
        if record.recurrent_rna and record.moderate_unique_protein:
            record.routes.append("RECURRENT_RNA_MODERATE_UNIQUE_PROTEIN")
        if record.predictor_exact and record.recurrent_rna:
            record.routes.append("EXACT_PREDICTOR_RECURRENT_RNA")
        if record.predictor_exact and record.strong_unique_protein:
            record.routes.append("EXACT_PREDICTOR_STRONG_UNIQUE_PROTEIN")
        if record.predictor_exact and record.moderate_unique_protein and bool(policy["predictor"].get("allow_exact_moderate_unique_route", False)):
            record.routes.append("EXACT_PREDICTOR_MODERATE_UNIQUE_PROTEIN")
        if record.strong_unique_protein and bool(protein_policy.get("allow_strong_unique_core_route", False)):
            record.routes.append("STRONG_UNIQUE_CORE_PROTEIN")
        long_gene = rna_policy.get("long_gene_evidence_rescue", {})
        if (
            bool(long_gene.get("enabled"))
            and intron_count >= int(long_gene["min_introns"])
            and recurrent_count / intron_count >= float(long_gene["min_recurrent_fraction"])
            and record.moderate_unique_protein
            and _numeric(record.row, "protein_chain_recall") >= float(long_gene["min_protein_chain_recall"])
            and (repeat_neutral(record) or record.repeat_risk == "LOW")
        ):
            record.routes.append("LONG_GENE_RNA_PROTEIN")
    if additional_routes is not None:
        record.routes.extend(additional_routes(record))
    if not record.routes:
        record.reasons.append("NO_INDEPENDENT_EVIDENCE_ROUTE")
        return
    if repeat_neutral(record):
        record.repeat_override_reason = "OVERLAP_REPORT_ONLY"
    elif record.repeat_risk == "LOW":
        record.repeat_override_reason = "LOW_CLASSIFIED_TE"
    elif record.strong_unique_full_protein:
        record.repeat_override_reason = "STRONG_UNIQUE_FULL_PROTEIN"
    elif record.strong_rna:
        # Restore the legacy RNA rescue without changing classified TE scope,
        # protein uniqueness, or downstream graph and isoform requirements.
        record.repeat_override_reason = "STRONG_RECURRENT_RNA"
    elif record.repeat_risk == "HIGH" and record.recurrent_rna and record.moderate_unique_protein:
        record.repeat_override_reason = "RECURRENT_RNA_MODERATE_UNIQUE_PROTEIN"
    else:
        record.repeat_override_reason = "INSUFFICIENT_INDEPENDENT_PROTEIN_OR_RNA"
    record.repeat_override = record.repeat_override_reason != "INSUFFICIENT_INDEPENDENT_PROTEIN_OR_RNA"
    if not record.repeat_override:
        record.reasons.append(f"{record.repeat_risk}_REPEAT_WITHOUT_OVERRIDE")
        return
    record.eligible = True
    record.decision = "ELIGIBLE_GRAPH"
    record.reasons.append("EVIDENCE_PASS" if repeat_neutral(record) else "EVIDENCE_AND_REPEAT_GATE_PASS")
    if record.rna_cds_conflict:
        record.reasons.append("RNA_CDS_CONFLICT_DOWNRANKED")
    if record.rna_cds_competition and json.loads(record.row["rna_cds_countered_conflicts_json"]):
        record.reasons.append("RNA_CDS_CONFLICT_COUNTERED_BY_OWN_RNA")
    if record.rna_cds_ranking and record.row["rna_cds_status"] in {"SUPPORTED", "MIXED"}:
        record.reasons.append("RNA_CDS_COMPLETE_STRUCTURE_SUPPORT")
    add_whole_model_ranking_reasons(record)


def _evidence_rank(record: Record) -> tuple[object, ...]:
    routes = set(record.routes)
    both = _dual_evidence_rank(record)
    route_rank = (
        10 if record.fusion_center_precedence else
        9 if both else 8 if routes & {"STRONG_RECURRENT_RNA", "STRONG_UNIQUE_FULL_PROTEIN"}
        else 7 if routes & {"RECURRENT_RNA_MODERATE_UNIQUE_PROTEIN", "LONG_GENE_RNA_PROTEIN"}
        else 6 if routes & {"STRONG_UNIQUE_CORE_PROTEIN", "EXACT_PREDICTOR_MODERATE_UNIQUE_PROTEIN"}
        else 5 if record.recurrent_rna else 4
    )
    protein_unique = record.row["protein_multiplicity_status"] in {
        "UNIQUE_PROTEIN_HIT", "UNIQUE_RETAINED_PROTEIN_HIT_SEED_MULTI_HIT_AUDIT"
    }
    rank = (
        -route_rank, -int(record.strong_rna), -int(record.strong_unique_full_protein),
        -int(record.predictor_exact), -_integer(record.row, "rna_tier"),
        -(_integer(record.row, "protein_tier") if protein_unique else 0),
        -((_integer(record.row, "protein_start_match") + _integer(record.row, "protein_stop_match")) if protein_unique else 0),
        -(_numeric(record.row, "protein_cds_f1") if protein_unique else 0.0),
        0.0 if repeat_neutral(record) else _numeric(record.row, "repeat_cds_fraction"),
        -_integer(record.row, "translated_protein_aa"), record.key,
    )
    if record.whole_model_ranking in {"minimal", "coherent", "protein_first"}:
        rank = (rank[0], rank[2], rank[1], *rank[3:])
    if record.rna_cds_ranking:
        # Explicit contradictions rank below otherwise eligible alternatives.
        # Preserve the established evidence-route hierarchy; a short complete
        # RNA is not by itself stronger than a better RNA/protein evidence class.
        complete = record.row["rna_cds_status"] in {"SUPPORTED", "MIXED"}
        if record.rna_cds_competition:
            # Evidence routes precede unresolved contradictions: absence of RNA
            # is not a stronger evidence class. Keep positive isoform evidence
            # and terminal conflicts separate, then use established tie-breaks.
            # No candidate ID, protein length or locus exceptions.
            return (
                rank[0], int(record.rna_cds_conflict), -int(complete),
                *rank[1:],
            )
        return (int(record.rna_cds_conflict), rank[0], -int(complete), *rank[1:])
    return rank


def _unique_isoform_proof(candidate: Record, anchors: Iterable[Record]) -> str | None:
    proofs: set[str] = set()
    for anchor in anchors:
        unique = set(candidate.introns) - set(anchor.introns)
        rna = bool(unique) and len(candidate.recurrent_introns) == len(candidate.introns)
        protein = (
            bool(unique) and candidate.strong_unique_protein
            and _numeric(candidate.row, "protein_intron_precision") >= 1.0
            and _numeric(candidate.row, "protein_intron_recall") >= 0.8
        )
        if not (rna or protein):
            return None
        if rna:
            proofs.add("RNA")
        if protein:
            proofs.add("PROTEIN")
    if proofs == {"RNA", "PROTEIN"}:
        return "RNA_AND_PROTEIN"
    return next(iter(proofs)) if proofs else None


def _detect_active_fusions(
    members: Sequence[Record], relations: Mapping[tuple[str, str], str]
) -> None:
    by_key = {record.key: record for record in members}
    for center in members:
        if not center.eligible:
            continue
        neighbors = [
            other for other in members
            if other.key != center.key and other.eligible
            and tuple(sorted((center.key, other.key))) in relations
        ]
        for index, left in enumerate(neighbors):
            for right_index in range(index + 1, len(neighbors)):
                right = neighbors[right_index]
                if tuple(sorted((left.key, right.key))) in relations:
                    continue
                first, second = sorted((left, right), key=lambda item: item.start)
                bridges = {
                    intron for intron in center.introns
                    if first.end < intron[0] and intron[1] < second.start
                }
                if center.start <= first.start and center.end >= second.end and bridges:
                    center.active_fusion = True
                    center.fusion_witnesses.update((left.key, right.key))
                    center.fusion_spanning_introns.update(bridges)
                    strong_bridge = bool(bridges & set(center.strong_introns))
                    if not (strong_bridge or center.strong_unique_full_protein):
                        center.eligible = False
                        center.decision = "REJECT_ACTIVE_FUSION_CENTER"
                        center.reasons.append("ACTIVE_FUSION_REQUIRES_STRONG_SPANNING_RNA_OR_FULL_PROTEIN")


def _frozen_base_alt_support(
    candidate: Record,
    core_alt_anchors: Sequence[Record],
    policy: Mapping[str, Any],
) -> str | None:
    """Return the bounded current-data exception for a base ALT model."""
    if candidate.strong_rna:
        return "STRONG_RNA"
    if candidate.strong_unique_full_protein:
        return "FULL_PROTEIN"
    support = policy["graph"].get("base_alt_current_support", {})
    intron_count = len(candidate.introns)
    recurrent = len(candidate.recurrent_introns)
    if (
        intron_count == 0
        or recurrent / intron_count
        < float(support.get("min_recurrent_fraction", 0.8))
        or intron_count - recurrent
        > int(support.get("max_nonrecurrent_introns", 2))
    ):
        return None
    minimum_unique = int(
        support.get("relative_rna_min_unique_junctions_per_core_alt", 1)
    )
    if all(
        len(set(candidate.recurrent_introns) - set(anchor.introns)) >= minimum_unique
        for anchor in core_alt_anchors
    ):
        return "RELATIVE_RECURRENT_RNA"
    terminal_count = (
        _integer(candidate.row, "protein_start_match")
        + _integer(candidate.row, "protein_stop_match")
    )
    protein_unique = candidate.row["protein_multiplicity_status"] in {
        "UNIQUE_PROTEIN_HIT",
        "UNIQUE_RETAINED_PROTEIN_HIT_SEED_MULTI_HIT_AUDIT",
    }
    target_linked = _integer(candidate.row, "protein_target_linked") == 1
    if (
        intron_count >= int(support.get("long_gene_min_introns", 10))
        and terminal_count
        >= int(support.get("long_gene_required_terminal_matches", 2))
        and _numeric(candidate.row, "protein_cds_f1")
        >= float(support.get("long_gene_min_protein_cds_f1", 0.85))
        and protein_unique
        and target_linked
    ):
        return "LONG_GENE_RNA_PROTEIN"
    return None


def evaluate_and_select(
    loci: Mapping[str, Sequence[Record]],
    edges: Sequence[Edge],
    policy: Mapping[str, Any],
    *,
    sample_count: int,
    replicates_per_condition: int,
    base_selected_keys: frozenset[str] | set[str] | None = None,
) -> list[Record]:
    neutral = confirmed_only(policy)
    if neutral:
        loci, edges = partition_pools(loci, edges)
    policy_id = str(policy["policy_id"])
    verified_base = frozenset(base_selected_keys or ())
    by_locus_edges: dict[str, list[Edge]] = defaultdict(list)
    for edge in edges:
        by_locus_edges[edge.locus_id].append(edge)
    selected_all: list[Record] = []
    for locus_id, members in sorted(loci.items()):
        for record in members:
            evaluate_record(
                record, policy, sample_count=sample_count,
                replicates_per_condition=replicates_per_condition,
            )
            record.row["selection_policy_id"] = policy_id
            record.balanced_base_selected = record.key in verified_base
            if (
                record.balanced_base_selected
                and record.hard_qc
                and (repeat_neutral(record) or record.repeat_risk == "LOW")
            ):
                record.routes.append("FROZEN_BALANCED_BASE")
                record.repeat_override = True
                record.repeat_override_reason = "OVERLAP_REPORT_ONLY" if neutral else "LOW_CLASSIFIED_TE"
                record.eligible = True
                record.decision = "ELIGIBLE_GRAPH"
                record.reasons.append("RECOMPUTED_BALANCED_BASE" if neutral else "FROZEN_BALANCED_BASE_LOW_REPEAT")
        relations = {
            tuple(sorted((edge.left, edge.right))): edge.relation
            for edge in by_locus_edges.get(locus_id, ())
        }
        _detect_active_fusions(members, relations)
        precedence = policy["graph"].get("fusion_center_precedence", {})
        for record in members:
            if (
                record.active_fusion
                and record.eligible
                and record.hard_qc
                and len(record.introns) >= int(precedence.get("minimum_intron_count", 2))
                and bool(record.fusion_spanning_introns)
                and record.fusion_spanning_introns <= set(record.strong_introns)
                and len(record.strong_introns) == len(record.introns)
                and len(record.recurrent_introns) == len(record.introns)
                and (repeat_neutral(record) or record.repeat_risk == "LOW")
                and record.balanced_base_selected
                and record.moderate_unique_protein
                and _numeric(record.row, "protein_chain_recall")
                >= float(precedence.get("minimum_protein_chain_recall", 0.95))
            ):
                if fusion_whole_model_support(record):
                    record.fusion_center_precedence = True
                    record.reasons.append("FUSION_CENTER_STRICT_PRECEDENCE")
                else:
                    record.reasons.append("FUSION_PRECEDENCE_WITHHELD_WHOLE_RNA_OR_CONFLICT")
        selected: list[Record] = []
        for candidate in sorted(
            (
                record for record in members
                if record.eligible
                and any(route != "FROZEN_BALANCED_BASE" for route in record.routes)
            ),
            key=_evidence_rank,
        ):
            incompatible = [
                anchor for anchor in selected
                if relations.get(tuple(sorted((candidate.key, anchor.key)))) in {"SAME_CHAIN", "HARD_OVERLAP"}
                or anchor.key in candidate.fusion_witnesses
                or candidate.key in anchor.fusion_witnesses
            ]
            if incompatible:
                candidate.decision = "REJECT_INCOMPATIBLE_GRAPH"
                candidate.blockers = sorted(anchor.key for anchor in incompatible)
                candidate.reasons.append("LOWER_RANK_INCOMPATIBLE_MODEL")
                continue
            alt_anchors = [
                anchor for anchor in selected
                if relations.get(tuple(sorted((candidate.key, anchor.key)))) == "ALT_SPLICE"
            ]
            if alt_anchors:
                proof = _unique_isoform_proof(candidate, alt_anchors)
                if proof is None:
                    candidate.decision = "REJECT_ALT_ISOFORM_NO_UNIQUE_EVIDENCE"
                    candidate.blockers = sorted(anchor.key for anchor in alt_anchors)
                    candidate.reasons.append("ALT_ISOFORM_WITHOUT_UNIQUE_EVIDENCE")
                    continue
                candidate.reasons.append(f"ALT_ISOFORM_UNIQUE_{proof}")
            candidate.selected = True
            candidate.decision = "SELECTED"
            selected.append(candidate)

        # The current evidence graph is the immutable core.  A same-run
        # balanced-base model may supplement it, but can never displace or
        # conflict with a core model.  Relations among verified base supplements
        # are retained as an explicit, audited grandfathered condition.
        core_selected = frozenset(record.key for record in selected)
        for candidate in sorted(members, key=_evidence_rank):
            if (
                candidate.selected
                or not candidate.eligible
                or not candidate.balanced_base_selected
                or "FROZEN_BALANCED_BASE" not in candidate.routes
            ):
                continue
            if candidate.active_fusion:
                candidate.reasons.append(
                    "ACTIVE_FUSION_NOT_ELIGIBLE_FOR_FROZEN_BASE_SUPPLEMENT"
                )
                continue
            core_blockers = [
                anchor for anchor in selected
                if anchor.key in core_selected
                and (
                    relations.get(tuple(sorted((candidate.key, anchor.key))))
                    in {"SAME_CHAIN", "HARD_OVERLAP"}
                    or anchor.key in candidate.fusion_witnesses
                    or candidate.key in anchor.fusion_witnesses
                )
            ]
            if core_blockers:
                candidate.decision = "REJECT_INCOMPATIBLE_GRAPH"
                candidate.blockers = sorted(anchor.key for anchor in core_blockers)
                candidate.reasons.append("FROZEN_BASE_BLOCKED_BY_EVIDENCE_CORE")
                continue
            alt_core = [
                anchor for anchor in selected
                if anchor.key in core_selected
                and relations.get(tuple(sorted((candidate.key, anchor.key)))) == "ALT_SPLICE"
            ]
            proof = _unique_isoform_proof(candidate, alt_core) if alt_core else None
            support = None
            if alt_core and proof is None:
                support = _frozen_base_alt_support(candidate, alt_core, policy)
                if support is None:
                    candidate.decision = "REJECT_ALT_ISOFORM_NO_UNIQUE_EVIDENCE"
                    candidate.blockers = sorted(anchor.key for anchor in alt_core)
                    candidate.reasons.append(
                        "FROZEN_BASE_ALT_WITHOUT_CURRENT_STRONG_SUPPORT"
                    )
                    continue
            grandfathered = [
                anchor.key for anchor in selected
                if anchor.key not in core_selected
                and relations.get(tuple(sorted((candidate.key, anchor.key))))
                in {"SAME_CHAIN", "HARD_OVERLAP"}
            ]
            if grandfathered:
                candidate.reasons.append(
                    "VERIFIED_BASE_INCOMPATIBLE_GRANDFATHER:"
                    + ",".join(sorted(grandfathered))
                )
            candidate.selected = True
            candidate.base_supplement = True
            candidate.decision = "SELECTED"
            candidate.blockers = []
            candidate.reasons.append("FROZEN_BASE_SUPPLEMENT_ADDED")
            if proof:
                candidate.reasons.append(f"FROZEN_BASE_ALT_UNIQUE_{proof}")
            elif support:
                candidate.reasons.append(f"FROZEN_BASE_ALT_{support}_EXCEPTION")
            selected.append(candidate)
        selected_all.extend(selected)
    if neutral:
        for record in selected_all:
            if record.row.get("candidate_pool") == "TE":
                record.row["te_selected"] = "1"
                record.selected = False
                record.decision = "TE_SELECTED"
        selected_all = [record for record in selected_all if record.selected]
    return sorted(selected_all, key=lambda item: (item.seqid, item.start, item.end, item.strand, item.key))


def audit_row(record: Record, policy_id: str) -> dict[str, object]:
    recurrent = record.recurrent_introns
    introns = record.introns
    return {
        **{key: record.row.get(key, "") for key in TE_COLUMNS},
        **{key: record.row.get(key, "NA") for key in RNA_CDS_COLUMNS + RNA_CDS_SAMPLE_COLUMNS},
        "candidate_id": record.key, "structure_id": record.row["structure_id"],
        "seqid": record.seqid, "strand": record.strand,
        "start": record.start, "end": record.end, "locus_id": record.locus_id,
        "policy_id": policy_id, "hard_qc_pass": int(record.hard_qc),
        "rna_calibration_mode": record.rna_mode,
        "strong_recurrent_rna": int(record.strong_rna),
        "recurrent_rna": int(record.recurrent_rna),
        "strong_unique_protein": int(record.strong_unique_protein),
        "strong_unique_full_protein": int(record.strong_unique_full_protein),
        "moderate_unique_protein": int(record.moderate_unique_protein),
        "predictor_exact": int(record.predictor_exact), "repeat_risk": record.repeat_risk,
        "repeat_scope": record.row["repeat_scope"],
        "repeat_cds_fraction": record.row["repeat_cds_fraction"],
        "all_repeat_cds_fraction": record.row.get("all_repeat_cds_fraction", "NA"),
        "strong_intron_count": len(record.strong_introns),
        "intron_count": len(record.introns),
        "rna_recurrent_intron_count": len(recurrent),
        "rna_recurrent_intron_fraction": format(len(recurrent) / len(introns), ".10g") if introns else "NA",
        "rna_nonrecurrent_introns_json": json.dumps([intron for intron in introns if intron not in recurrent]),
        "repeat_override_reason": record.repeat_override_reason,
        "repeat_override": int(record.repeat_override),
        "balanced_base_selected": int(record.balanced_base_selected),
        "base_supplement": int(record.base_supplement),
        "eligible": int(record.eligible),
        "rna_cds_conflict_penalty": int(record.rna_cds_conflict),
        "rna_cds_complete_support": int(
            record.rna_cds_ranking and record.row["rna_cds_status"] in {"SUPPORTED", "MIXED"}
        ),
        "selection_rank_json": json.dumps(_evidence_rank(record), separators=(",", ":")),
        "selected": int(record.selected), "active_fusion": int(record.active_fusion),
        "fusion_center_precedence": int(record.fusion_center_precedence),
        "routes": ",".join(record.routes) or "NA", "decision": record.decision,
        "blockers": ",".join(record.blockers) or "NA",
        "reasons": ",".join(record.reasons) or "NA",
    }


def edge_row(edge: Edge) -> dict[str, object]:
    return {
        "locus_id": edge.locus_id, "seqid": edge.seqid, "strand": edge.strand,
        "left_candidate_id": edge.left, "right_candidate_id": edge.right,
        "cds_overlap_bp": edge.cds_overlap_bp,
        "shorter_cds_overlap_fraction": format(edge.shorter_fraction, ".10g"),
        "shared_intron_count": edge.shared_intron_count,
        "exact_intron_chain": int(edge.exact_intron_chain),
        "hard_overlap": int(edge.hard_overlap), "relation": edge.relation,
    }


def write_dict_table(
    path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
