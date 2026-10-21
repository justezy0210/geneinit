"""One admission policy and one conflict resolver, with local invalidation.

This is an experimental successor to v2, not a promise of identical membership.
The supplemental evidence predicates are retained; the second selector and its
grandfathered conflicts are not. All candidates, including structural repairs,
pass the same admission and conflict rules.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
from typing import Mapping, Sequence

from geneinit.evidence.metrics import IntervalIndex, overlap_bp
from geneinit.evidence.te_classification import confirmed_only, partition_pools
from geneinit.selection.balanced import (
    _balanced_routes,
    _core_metrics,
    _full_support,
    _protein_predictor_rescue,
    _protein_strong,
    _protein_supported_terminal,
    _rna_predictor_complement,
    _rna_protein_complement,
)
from geneinit.selection.core import (
    Edge, Record, SelectionError, _evidence_rank, _frozen_base_alt_support,
    _integer, _numeric, _unique_isoform_proof, build_loci_and_edges,
    add_whole_model_ranking_reasons, evaluate_record, fusion_whole_model_support,
    directional_protein_fit_for_ranking, multisample_rna_support_for_ranking,
    protein_fit_for_ranking, rna_cds_conflict_penalty,
    whole_model_ranking_mode, whole_model_rna_support,
)


ENGINE = "unified_v1"


def _fingerprint(record: Record) -> str:
    row = dict(record.row)
    row.pop("selection_policy_id", None)
    row.pop("te_selected", None)
    return hashlib.sha256(json.dumps(
        [row, [(p.start, p.end, p.phase) for p in record.cds]],
        sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def supplemental_routes(record: Record) -> list[str]:
    """Evidence predicates only: no second QC, graph, selection, or fallback."""
    if not record.introns:
        terminals = _integer(record.row, "protein_start_match") + _integer(
            record.row, "protein_stop_match"
        )
        if _protein_strong(record) and terminals >= 1:
            return ["SUPPORTED_SINGLE_PROTEIN"]
        if _core_metrics(record, 0.9) and _integer(record.row, "predictor_tier") == 3:
            return ["SUPPORTED_SINGLE_PREDICTOR"]
        quality = min(
            _numeric(record.row, key, -1.0)
            for key in ("protein_model_precision", "protein_chain_recall",
                        "protein_cds_f1", "protein_phase_fraction")
        )
        if (
            (terminals >= 1 and quality >= 0.7)
            or (terminals == 0 and quality >= 0.85
                and _integer(record.row, "predictor_tier") == 3)
        ):
            return ["SUPPORTED_SINGLE_TERMINAL"]
        return []
    routes = []
    if _full_support(record):
        routes.append("SUPPORTED_FULL_RNA")
    if _protein_strong(record):
        routes.append("SUPPORTED_CORE_PROTEIN")
    elif _protein_supported_terminal(record):
        routes.append("SUPPORTED_TERMINAL_PROTEIN")
    elif _protein_predictor_rescue(record):
        routes.append("SUPPORTED_PREDICTOR_PROTEIN")
    if _rna_predictor_complement(record):
        routes.append("SUPPORTED_RNA_PREDICTOR")
    if _rna_protein_complement(record):
        routes.append("SUPPORTED_RNA_PROTEIN")
    routes.extend("SUPPORTED_" + route for route in _balanced_routes(record))
    return routes


def _mark_fusions(members: Sequence[Record], relations: Mapping) -> None:
    """Evaluate witnesses from a fixed eligible snapshot, not mutation order."""
    eligible = sorted((record for record in members if record.eligible),
                      key=lambda record: record.key)
    rejected = []
    for center in eligible:
        neighbors = [
            other for other in eligible
            if other.key != center.key
            and tuple(sorted((center.key, other.key))) in relations
        ]
        for index, left in enumerate(neighbors):
            for right in neighbors[index + 1:]:
                if tuple(sorted((left.key, right.key))) in relations:
                    continue
                first, second = sorted((left, right), key=lambda item: item.start)
                bridges = {
                    intron for intron in center.introns
                    if first.end < intron[0] and intron[1] < second.start
                }
                if center.start > first.start or center.end < second.end or not bridges:
                    continue
                center.active_fusion = True
                center.fusion_witnesses.update((left.key, right.key))
                center.fusion_spanning_introns.update(bridges)
                if not (bridges & set(center.strong_introns)
                        or center.strong_unique_full_protein):
                    rejected.append(center)
    for center in rejected:
        center.eligible = False
        center.decision = "REJECT_ACTIVE_FUSION_CENTER"
        center.reasons.append(
            "ACTIVE_FUSION_REQUIRES_STRONG_SPANNING_RNA_OR_FULL_PROTEIN"
        )


def resolve_loci(loci: Mapping[str, Sequence[Record]], edges: Sequence[Edge],
                 policy: Mapping) -> list[Record]:
    """Resolve every candidate under the same ranking and conflict rules."""
    loci, edges = partition_pools(loci, edges)
    by_locus = {}
    for edge in edges:
        by_locus.setdefault(edge.locus_id, []).append(edge)
    selected_all = []
    for locus_id, members in sorted(loci.items()):
        relations = {
            tuple(sorted((edge.left, edge.right))): edge.relation
            for edge in by_locus.get(locus_id, ())
        }
        _mark_fusions(members, relations)
        precedence = policy["graph"]["fusion_center_precedence"]
        for record in members:
            if (
                record.active_fusion and record.eligible
                and len(record.introns) >= precedence["minimum_intron_count"]
                and record.fusion_spanning_introns
                and record.fusion_spanning_introns <= set(record.strong_introns)
                and len(record.strong_introns) == len(record.introns)
                and len(record.recurrent_introns) == len(record.introns)
                and record.moderate_unique_protein
                and _numeric(record.row, "protein_chain_recall")
                >= precedence["minimum_protein_chain_recall"]
            ):
                if fusion_whole_model_support(record):
                    record.fusion_center_precedence = True
                    record.reasons.append("FUSION_CENTER_STRICT_PRECEDENCE")
                else:
                    record.reasons.append("FUSION_PRECEDENCE_WITHHELD_WHOLE_RNA_OR_CONFLICT")
        selected = []
        for candidate in sorted((r for r in members if r.eligible),
                                key=_evidence_rank):
            blockers = [
                anchor for anchor in selected
                if relations.get(tuple(sorted((candidate.key, anchor.key))))
                in {"SAME_CHAIN", "HARD_OVERLAP"}
                or anchor.key in candidate.fusion_witnesses
                or candidate.key in anchor.fusion_witnesses
            ]
            if blockers:
                candidate.decision = "REJECT_INCOMPATIBLE_GRAPH"
                candidate.blockers = sorted(anchor.key for anchor in blockers)
                candidate.reasons.append("LOWER_RANK_INCOMPATIBLE_MODEL")
                continue
            anchors = [
                anchor for anchor in selected
                if relations.get(tuple(sorted((candidate.key, anchor.key)))) == "ALT_SPLICE"
            ]
            if anchors:
                proof = _unique_isoform_proof(candidate, anchors)
                # The former supplement support criterion is now available
                # uniformly, rather than being conditional on a second selector.
                if proof is None:
                    proof = _frozen_base_alt_support(candidate, anchors, policy)
                if proof is None:
                    candidate.decision = "REJECT_ALT_ISOFORM_NO_UNIQUE_EVIDENCE"
                    candidate.blockers = sorted(anchor.key for anchor in anchors)
                    candidate.reasons.append("ALT_ISOFORM_WITHOUT_UNIQUE_EVIDENCE")
                    continue
                candidate.reasons.append("ALT_ISOFORM_UNIQUE_" + proof)
            candidate.selected = True
            candidate.decision = "SELECTED"
            selected.append(candidate)
        for record in selected:
            if record.rna_cds_conflict:
                record.reasons.append("SELECTED_WITH_UNRESOLVED_RNA_CDS_CONFLICT")
            if record.row.get("candidate_pool") == "TE":
                record.row["te_selected"] = "1"
                record.selected = False
                record.decision = "TE_SELECTED"
            else:
                selected_all.append(record)
    return sorted(selected_all, key=lambda r: (r.seqid, r.start, r.end, r.strand, r.key))


@dataclass(frozen=True)
class UpdateStats:
    evaluated_candidates: int
    reused_candidates: int
    affected_candidates: int
    resolved_loci: int
    untouched_candidates: int


class UnifiedSelector:
    """Cache candidate evidence; invalidate complete old/new competition loci.

    Updates must provide fresh evidence rows for every changed structure or
    changed evidence source. A policy or global RNA-design change requires a
    new selector. A local update cannot reuse stale upstream evidence.
    """

    def __init__(self, assembly_id: str, policy: Mapping, *,
                 sample_count: int, replicates_per_condition: int):
        whole_model_ranking_mode(policy)
        if "rna_cds_boundary_review" in policy:
            raise SelectionError("obsolete RNA CDS review-only policy; use boundary ranking policy")
        if policy.get("rna_cds_competition_ranking") and not policy.get("rna_cds_boundary_ranking"):
            raise SelectionError("RNA CDS competition requires boundary ranking evidence")
        if not confirmed_only(policy):
            raise SelectionError("unified_v1 requires confirmed-only TE policy")
        if policy.get("engine") != ENGINE:
            raise SelectionError("unified_v1 requires an explicit experimental policy")
        self.assembly_id = assembly_id
        self._policy = deepcopy(dict(policy))
        self.sample_count = sample_count
        self.replicates_per_condition = replicates_per_condition
        self.records: dict[str, Record] = {}
        self._evidence: dict[str, tuple[str, Record]] = {}
        self._members: dict[str, set[str]] = {}
        self._index = None
        self.last_loci = {}
        self.last_edges = []
        self.stats = UpdateStats(0, 0, 0, 0, 0)

    @property
    def selected(self) -> list[Record]:
        return sorted((r for r in self.records.values() if r.selected),
                      key=lambda r: (r.seqid, r.start, r.end, r.strand, r.key))

    @property
    def selected_keys(self) -> set[str]:
        """Return membership without sorting records needed only for output."""
        return {record.key for record in self.records.values() if record.selected}

    def restore(self, records: Sequence[Record], audits: Mapping[str, Mapping]) -> None:
        """Restore this policy's completed decisions without reevaluating QC."""
        if self.records:
            raise SelectionError("restore requires an empty selector")
        if len({r.key for r in records}) != len(records):
            raise SelectionError("duplicate candidates in selector checkpoint")
        if {r.key for r in records} != set(audits):
            raise SelectionError("selector checkpoint candidate sets differ")
        boolean_fields = {
            "hard_qc": "hard_qc_pass", "strong_rna": "strong_recurrent_rna",
            "recurrent_rna": "recurrent_rna",
            "strong_unique_protein": "strong_unique_protein",
            "strong_unique_full_protein": "strong_unique_full_protein",
            "moderate_unique_protein": "moderate_unique_protein",
            "predictor_exact": "predictor_exact", "repeat_override": "repeat_override",
        }
        for source in records:
            audit = audits[source.key]
            if audit["policy_id"] != self._policy["policy_id"]:
                raise SelectionError("selector checkpoint policy mismatch")
            row = dict(source.row)
            row["te_policy_id"] = self._policy["te_policy_id"]
            row.pop("te_selected", None)
            row["selection_policy_id"] = self._policy["policy_id"]
            pristine = Record(row=row, cds=source.cds)
            pristine.whole_model_ranking = whole_model_ranking_mode(self._policy)
            pristine.ranking_multisample_rna = multisample_rna_support_for_ranking(
                pristine, self._policy,
            )
            pristine.ranking_protein_fit = protein_fit_for_ranking(pristine, self._policy)
            pristine.ranking_directional_protein_fit = directional_protein_fit_for_ranking(
                pristine, self._policy,
            )
            pristine.rna_cds_ranking = bool(self._policy.get("rna_cds_boundary_ranking", False))
            pristine.rna_cds_competition = bool(self._policy.get("rna_cds_competition_ranking", False))
            pristine.rna_cds_conflict = rna_cds_conflict_penalty(row, self._policy)
            if pristine.whole_model_ranking:
                whole_model_rna_support(pristine)
            for attribute, column in boolean_fields.items():
                setattr(pristine, attribute, bool(int(audit[column])))
            pristine.rna_mode = audit["rna_calibration_mode"]
            pristine.repeat_risk = audit["repeat_risk"]
            pristine.repeat_override_reason = audit["repeat_override_reason"]
            pristine.routes = [] if audit["routes"] == "NA" else audit["routes"].split(",")
            pristine.eligible = (
                pristine.hard_qc and bool(pristine.routes) and pristine.repeat_override
            )
            pristine.decision = "ELIGIBLE_GRAPH" if pristine.eligible else "REJECT_EVIDENCE"
            pristine.reasons = [
                "HARD_QC_FAIL" if not pristine.hard_qc
                else "EVIDENCE_PASS" if pristine.eligible
                else "NO_INDEPENDENT_EVIDENCE_ROUTE"
            ]
            if pristine.eligible and pristine.rna_cds_conflict:
                pristine.reasons.append("RNA_CDS_CONFLICT_DOWNRANKED")
            if pristine.eligible and pristine.rna_cds_competition and json.loads(row["rna_cds_countered_conflicts_json"]):
                pristine.reasons.append("RNA_CDS_CONFLICT_COUNTERED_BY_OWN_RNA")
            if pristine.eligible and pristine.rna_cds_ranking and row["rna_cds_status"] in {"SUPPORTED", "MIXED"}:
                pristine.reasons.append("RNA_CDS_COMPLETE_STRUCTURE_SUPPORT")
            add_whole_model_ranking_reasons(pristine)
            self._evidence[source.key] = (_fingerprint(pristine), deepcopy(pristine))
            record = deepcopy(pristine)
            record.locus_id = audit["locus_id"]
            record.selected = bool(int(audit["selected"]))
            record.eligible = bool(int(audit["eligible"]))
            record.active_fusion = bool(int(audit["active_fusion"]))
            record.fusion_center_precedence = bool(int(audit["fusion_center_precedence"]))
            record.decision = audit["decision"]
            record.blockers = [] if audit["blockers"] == "NA" else audit["blockers"].split(",")
            record.reasons = [] if audit["reasons"] == "NA" else audit["reasons"].split(",")
            if record.decision == "TE_SELECTED":
                record.row["te_selected"] = "1"
            self.records[record.key] = record
            self._members.setdefault(record.locus_id, set()).add(record.key)
        self._index = IntervalIndex(
            self.records.values(), seqid=lambda r: r.seqid,
            start=lambda r: r.start, end=lambda r: r.end, strand=lambda r: r.strand,
        )

    def update(self, upserts: Sequence[Record], removed: Sequence[str] = ()) -> UpdateStats:
        new = {record.key: record for record in upserts}
        if len(new) != len(upserts):
            raise SelectionError("duplicate candidate IDs in unified update")
        removed = set(removed)
        if removed & new.keys():
            raise SelectionError("candidate cannot be removed and updated together")
        if removed - self.records.keys():
            raise SelectionError("cannot remove unknown candidate")
        affected = set(new) | removed
        for key in list(affected):
            if key in self.records:
                affected.update(self._members[self.records[key].locus_id])
        if self._index is not None:
            for record in new.values():
                for other in self._index.query(record.seqid, record.start,
                                               record.end, record.strand):
                    if other.key in removed:
                        continue
                    if (set(record.introns) & set(other.introns)
                            or overlap_bp(record.intervals, other.intervals)):
                        affected.update(self._members[other.locus_id])
        evaluated = reused = 0
        working = []
        for key in sorted(affected - removed):
            source = new.get(key, self.records.get(key))
            row = dict(source.row)
            # These are downstream decisions, never evidence-cache inputs.
            row.pop("selection_policy_id", None)
            row.pop("te_selected", None)
            row["te_policy_id"] = self._policy["te_policy_id"]
            fingerprint = _fingerprint(Record(row=row, cds=source.cds))
            cached = self._evidence.get(key)
            if cached is not None and cached[0] == fingerprint:
                record = deepcopy(cached[1])
                reused += 1
            else:
                record = Record(row=row, cds=source.cds)
                evaluate_record(
                    record, self._policy, sample_count=self.sample_count,
                    replicates_per_condition=self.replicates_per_condition,
                    additional_routes=supplemental_routes,
                )
                record.row["selection_policy_id"] = self._policy["policy_id"]
                self._evidence[key] = (fingerprint, deepcopy(record))
                evaluated += 1
            working.append(record)
        loci, edges = build_loci_and_edges(self.assembly_id, working)
        resolve_loci(loci, edges, self._policy)
        self.last_loci, self.last_edges = loci, edges
        for key in affected:
            if key in self.records:
                self._members.pop(self.records[key].locus_id, None)
        for key in removed:
            self.records.pop(key)
            self._evidence.pop(key, None)
        self.records.update((record.key, record) for record in working)
        self._members.update((locus, {r.key for r in members})
                             for locus, members in loci.items())
        if upserts or removed:
            self._index = IntervalIndex(
                self.records.values(), seqid=lambda r: r.seqid,
                start=lambda r: r.start, end=lambda r: r.end,
                strand=lambda r: r.strand,
            )
        self.stats = UpdateStats(evaluated, reused, len(working), len(loci),
                                len(self.records) - len(working))
        return self.stats
