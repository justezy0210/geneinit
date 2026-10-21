"""Score new structural proposals and pass them through the shared selector."""

from __future__ import annotations

from dataclasses import asdict
from typing import Mapping

from geneinit.evidence.components import (
    CandidateRecord, EvidenceContext, format_value,
)
from geneinit.evidence.te_classification import TE_COLUMNS
from geneinit.models.candidate import CdsPart, audit_orf
from geneinit.selection.core import Record, _parse_cds, audit_row
from geneinit.selection.unified import UnifiedSelector


class RepairEvidence(EvidenceContext):
    """Reuse existing alignments and indices; never inherit parental scores."""

    def score(self, model) -> Record:
        parts = tuple(CdsPart(f.start, f.end, int(f.phase))
                      for f in sorted(model.cds, key=lambda f: f.start))
        orf = audit_orf(model.seqid, model.strand, parts, self.genome)
        candidate = CandidateRecord(
            model.candidate_id, model.candidate_id, model.seqid, model.strand,
            parts, (), orf.translated_cds_length, orf.status,
        )
        row = self.score_candidate(candidate)
        return Record({key: format_value(value) for key, value in row.items()}, parts)


class UnifiedRepair:
    def __init__(self, *, policy, assembly_id, sample_count, replicates_per_condition,
                 evidence_by_candidate, audit_by_candidate, scorer):
        self.selector = UnifiedSelector(
            assembly_id, policy, sample_count=sample_count,
            replicates_per_condition=replicates_per_condition,
        )
        records = []
        for key, evidence in evidence_by_candidate.items():
            row = dict(evidence)
            row.update({name: audit_by_candidate[key].get(name, "")
                        for name in TE_COLUMNS if name != "te_selected"})
            records.append(Record(row, _parse_cds(row["cds_signature"])))
        self.selector.restore(records, audit_by_candidate)
        self.scorer = scorer
        self.models = {}
        self.seen = set(self.selector.records)
        self.rounds = []
        self.decisions = []
        self.policy_id = policy["policy_id"]

    def apply(self, nodes, proposals, promotion_rows, generated_classifications):
        """Adopt only selected proposals; re-evaluate entire affected loci."""
        from geneinit.postqc.core import Node, _candidate_model

        updates = []
        for row in promotion_rows:
            candidate_id = row["candidate_id"]
            if candidate_id in self.seen:
                continue
            model = proposals[row["transcript_id"]].model
            record = self.scorer.score(model)
            classification = generated_classifications.get(candidate_id)
            if classification is None:
                raise ValueError("unified repair requires sequence-bound TE classification")
            record.row.update({key: format_value(value)
                               for key, value in classification.items()
                               if key in TE_COLUMNS and key != "te_selected"})
            self.models[candidate_id] = model
            self.seen.add(candidate_id)
            updates.append(record)
        if not updates:
            return [], [], []
        before = set()
        by_candidate = {}
        for node in nodes.values():
            by_candidate[node.candidate_id] = node
            if not node.retired:
                before.add(node.candidate_id)
        stats = self.selector.update(updates)
        round_number = len(self.rounds) + 1
        self.rounds.append({"round": round_number, **asdict(stats)})
        changed_keys = {
            record.key for members in self.selector.last_loci.values()
            for record in members
        }
        for key in sorted(changed_keys):
            self.decisions.append({
                "round": round_number, **audit_row(self.selector.records[key], self.policy_id)
            })
        chosen = self.selector.selected_keys
        retirements = []
        for key in sorted(before - chosen):
            node = by_candidate[key]
            node.retired = True
            node.operation = "RETIRE_UNIFIED_CONFLICT"
            node.reason = self.selector.records[key].decision
            retirements.append({
                "transcript_id": node.transcript_id, "candidate_id": key,
                "input_gene_id": node.input_gene_id, "seqid": node.model.seqid,
                "strand": node.model.strand, "start": node.start, "end": node.end,
                "neighbor_component_count": "NA", "protected_neighbor_component_count": "NA",
                "status": "RETIRED_UNIFIED_SELECTION",
                "reason": node.reason,
            })
        promotions = []
        new_rna = []
        for key in sorted(chosen - before):
            record = self.selector.records[key]
            synthetic = key in self.models
            if key in by_candidate:
                node = by_candidate[key]
                node.retired = False
                node.operation = ""
                node.reason = ""
            else:
                model = self.models[key] if synthetic else _candidate_model(record.row)
                node = Node(
                    model=model, input_gene_id=".",
                    source_status=("PROMOTED_STRONG_RNA_BRIDGE" if synthetic
                                   else "PROMOTED_UNIFIED_RESELECTION"),
                    selection_locus_id=record.locus_id, routes=",".join(record.routes),
                    external_support=True, rna_tier=int(record.row["rna_tier"]),
                    protein_tier=int(record.row["protein_tier"]),
                    active_fusion=record.active_fusion,
                )
                nodes[model.transcript_id] = node
                by_candidate[key] = node
            result = {
                "candidate_id": key, "transcript_id": node.transcript_id,
                "seqid": record.seqid, "strand": record.strand,
                "start": record.start, "end": record.end,
                "initial_component_count": "NA", "hit_transcript_count": "NA",
                "hit_transcripts": "NA", "rna_tier": node.rna_tier,
                "protein_tier": node.protein_tier,
                "status": "PROMOTED_UNIFIED_RNA" if synthetic else "RESTORED_UNIFIED_SELECTION",
            }
            (new_rna if synthetic else promotions).append(result)
        # A new competitor can also change a retained model's fusion status.
        for key in changed_keys:
            node = by_candidate.get(key)
            if node is not None and not node.retired:
                record = self.selector.records[key]
                node.active_fusion = record.active_fusion
                node.selection_locus_id = record.locus_id
                node.routes = ",".join(record.routes)
        return promotions, new_rna, retirements
