import csv
from pathlib import Path
import tempfile
import unittest

from geneinit.models.candidate import CdsPart
from geneinit.selection.core import (
    COMPONENT_REQUIRED_COLUMNS,
    Record,
    SelectionError,
    build_loci_and_edges,
    evaluate_and_select,
    iter_component_blocks,
)


POLICY = {
    "policy_id": "test_policy",
    "hard_qc": {"common_qc_pass": 1, "min_protein_aa": 50, "single_exon_min_protein_aa": 100},
    "rna": {
        "min_recurrent_fraction": 0.8,
        "max_nonrecurrent_introns": 1,
        "strong_min_tier": 3,
        "moderate_min_tier": 2,
        "min_samples_recurrent": 2,
        "min_replicates_strong": 2,
        "long_gene_evidence_rescue": {
            "enabled": True,
            "min_introns": 20,
            "min_recurrent_fraction": 0.8,
            "min_protein_tier": 2,
            "min_protein_chain_recall": 0.9,
        },
    },
    "protein": {
        "moderate_min_tier": 2,
        "strong_min_tier": 3,
        "min_core_metric": 0.8,
        "required_terminal_matches": 2,
        "allow_strong_unique_core_route": True,
    },
    "predictor": {
        "exact_min_tier": 3,
        "require_both_predictors": True,
        "allow_exact_moderate_unique_route": True,
    },
    "repeat": {"high_cds_fraction": 0.5, "severe_cds_fraction": 0.8},
    "graph": {},
}


def make_record(
    key,
    cds,
    *,
    hq=None,
    samples=None,
    replicate=None,
    rna_tier=3,
    protein_tier=0,
    protein_unique=False,
    protein_precision=0.0,
    protein_intron_precision=0.0,
    protein_intron_recall=0.0,
    terminals=0,
    protein_aa=200,
):
    parts = tuple(CdsPart(start, end, phase) for start, end, phase in cds)
    introns = tuple(
        (left.end + 1, right.start - 1)
        for left, right in zip(parts, parts[1:])
    )
    hq = list(hq if hq is not None else [5] * len(introns))
    samples = list(samples if samples is not None else [2] * len(introns))
    replicate = list(replicate if replicate is not None else [1] * len(introns))
    row = {
        "candidate_id": key,
        "structure_id": "S-" + key,
        "seqid": "chr1",
        "strand": "+",
        "start": str(parts[0].start),
        "end": str(parts[-1].end),
        "cds_signature": ",".join(
            f"{part.start}-{part.end}:{part.phase}" for part in parts
        ),
        "introns_json": __import__("json").dumps(introns),
        "intron_count": str(len(introns)),
        "orf_status": "complete",
        "common_qc_pass": "1",
        "translated_protein_aa": str(protein_aa),
        "ambiguous_splice_count": "0",
        "noncanonical_splice_count": "0",
        "splice_motifs": __import__("json").dumps(["GT-AG"] * len(introns)),
        "rna_tier": str(rna_tier),
        "rna_hq_read_counts": __import__("json").dumps(hq),
        "rna_sample_hq_ge1_counts": __import__("json").dumps(samples),
        "rna_condition_ge2rep_hq_ge1_counts": __import__("json").dumps(replicate),
        "protein_tier": str(protein_tier),
        "protein_multiplicity_status": "UNIQUE_PROTEIN_HIT" if protein_unique else "NO_PROTEIN_CHAIN",
        "protein_target_linked": "1" if protein_unique else "0",
        "protein_model_precision": str(protein_precision),
        "protein_chain_recall": str(protein_precision),
        "protein_cds_f1": str(protein_precision),
        "protein_phase_fraction": str(protein_precision),
        "protein_intron_precision": str(protein_intron_precision),
        "protein_intron_recall": str(protein_intron_recall),
        "protein_start_match": str(min(terminals, 1)),
        "protein_stop_match": str(max(0, terminals - 1)),
        "predictor_tier": "0",
        "predictors": "augustus",
        "repeat_cds_fraction": "0",
        "repeat_scope": "classified_te_v1",
    }
    return Record(row=row, cds=parts)


class SelectionCoreTests(unittest.TestCase):
    def run_selection(self, records):
        loci, edges = build_loci_and_edges("ASM", records)
        selected = evaluate_and_select(
            loci, edges, POLICY, sample_count=6, replicates_per_condition=3
        )
        return {record.key for record in selected}, records, edges

    def test_singleton_condition_does_not_disable_candidate_replicated_rna(self):
        from geneinit.selection.core import evaluate_record, audit_row
        for global_min in (1, 3):
            record = make_record("supported", ((100, 150, 0), (201, 300, 0)))
            evaluate_record(record, POLICY, sample_count=25, replicates_per_condition=global_min)
            self.assertTrue(record.strong_rna)
            self.assertEqual(record.rna_mode, "REPLICATED_STRONG")
            self.assertEqual(audit_row(record, "test")["strong_intron_count"], 1)

    def test_other_replicated_conditions_do_not_make_this_candidate_strong(self):
        from geneinit.selection.core import evaluate_record
        record = make_record("unreplicated", ((100, 150, 0), (201, 300, 0)), replicate=(0,))
        evaluate_record(record, POLICY, sample_count=25, replicates_per_condition=3)
        self.assertFalse(record.strong_rna)
        self.assertTrue(record.recurrent_rna)
        self.assertEqual(record.rna_mode, "MULTISAMPLE_RECURRENT")

    def test_classified_te_strong_rna_only_restores_legacy_repeat_rescue(self):
        for fraction in (0.5, 0.8, 1.0):
            record = make_record("rna_only", ((100, 150, 0), (201, 300, 0)))
            record.row["repeat_cds_fraction"] = str(fraction)
            selected, _, _ = self.run_selection([record])
            self.assertTrue(record.strong_rna)
            self.assertEqual(selected, {"rna_only"})
            self.assertTrue(record.repeat_override)
            self.assertEqual(record.repeat_override_reason, "STRONG_RECURRENT_RNA")
            self.assertFalse(record.moderate_unique_protein)

    def test_classified_te_strong_rna_rescue_preserves_protein_audit(self):
        for status, tier, linked in (
            ("NO_PROTEIN_CHAIN", 0, 0),
            ("UNIQUE_PROTEIN_HIT", 1, 1),
            ("AMBIGUOUS_PROTEIN_MULTI_HIT", 1, 1),
            ("AMBIGUOUS_PROTEIN_MULTI_HIT", 2, 0),
            ("UNRESOLVED_PROTEIN_TARGET", 2, 0),
            ("UNIQUE_PROTEIN_HIT", 2, 1),
        ):
            with self.subTest(status=status, tier=tier, linked=linked):
                record = make_record("joint", ((100, 150, 0), (201, 300, 0)),
                                     protein_tier=tier)
                record.row.update(protein_multiplicity_status=status,
                                  protein_target_linked=str(linked), repeat_cds_fraction="0.9")
                selected, _, _ = self.run_selection([record])
                self.assertEqual(selected, {"joint"})
                self.assertEqual(record.repeat_override_reason, "STRONG_RECURRENT_RNA")
                self.assertEqual(record.row["protein_multiplicity_status"], status)
                self.assertEqual(record.row["protein_target_linked"], str(linked))
                self.assertEqual(record.moderate_unique_protein,
                                 status == "UNIQUE_PROTEIN_HIT" and tier >= 2)

    def test_classified_te_strong_rna_rescues_linked_multihit_protein(self):
        from geneinit.selection.core import audit_row
        for fraction in (0.5, 0.8, 1.0):
            with self.subTest(fraction=fraction):
                record = make_record("multihit", ((100, 150, 0), (201, 300, 0)),
                                     protein_tier=2)
                record.row.update(protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
                                  protein_target_linked="1", repeat_cds_fraction=str(fraction))
                selected, _, _ = self.run_selection([record])
                self.assertEqual(selected, {"multihit"})
                self.assertEqual(audit_row(record, "test")["repeat_override_reason"],
                                 "STRONG_RECURRENT_RNA")
                self.assertFalse(record.moderate_unique_protein)
                self.assertFalse(record.strong_unique_protein)
                self.assertFalse(record.strong_unique_full_protein)
                self.assertEqual(record.row["protein_multiplicity_status"],
                                 "AMBIGUOUS_PROTEIN_MULTI_HIT")

    def test_rna_te_rescue_requires_strong_replicated_rna_and_hard_qc(self):
        failures = (
            {"rna_condition_ge2rep_hq_ge1_counts": "[0]"},
            {"rna_hq_read_counts": "[2]"},
            {"rna_tier": "2"},
            {"orf_status": "partial"},
            {"common_qc_pass": "0"},
            {"translated_protein_aa": "49"},
        )
        for fraction in (0.5, 0.8):
            for failure in failures:
                with self.subTest(fraction=fraction, failure=failure):
                    record = make_record("unsupported", ((100, 150, 0), (201, 300, 0)),
                                         protein_tier=2)
                    record.row.update(protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
                                      protein_target_linked="1", repeat_cds_fraction=str(fraction))
                    record.row.update(failure)
                    selected, _, _ = self.run_selection([record])
                    self.assertFalse(selected)

    def test_multihit_te_rescue_preserves_single_exon_and_graph_restrictions(self):
        single = make_record("single", ((100, 400, 0),),
                             protein_tier=4, protein_precision=1, terminals=2)
        single.row.update(protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
                          protein_target_linked="1", repeat_cds_fraction="0.9")
        selected, _, _ = self.run_selection([single])
        self.assertFalse(selected)

        records = [make_record(key, ((100, 150, 0), (201, 300, 0)), protein_tier=2)
                   for key in ("first", "second")]
        for record in records:
            record.row.update(protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
                              protein_target_linked="1", repeat_cds_fraction="0.9")
        selected, records, _ = self.run_selection(records)
        self.assertEqual(len(selected), 1)
        self.assertTrue(all(record.repeat_override for record in records))
        self.assertEqual(sum(record.decision == "REJECT_INCOMPATIBLE_GRAPH" for record in records), 1)

    def test_severe_te_with_recurrent_but_not_strong_rna_stays_rejected(self):
        record = make_record("recurrent", ((100, 150, 0), (201, 300, 0)),
                             protein_tier=2, protein_unique=True, replicate=(0,))
        record.row["repeat_cds_fraction"] = "0.8"
        selected, _, _ = self.run_selection([record])
        self.assertFalse(selected)

    def test_high_te_recurrent_unique_protein_rescue_is_preserved(self):
        record = make_record("recurrent", ((100, 150, 0), (201, 300, 0)),
                             protein_tier=2, protein_unique=True, replicate=(0,))
        record.row["repeat_cds_fraction"] = "0.5"
        selected, _, _ = self.run_selection([record])
        self.assertEqual(selected, {"recurrent"})
        self.assertFalse(record.strong_rna)
        self.assertEqual(record.repeat_override_reason, "RECURRENT_RNA_MODERATE_UNIQUE_PROTEIN")

    def test_full_unique_protein_repeat_rescue_is_preserved(self):
        record = make_record("full", ((100, 400, 0),), rna_tier=0,
                             protein_tier=3, protein_unique=True, protein_precision=1, terminals=2)
        record.row["repeat_cds_fraction"] = "1"
        selected, _, _ = self.run_selection([record])
        self.assertEqual(selected, {"full"})
        self.assertEqual(record.repeat_override_reason, "STRONG_UNIQUE_FULL_PROTEIN")

    def test_legacy_repeat_scope_is_rejected(self):
        record = make_record("legacy", ((100, 150, 0), (201, 300, 0)))
        for scope in (None, "all_repeats", ""):
            row = dict(record.row)
            if scope is None:
                del row["repeat_scope"]
            else:
                row["repeat_scope"] = scope
            with self.assertRaisesRegex(SelectionError, "rerun score_component_evidence"):
                Record(row, record.cds)

    def test_strong_rna_fraction_and_missing_intron_limits_are_preserved(self):
        from geneinit.selection.core import evaluate_record
        for n, supported, expected in ((5, 4, True), (6, 4, False), (10, 8, False)):
            record = make_record("fraction", tuple((100*i+1, 100*i+50, 0) for i in range(n+1)),
                                 replicate=[1]*supported+[0]*(n-supported))
            evaluate_record(record, POLICY, sample_count=25, replicates_per_condition=1)
            self.assertEqual(record.strong_rna, expected)

    def test_same_intron_chain_keeps_one_representative(self):
        left = make_record("left", ((100, 150, 0), (201, 300, 0)), protein_aa=180)
        right = make_record("right", ((90, 150, 0), (201, 320, 0)), protein_aa=170)
        selected, records, edges = self.run_selection([left, right])
        self.assertEqual(len(selected), 1)
        self.assertEqual(edges[0].relation, "SAME_CHAIN")
        self.assertEqual(sum(record.decision == "REJECT_INCOMPATIBLE_GRAPH" for record in records), 1)

    def test_alt_isoform_requires_unique_recurrent_intron(self):
        anchor = make_record("anchor", ((100, 150, 0), (201, 400, 0)), protein_aa=220)
        alternative = make_record(
            "alternative",
            ((100, 150, 0), (201, 250, 0), (301, 400, 0)),
            hq=(5, 5), samples=(2, 2), replicate=(1, 1), protein_aa=180,
        )
        selected, _records, edges = self.run_selection([anchor, alternative])
        self.assertEqual(selected, {"anchor", "alternative"})
        self.assertIn("ALT_SPLICE", {edge.relation for edge in edges})

    def test_alt_isoform_without_unique_evidence_is_rejected(self):
        anchor = make_record("anchor", ((100, 150, 0), (201, 400, 0)), protein_aa=220)
        alternative = make_record(
            "alternative",
            ((100, 150, 0), (201, 250, 0), (301, 400, 0)),
            hq=(0, 0), samples=(0, 0), replicate=(0, 0), rna_tier=0,
            protein_tier=3, protein_unique=True, protein_precision=0.8,
            protein_intron_precision=0.8, protein_intron_recall=0.8,
            protein_aa=180,
        )
        selected, records, _edges = self.run_selection([anchor, alternative])
        self.assertEqual(selected, {"anchor"})
        self.assertEqual(alternative.decision, "REJECT_ALT_ISOFORM_NO_UNIQUE_EVIDENCE")

    def test_recurrent_rna_route_is_disabled_for_one_sample(self):
        record = make_record(
            "weak_design", ((100, 150, 0), (201, 300, 0)),
            hq=(2,), samples=(1,), replicate=(0,), rna_tier=2,
        )
        loci, edges = build_loci_and_edges("ASM", [record])
        selected = evaluate_and_select(
            loci, edges, POLICY, sample_count=1, replicates_per_condition=1
        )
        self.assertEqual(selected, [])
        self.assertEqual(record.rna_mode, "DISABLED_WEAK_DESIGN")

    def test_fusion_center_requires_strong_bridge_or_full_protein(self):
        left = make_record("left", ((100, 250, 0),), protein_tier=4, protein_unique=True, protein_precision=1.0, terminals=2)
        right = make_record("right", ((401, 550, 0),), protein_tier=4, protein_unique=True, protein_precision=1.0, terminals=2)
        center = make_record(
            "center", ((100, 300, 0), (351, 550, 0)),
            hq=(2,), samples=(2,), replicate=(0,), rna_tier=2, protein_aa=250,
        )
        selected, _records, _edges = self.run_selection([left, center, right])
        self.assertTrue(center.active_fusion)
        self.assertNotIn("center", selected)
        self.assertEqual(center.decision, "REJECT_ACTIVE_FUSION_CENTER")

    def test_streamed_blocks_are_exact_selection_graph_partitions(self):
        source = [
            make_record("a", ((100, 180, 0),)),
            make_record("b", ((151, 230, 0),)),
            make_record("c", ((400, 500, 0),)),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "components.tsv"
            columns = sorted(COMPONENT_REQUIRED_COLUMNS)
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=columns, delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for record in source:
                    writer.writerow(
                        {column: record.row.get(column, "NA") for column in columns}
                    )
            blocks = list(iter_component_blocks(path))
        self.assertEqual([len(block) for block in blocks], [2, 1])
        all_records = [record for block in blocks for record in block]
        all_loci, all_edges = build_loci_and_edges("ASM", all_records)
        block_locus_count = 0
        block_edges = []
        for block in blocks:
            loci, edges = build_loci_and_edges("ASM", block)
            block_locus_count += len(loci)
            block_edges.extend(edges)
        self.assertEqual(block_locus_count, len(all_loci))
        self.assertEqual(
            [(edge.left, edge.right, edge.relation) for edge in block_edges],
            [(edge.left, edge.right, edge.relation) for edge in all_edges],
        )

    def test_nonadjacent_duplicate_candidate_id_fails_closed(self):
        source = [
            make_record("a", ((100, 102, 0),)),
            make_record("b", ((200, 202, 0),)),
            make_record("a", ((300, 302, 0),)),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "components.tsv"
            columns = sorted(COMPONENT_REQUIRED_COLUMNS)
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=columns, delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for record in source:
                    writer.writerow(
                        {column: record.row.get(column, "NA") for column in columns}
                    )
            with self.assertRaisesRegex(
                SelectionError, "duplicate/empty candidate"
            ):
                list(iter_component_blocks(path))


if __name__ == "__main__":
    unittest.main()
