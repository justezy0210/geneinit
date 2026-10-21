import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from geneinit.postqc.core import (
    Node,
    build_edges,
    coding_components,
    derive_rna_bridge_thresholds,
    edge_audit_rows,
    iterate_structural_qc,
    PostQcError,
    recover_supported_bridges,
    recover_strong_rna_bridges,
    retire_unsupported_connectors,
    run_postqc,
    same_frame_overlap_bp,
)
from geneinit.evidence.junctions import AGGREGATE_COLUMNS
from geneinit.evidence.components import JunctionAggregate
from geneinit.utr.core import Feature, TranscriptModel


def transcript(transcript_id, candidate_id, start, end, phase=0, gene_id="g1"):
    cds = Feature("chr1", "test", "CDS", start, end, "+", str(phase), {})
    exon = Feature("chr1", "test", "exon", start, end, "+", ".", {})
    return TranscriptModel(
        gene_id=gene_id, transcript_id=transcript_id, candidate_id=candidate_id,
        seqid="chr1", strand="+", source="test", exons=[exon], cds=[cds],
    )


def node(transcript_id, candidate_id, start, end, *, external=True, gene_id="g1"):
    return Node(
        model=transcript(transcript_id, candidate_id, start, end, gene_id=gene_id),
        input_gene_id=gene_id, source_status="RELEASED", selection_locus_id="L1",
        routes="RNA" if external else "FROZEN_BALANCED_BASE",
        external_support=external, rna_tier=3 if external else 0, protein_tier=0,
    )


class PostQcCoreTests(unittest.TestCase):
    def test_structural_qc_iterates_until_a_no_change_round(self):
        nodes = {"t1": node("t1", "C1", 100, 180)}
        calls = 0

        def apply_iteration():
            nonlocal calls
            calls += 1
            if calls == 1:
                nodes["t1"].retired = True
                return [], [], [], [{"transcript_id": "t1"}]
            return [], [], [], []

        _promotions, _reviews, retirements, rows, converged = iterate_structural_qc(
            nodes, 5, apply_iteration
        )
        self.assertEqual(converged, 2)
        self.assertEqual([row["changed"] for row in rows], [1, 0])
        self.assertEqual(len(retirements), 1)

    def test_structural_qc_fails_when_maximum_iterations_do_not_converge(self):
        nodes = {"t1": node("t1", "C1", 100, 180)}

        def apply_iteration():
            nodes["t1"].retired = not nodes["t1"].retired
            return [], [], [], []

        with self.assertRaisesRegex(PostQcError, "maximum_iterations=2"):
            iterate_structural_qc(nodes, 2, apply_iteration)

    def test_adaptive_intron_thresholds_use_quantiles_and_fallback(self):
        policy = {
            "mode": "adaptive",
            "minimum_training_intron_count": 3,
            "automatic_quantile": 0.5,
            "review_quantile": 1.0,
            "automatic_floor_bp": 100,
            "automatic_ceiling_bp": 50000,
            "review_multiplier": 4.0,
            "review_ceiling_bp": 500000,
            "fallback_automatic_max_intron_bp": 5000,
            "fallback_review_max_intron_bp": 20000,
        }
        adaptive = derive_rna_bridge_thresholds([100, 200, 6000], policy)
        self.assertEqual(adaptive.effective_mode, "adaptive")
        self.assertEqual(adaptive.automatic_max_intron_bp, 200)
        self.assertEqual(adaptive.review_max_intron_bp, 6000)
        fallback = derive_rna_bridge_thresholds([100, 200], policy)
        self.assertEqual(fallback.effective_mode, "adaptive_fallback")
        self.assertEqual(fallback.automatic_max_intron_bp, 5000)
        self.assertEqual(fallback.review_max_intron_bp, 20000)
        fixed = derive_rna_bridge_thresholds(
            [100, 200],
            {
                "mode": "fixed",
                "fixed_automatic_max_intron_bp": 5000,
                "fixed_review_max_intron_bp": 20000,
            },
        )
        self.assertEqual(fixed.effective_mode, "fixed")
        self.assertEqual(fixed.automatic_max_intron_bp, 5000)

    def test_same_frame_overlap_uses_local_cds_phase(self):
        left = transcript("t1", "C1", 100, 180, phase=0)
        same = transcript("t2", "C2", 151, 230, phase=0)
        different = transcript("t3", "C3", 152, 230, phase=0)
        self.assertEqual(same_frame_overlap_bp(left, same), 30)
        self.assertEqual(same_frame_overlap_bp(left, different), 0)

    def test_components_ignore_old_gene_boundaries(self):
        nodes = {
            "t1": node("t1", "C1", 100, 180, gene_id="g1"),
            "t2": node("t2", "C2", 151, 230, gene_id="g2"),
        }
        components, _roots = coding_components(nodes, build_edges(nodes))
        self.assertEqual(len(components), 1)

    def test_active_fusion_is_preserved_as_a_hard_postqc_boundary(self):
        nodes = {
            "t1": node("t1", "C1", 100, 180, gene_id="g1"),
            "t2": node("t2", "C2", 151, 230, gene_id="g2"),
        }
        nodes["t1"].active_fusion = True
        edges = build_edges(nodes)
        self.assertEqual(len(edges), 1)
        self.assertTrue(edges[0].active_fusion_boundary)
        self.assertFalse(edges[0].linked(3))
        components, _roots = coding_components(nodes, edges)
        self.assertEqual(len(components), 2)

    def test_weak_articulation_between_supported_components_is_retired(self):
        nodes = {
            "left": node("left", "C1", 100, 180, external=True, gene_id="g1"),
            "connector": node("connector", "C2", 151, 250, external=False, gene_id="g1"),
            "right": node("right", "C3", 220, 300, external=True, gene_id="g1"),
        }
        rows = retire_unsupported_connectors(nodes)
        self.assertTrue(nodes["connector"].retired)
        self.assertEqual(rows[0]["status"], "RETIRED_UNSUPPORTED_FUSION_CONNECTOR")
        components, _roots = coding_components(nodes, build_edges(nodes))
        self.assertEqual(len(components), 2)

    def test_active_fusion_connector_is_never_retired(self):
        nodes = {
            "left": node("left", "C1", 100, 180, external=True),
            "connector": node("connector", "C2", 151, 250, external=False),
            "right": node("right", "C3", 220, 300, external=True),
        }
        nodes["connector"].active_fusion = True
        self.assertEqual(retire_unsupported_connectors(nodes), [])
        self.assertFalse(nodes["connector"].retired)

    def test_connector_round_does_not_cascade_over_retirements(self):
        nodes = {
            "left": node("left", "C1", 100, 180, external=True),
            "x": node("x", "C2", 151, 231, external=False),
            "y": node("y", "C3", 202, 282, external=False),
            "right": node("right", "C4", 253, 333, external=True),
        }
        rows = retire_unsupported_connectors(nodes)
        self.assertEqual([row["transcript_id"] for row in rows], ["x"])
        self.assertTrue(nodes["x"].retired)
        self.assertFalse(nodes["y"].retired)

    def test_supported_unselected_candidate_bridging_two_components_is_promoted(self):
        nodes = {
            "left": node("left", "C1", 100, 180, gene_id="g1"),
            "right": node("right", "C2", 220, 300, gene_id="g2"),
        }
        candidate = {
            "candidate_id": "C3", "seqid": "chr1", "strand": "+",
            "start": "151", "end": "249", "cds_signature": "151-249:0",
        }
        evidence = {
            "C3": {
                "candidate_id": "C3", "orf_status": "complete", "common_qc_pass": "1",
                "rna_tier": "2", "protein_tier": "0", "repeat_cds_fraction": "0",
            }
        }
        audit = {
            "C3": {
                "candidate_id": "C3", "locus_id": "L3", "selected": "0",
                "hard_qc_pass": "1", "active_fusion": "0",
            }
        }
        rows = recover_supported_bridges(nodes, [candidate], evidence, audit)
        self.assertEqual(len(rows), 1)
        self.assertIn("GIPROM-C3", nodes)
        components, _roots = coding_components(nodes, build_edges(nodes))
        self.assertEqual(len(components), 1)

    def test_promoted_bridge_is_not_an_anchor_for_later_candidates(self):
        nodes = {
            "left": node("left", "C1", 100, 180, gene_id="g1"),
            "right": node("right", "C2", 220, 300, gene_id="g2"),
        }
        candidates = [
            {
                "candidate_id": "C3", "seqid": "chr1", "strand": "+",
                "start": "151", "end": "249", "cds_signature": "151-249:0",
            },
            {
                "candidate_id": "C4", "seqid": "chr1", "strand": "+",
                "start": "160", "end": "219", "cds_signature": "160-219:0",
            },
        ]
        evidence = {
            candidate_id: {
                "candidate_id": candidate_id, "orf_status": "complete",
                "common_qc_pass": "1", "rna_tier": "2", "protein_tier": "0",
                "repeat_cds_fraction": "0",
            }
            for candidate_id in ("C3", "C4")
        }
        audit = {
            candidate_id: {
                "candidate_id": candidate_id, "locus_id": candidate_id,
                "selected": "0", "hard_qc_pass": "1", "active_fusion": "0",
            }
            for candidate_id in ("C3", "C4")
        }
        rows = recover_supported_bridges(nodes, candidates, evidence, audit)
        self.assertEqual([row["candidate_id"] for row in rows], ["C3"])
        self.assertNotIn("GIPROM-C4", nodes)

    def test_promotion_created_final_edges_are_in_the_audit(self):
        nodes = {
            "left": node("left", "C1", 100, 180, gene_id="g1"),
            "right": node("right", "C2", 220, 300, gene_id="g2"),
        }
        initial_edges = build_edges(nodes)
        promoted = node("promoted", "C3", 151, 249, gene_id=".")
        promoted.source_status = "PROMOTED_SUPPORTED_BRIDGE"
        nodes["promoted"] = promoted
        for value in nodes.values():
            value.postqc_gene_id = "GPQC-1"
        rows = edge_audit_rows(nodes, initial_edges, build_edges(nodes), 3)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["present_initial"] == 0 for row in rows))
        self.assertTrue(all(row["present_final"] == 1 for row in rows))
        self.assertTrue(all(row["introduced_by_promotion"] == 1 for row in rows))

    def test_strong_exact_junction_reconstructs_only_complete_orf_bridge(self):
        nodes = {
            "left": node("left", "C1", 100, 150, gene_id="g1"),
            "right": node("right", "C2", 201, 251, gene_id="g2"),
        }
        sequence = list("C" * 300)
        sequence[99:102] = "ATG"
        sequence[248:251] = "TAA"
        rows = recover_strong_rna_bridges(
            nodes,
            {
                ("chr1", 151, 200, "+"): JunctionAggregate(
                    total_hq_reads=5,
                    sample_count_hq_ge1=2,
                    condition_count_ge2rep_hq_ge1=1,
                    condition_count_ge2rep_hq_ge3=1,
                )
            },
            {"chr1": "".join(sequence)},
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["initial_component_count"], 2)
        self.assertEqual(rows[0]["hit_transcript_count"], 2)
        self.assertEqual(rows[0]["hit_transcripts"], "left,right")
        self.assertTrue(any(value.startswith("GIRNA-") for value in nodes))
        components, _roots = coding_components(nodes, build_edges(nodes))
        self.assertEqual(len(components), 1)

    def test_strong_rna_microexon_candidate_must_link_two_actual_roots(self):
        left = node("left", "C1", 100, 101, gene_id="g1")
        right = Node(
            model=transcript("right", "C2", 151, 202, phase=1, gene_id="g2"),
            input_gene_id="g2", source_status="RELEASED",
            selection_locus_id="L2", routes="RNA", external_support=True,
            rna_tier=3, protein_tier=0,
        )
        nodes = {"left": left, "right": right}
        sequence = list("C" * 250)
        sequence[99:101] = "AT"
        sequence[150] = "G"
        sequence[199:202] = "TAA"
        rows = recover_strong_rna_bridges(
            nodes,
            {
                ("chr1", 102, 150, "+"): JunctionAggregate(
                    total_hq_reads=5,
                    sample_count_hq_ge1=2,
                    condition_count_ge2rep_hq_ge1=1,
                    condition_count_ge2rep_hq_ge3=1,
                )
            },
            {"chr1": "".join(sequence)},
        )
        self.assertEqual(rows, [])
        self.assertFalse(any(value.startswith("GIRNA-") for value in nodes))

    def test_strong_rna_bridge_cannot_use_active_fusion_anchor(self):
        nodes = {
            "left": node("left", "C1", 100, 150, gene_id="g1"),
            "right": node("right", "C2", 201, 251, gene_id="g2"),
        }
        nodes["left"].active_fusion = True
        sequence = list("C" * 300)
        sequence[99:102] = "ATG"
        sequence[248:251] = "TAA"
        rows = recover_strong_rna_bridges(
            nodes,
            {
                ("chr1", 151, 200, "+"): JunctionAggregate(
                    total_hq_reads=5,
                    sample_count_hq_ge1=2,
                    condition_count_ge2rep_hq_ge1=1,
                    condition_count_ge2rep_hq_ge3=1,
                )
            },
            {"chr1": "".join(sequence)},
        )
        self.assertEqual(rows, [])
        self.assertFalse(any(value.startswith("GIRNA-") for value in nodes))

    def test_strong_rna_bridge_cannot_span_active_fusion_model(self):
        nodes = {
            "left": node("left", "C1", 100, 150, gene_id="g1"),
            "barrier": node("barrier", "CF", 160, 190, gene_id="gf"),
            "right": node("right", "C2", 201, 251, gene_id="g2"),
        }
        nodes["barrier"].active_fusion = True
        sequence = list("C" * 300)
        sequence[99:102] = "ATG"
        sequence[248:251] = "TAA"
        rows = recover_strong_rna_bridges(
            nodes,
            {
                ("chr1", 151, 200, "+"): JunctionAggregate(
                    total_hq_reads=5,
                    sample_count_hq_ge1=2,
                    condition_count_ge2rep_hq_ge1=1,
                    condition_count_ge2rep_hq_ge3=1,
                )
            },
            {"chr1": "".join(sequence)},
        )
        self.assertEqual(rows, [])
        self.assertFalse(any(value.startswith("GIRNA-") for value in nodes))

    def test_long_junction_does_not_create_synthetic_rna_bridge(self):
        nodes = {
            "left": node("left", "C1", 100, 150, gene_id="g1"),
            "right": node("right", "C2", 5201, 5251, gene_id="g2"),
        }
        sequence = list("C" * 6000)
        sequence[99:102] = "ATG"
        sequence[5248:5251] = "TAA"
        review_rows = []
        rows = recover_strong_rna_bridges(
            nodes,
            {
                ("chr1", 151, 5200, "+"): JunctionAggregate(
                    total_hq_reads=10,
                    sample_count_hq_ge1=2,
                    condition_count_ge2rep_hq_ge1=1,
                    condition_count_ge2rep_hq_ge3=1,
                )
            },
            {"chr1": "".join(sequence)}, review_rows=review_rows,
        )
        self.assertEqual(rows, [])
        self.assertFalse(any(value.startswith("GIRNA-") for value in nodes))
        self.assertEqual(len(review_rows), 1)
        self.assertEqual(review_rows[0]["threshold_class"], "ADDITIONAL_EVIDENCE")

    def test_end_to_end_cross_gene_link_becomes_one_fixed_point_gene(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            grouped = root / "grouped.gff3"
            grouped.write_text(
                "##gff-version 3\n"
                "chr1\ttest\tgene\t100\t180\t.\t+\t.\tID=G1;\n"
                "chr1\ttest\tmRNA\t100\t180\t.\t+\t.\tID=t1;Parent=G1;candidate_id=C1;\n"
                "chr1\ttest\texon\t100\t180\t.\t+\t.\tID=e1;Parent=t1;\n"
                "chr1\ttest\tCDS\t100\t180\t.\t+\t0\tID=c1;Parent=t1;\n"
                "chr1\ttest\tgene\t151\t230\t.\t+\t.\tID=G2;\n"
                "chr1\ttest\tmRNA\t151\t230\t.\t+\t.\tID=t2;Parent=G2;candidate_id=C2;\n"
                "chr1\ttest\texon\t151\t230\t.\t+\t.\tID=e2;Parent=t2;\n"
                "chr1\ttest\tCDS\t151\t230\t.\t+\t0\tID=c2;Parent=t2;\n",
                encoding="utf-8",
            )

            def write_table(name, columns, rows):
                path = root / name
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(columns)
                    writer.writerows(rows)
                return path

            grouping = write_table(
                "grouping.tsv",
                (
                    "candidate_id", "temporary_transcript_id",
                    "provisional_gene_id", "selection_locus_id", "routes",
                    "active_fusion",
                ),
                (("C1", "t1", "G1", "L1", "STRONG_RECURRENT_RNA", 0),
                 ("C2", "t2", "G2", "L2", "STRONG_UNIQUE_FULL_PROTEIN", 0)),
            )
            candidates = write_table(
                "candidates.tsv",
                ("candidate_id", "seqid", "strand", "start", "end", "cds_signature"),
                (("C1", "chr1", "+", 100, 180, "100-180:0"),
                 ("C2", "chr1", "+", 151, 230, "151-230:0")),
            )
            evidence = write_table(
                "evidence.tsv",
                ("candidate_id", "orf_status", "common_qc_pass", "rna_tier", "protein_tier", "repeat_cds_fraction", "repeat_scope"),
                (("C1", "complete", 1, 3, 0, 0, "classified_te_v1"), ("C2", "complete", 1, 0, 3, 0, "classified_te_v1")),
            )
            audit = write_table(
                "audit.tsv",
                ("candidate_id", "locus_id", "selected", "hard_qc_pass", "active_fusion", "routes"),
                (("C1", "L1", 1, 1, 0, "STRONG_RECURRENT_RNA"),
                 ("C2", "L2", 1, 1, 0, "STRONG_UNIQUE_FULL_PROTEIN")),
            )
            output = root / "postqc.gff3"
            transcript_decisions = root / "transcripts.tsv"
            gene_decisions = root / "genes.tsv"
            edges = root / "edges.tsv"
            promotions = root / "promotions.tsv"
            bridge_reviews = root / "bridge_reviews.tsv"
            retirements = root / "retirements.tsv"
            summary = root / "summary.tsv"
            complete = root / "complete.tsv"
            iterations = root / "iterations.tsv"
            policy = root / "policy.json"
            policy.write_text(
                json.dumps(
                    {
                        "policy_id": "fixture",
                        "postqc_schema": "geneinit_postqc_policy_v2",
                        "postqc": {
                            "minimum_same_frame_overlap_bp": 3,
                            "bridge_recovery": {
                                "minimum_rna_tier": 2,
                                "minimum_protein_tier": 2,
                                "maximum_repeat_cds_fraction": 0.8,
                                "exclude_active_fusion": True,
                            },
                            "rna_bridge_reconstruction": {
                                "minimum_hq_reads": 3,
                                "minimum_replicated_conditions": 1,
                                "intron_length_policy": {
                                    "mode": "adaptive",
                                    "minimum_training_intron_count": 2000,
                                    "automatic_quantile": 0.995,
                                    "review_quantile": 0.999,
                                    "automatic_floor_bp": 5000,
                                    "automatic_ceiling_bp": 50000,
                                    "review_multiplier": 4.0,
                                    "review_ceiling_bp": 500000,
                                    "fallback_automatic_max_intron_bp": 5000,
                                    "fallback_review_max_intron_bp": 20000,
                                },
                                "require_exact_exon_and_cds_boundaries": True,
                                "require_complete_orf": True,
                            },
                            "unsupported_connector": {
                                "minimum_protected_neighbor_components": 2,
                                "external_support_basis": "fixture",
                            },
                            "maximum_iterations": 5,
                        },
                    }
                ),
                encoding="utf-8",
            )
            with gzip.open(root / "junctions.tsv.gz", "wt", encoding="utf-8") as handle:
                handle.write("#" + "\t".join(AGGREGATE_COLUMNS) + "\n")
            (root / "genome.fa").write_text(">chr1\n" + "A" * 500 + "\n", encoding="utf-8")
            run_postqc(
                grouped, grouping, candidates, evidence, audit,
                root / "junctions.tsv.gz", root / "genome.fa", policy,
                ["chr1"], output,
                transcript_decisions, gene_decisions, edges, promotions,
                bridge_reviews, retirements, summary, complete, iterations,
            )
            gene_lines = [line for line in output.read_text().splitlines() if "\tgene\t" in line]
            self.assertEqual(len(gene_lines), 1)
            self.assertIn("CONSOLIDATE_GENE_IDS_AS_ISOFORMS", output.read_text())
            self.assertIn("fixed_point_status\tPASS", complete.read_text())
            self.assertIn("iterations_executed\t1", complete.read_text())
            self.assertTrue(iterations.read_text().rstrip().endswith("\t0\t1"))
            self.assertIn("rna_bridge_threshold_mode_effective\tadaptive_fallback", complete.read_text())


if __name__ == "__main__":
    unittest.main()
