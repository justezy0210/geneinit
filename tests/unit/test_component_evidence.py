import tempfile
from pathlib import Path
import unittest

from geneinit.evidence.components import (
    CandidateRecord,
    EvidenceError,
    JunctionAggregate,
    ProteinChain,
    RepeatFeature,
    RnaTranscript,
    iter_candidate_blocks,
    load_diamond_targets,
    load_repeats,
    repeat_metrics,
    load_protein_chains,
    load_stringtie_transcripts,
    predictor_metrics,
    protein_phase_metrics,
    retained_protein_targets,
    score_candidates,
    score_candidate_blocks,
)
from geneinit.evidence.metrics import IntervalIndex
from geneinit.models.candidate import CdsPart


class ComponentEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def candidate(self):
        return CandidateRecord(
            candidate_id="GCAN-test",
            structure_id="GSTR-test",
            seqid="chr1",
            strand="+",
            cds=(CdsPart(1, 4, 0), CdsPart(11, 18, 2)),
            predictors=("augustus", "genemark_etp"),
            translated_cds_length=12,
            orf_status="complete",
        )

    def chain(self):
        return ProteinChain(
            source="direct",
            seqid="chr1",
            strand="+",
            query_id="PROT1",
            target_id="MSTRG.1",
            cds_parts=(CdsPart(1, 4, 0), CdsPart(11, 18, 2)),
            introns=((5, 10),),
            starts=((1, 3),),
            stops=((16, 18),),
        )

    def test_repeat_union_uses_only_classified_te_and_preserves_total_audit(self):
        repeats = (
            RepeatFeature("chr1", 1, 18, "unknown"),
            RepeatFeature("chr1", 11, 14, "LTR/unknown"),
            RepeatFeature("chr1", 13, 16, "DNA/Helitron"),
            RepeatFeature("chr1", 13, 16, "DNA/Helitron"),
            RepeatFeature("chr1", 1, 18, "rRNA"),
        )
        index = IntervalIndex(repeats, seqid=lambda r: r.seqid, start=lambda r: r.start, end=lambda r: r.end)
        result = repeat_metrics(self.candidate(), index)
        self.assertEqual(result["repeat_scope"], "classified_te_v1")
        self.assertEqual(result["repeat_overlap_bp"], 6)
        self.assertEqual(result["repeat_cds_fraction"], 0.5)
        self.assertEqual(result["all_repeat_overlap_bp"], 12)
        self.assertEqual(result["all_repeat_cds_fraction"], 1.0)
        self.assertNotIn("unknown", __import__("json").loads(result["repeat_class_bp_json"]))
        self.assertIn("unknown", __import__("json").loads(result["all_repeat_class_bp_json"]))

    def test_repeat_classification_does_not_infer_te_from_arbitrary_names(self):
        from geneinit.evidence.repeats import is_classified_te
        path = self.root / "repeats.gff3"
        attributes = (
            "Classification=unknown;Name=LTR", "Classification=rRNA", "Class=Simple_repeat",
            "Classification=Low_complexity", "Classification=SINE?", "Name=TE_123",
            "Classification=lTr%2FGypsy", "class=DNA/unknown", "Name=family#LINE/L1",
        )
        path.write_text("".join(
            f"chr1\tx\trepeat_region\t1\t18\t.\t+\t.\t{attrs}\n"
            for attrs in attributes
        ))
        states = {r.classification: is_classified_te(r.classification) for r in load_repeats(path)}
        self.assertEqual(sum(states.values()), 3)
        for label in ("unknown", "rRNA", "Simple_repeat", "Low_complexity", "SINE?", "Unknown"):
            self.assertFalse(states[label])

    def test_unknown_only_repeats_have_zero_te_risk(self):
        index = IntervalIndex([RepeatFeature("chr1", 1, 18, "Unknown")],
                              seqid=lambda r: r.seqid, start=lambda r: r.start, end=lambda r: r.end)
        result = repeat_metrics(self.candidate(), index)
        self.assertEqual(result["repeat_cds_fraction"], 0)
        self.assertEqual(result["repeat_top_class"], "NA")
        self.assertEqual(result["all_repeat_cds_fraction"], 1)

    def test_phase_metrics_follow_transcript_cumulative_frame(self):
        self.assertEqual(protein_phase_metrics(self.candidate(), self.chain()), (2, 2, 1.0))

    def test_retained_chain_must_exist_in_raw_diamond_pairs(self):
        direct = self.root / "direct.tsv"
        nonhc = self.root / "nonhc.tsv"
        direct.write_text("MSTRG.1\tPROT1\n")
        nonhc.write_text("7_g\tPROT2\n")
        seeds = load_diamond_targets(direct, nonhc)
        retained = retained_protein_targets((self.chain(),), seeds)
        self.assertEqual(retained["PROT1"], frozenset({"direct:MSTRG.1"}))

    def test_nonhc_trace_is_applied_to_chain_coordinates(self):
        gff = self.root / "nonhc.gff"
        gff.write_text(
            "nonhc_1\tProtHint\tCDSpart\t3\t8\t1\t-\t0\tgrp=PROT2_7_g;\n"
        )
        chains = load_protein_chains(
            gff, "nonhc", trace={"nonhc_1": ("chr7", 101, 200)}
        )
        self.assertEqual(chains[0].seqid, "chr7")
        self.assertEqual(chains[0].cds_parts[0], CdsPart(103, 108, 0))

    def test_strandless_stringtie_model_is_not_assigned_to_both_strands(self):
        gtf = self.root / "merged.gtf"
        gtf.write_text(
            'chr1\tStringTie\texon\t1\t10\t.\t.\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n'
            'chr1\tStringTie\texon\t20\t30\t.\t+\t.\tgene_id "MSTRG.2"; transcript_id "MSTRG.2.1";\n'
        )
        transcripts = load_stringtie_transcripts(gtf)
        self.assertEqual([item.transcript_id for item in transcripts], ["MSTRG.2.1"])

    def test_component_score_keeps_independent_evidence_layers(self):
        candidate = self.candidate()
        chain = self.chain()
        transcript = RnaTranscript("MSTRG.1", "chr1", "+", ((1, 4), (11, 18)))
        junction = JunctionAggregate(
            total_hq_reads=6,
            sample_count_hq_ge1=2,
            condition_count_ge2rep_hq_ge1=1,
            condition_count_ge2rep_hq_ge3=1,
        )
        seeds = {"PROT1": frozenset({"direct:MSTRG.1"})}
        rows = score_candidates(
            (candidate,),
            (transcript,),
            (chain,),
            (RepeatFeature("chr1", 11, 12, "LTR"),),
            {("chr1", 5, 10, "+"): junction},
            {"chr1": "ATGAGTNNAGAACCCTAA"},
            seeds,
            seeds,
            1,
        )
        row = rows[0]
        self.assertEqual(row["rna_tier"], 4)
        self.assertEqual(row["protein_tier"], 4)
        self.assertEqual(row["predictor_tier"], 3)
        self.assertEqual(row["protein_multiplicity_status"], "UNIQUE_PROTEIN_HIT")
        self.assertAlmostEqual(row["repeat_cds_fraction"], 2 / 12)
        self.assertEqual(row["major_splice_count"], 1)

    def test_concordance_uses_new_lineage_in_collapsed_partner(self):
        candidate = CandidateRecord(
            "GCAN-a", "GSTR-a", "chr1", "+",
            (CdsPart(1, 30, 0), CdsPart(41, 70, 0)),
            ("augustus",), 60, "complete",
        )
        partner = CandidateRecord(
            "GCAN-b", "GSTR-b", "chr1", "+",
            (CdsPart(1, 30, 0), CdsPart(41, 67, 0)),
            ("augustus", "genemark_etp"), 57, "complete",
        )
        index = IntervalIndex(
            (candidate, partner), seqid=lambda item: item.seqid,
            start=lambda item: item.start, end=lambda item: item.end,
            strand=lambda item: item.strand,
        )
        metrics = predictor_metrics(candidate, index)
        self.assertEqual(metrics["predictor_other_candidate_id"], "GCAN-b")
        self.assertGreaterEqual(metrics["predictor_tier"], 1)

    def test_tied_best_protein_is_conservatively_multi_hit(self):
        candidate = self.candidate()
        direct = self.chain()
        nonhc = ProteinChain(
            source="nonhc", seqid="chr1", strand="+",
            query_id="PROT2", target_id="7_g",
            cds_parts=direct.cds_parts, introns=direct.introns,
            starts=direct.starts, stops=direct.stops,
        )
        retained = {
            "PROT1": frozenset({"direct:MSTRG.1"}),
            "PROT2": frozenset({"nonhc:7_g", "nonhc:8_g"}),
        }
        rows = score_candidates(
            (candidate,), (), (direct, nonhc), (), {},
            {"chr1": "ATGAGTNNAGAACCCTAA"}, retained, retained, 1,
        )
        self.assertEqual(rows[0]["protein_tied_best_chain_count"], 2)
        self.assertEqual(
            rows[0]["protein_multiplicity_status"],
            "AMBIGUOUS_PROTEIN_MULTI_HIT",
        )

    def test_streamed_overlap_blocks_preserve_scoring_results(self):
        candidate = self.candidate()
        distant = CandidateRecord(
            "GCAN-z", "GSTR-z", "chr1", "+",
            (CdsPart(40, 48, 0),), ("augustus",), 9, "complete",
        )
        table = self.root / "candidates.tsv"
        table.write_text(
            "candidate_id\tstructure_id\tseqid\tstrand\tpredictors\t"
            "translated_cds_length\torf_status\tcds_signature\n"
            "GCAN-test\tGSTR-test\tchr1\t+\taugustus,genemark_etp\t12\t"
            "complete\t1-4:0,11-18:2\n"
            "GCAN-z\tGSTR-z\tchr1\t+\taugustus\t9\tcomplete\t40-48:0\n",
            encoding="utf-8",
        )
        blocks = list(iter_candidate_blocks(table))
        self.assertEqual([len(block) for block in blocks], [1, 1])
        expected = score_candidates(
            (candidate, distant), (), (), (), {},
            {"chr1": "ATGAGTNNAGAACCCTAANNNNNNNNNNNNNNNNNNNNNATGAAATAA"},
            {}, {}, 1,
        )
        observed = list(
            score_candidate_blocks(
                blocks, (), (), (), {},
                {"chr1": "ATGAGTNNAGAACCCTAANNNNNNNNNNNNNNNNNNNNNATGAAATAA"},
                {}, {}, 1,
            )
        )
        self.assertEqual(observed, expected)

    def test_many_sparse_candidates_have_bounded_streaming_blocks(self):
        table = self.root / "many_candidates.tsv"
        with table.open("w", encoding="utf-8") as handle:
            handle.write(
                "candidate_id\tstructure_id\tseqid\tstrand\tpredictors\t"
                "translated_cds_length\torf_status\tcds_signature\n"
            )
            for index in range(10_000):
                start = index * 10 + 1
                handle.write(
                    f"GCAN-{index:05d}\tGSTR-{index:05d}\tchr1\t+\taugustus\t"
                    f"3\tcomplete\t{start}-{start + 2}:0\n"
                )
        count = 0
        maximum_block = 0
        for block in iter_candidate_blocks(table):
            count += len(block)
            maximum_block = max(maximum_block, len(block))
        self.assertEqual(count, 10_000)
        self.assertEqual(maximum_block, 1)

    def test_nonadjacent_duplicate_candidate_id_fails_closed(self):
        table = self.root / "duplicate_candidates.tsv"
        table.write_text(
            "candidate_id\tstructure_id\tseqid\tstrand\tpredictors\t"
            "translated_cds_length\torf_status\tcds_signature\n"
            "GCAN-a\tGSTR-a1\tchr1\t+\taugustus\t3\tcomplete\t1-3:0\n"
            "GCAN-b\tGSTR-b\tchr1\t+\taugustus\t3\tcomplete\t11-13:0\n"
            "GCAN-a\tGSTR-a2\tchr1\t+\taugustus\t3\tcomplete\t21-23:0\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(EvidenceError, "duplicate/empty candidate ID"):
            list(iter_candidate_blocks(table))


if __name__ == "__main__":
    unittest.main()
