import unittest

from geneinit.utr.core import Feature, TranscriptModel
from geneinit.utr.grouping import assess_pair


def model(candidate_id, cds, *, exons=None, start_evidence="."):
    features = [
        Feature("chr1", "test", "CDS", start, end, "+", str(phase), {})
        for start, end, phase in cds
    ]
    exon_features = [
        Feature("chr1", "test", "exon", start, end, "+", ".", {})
        for start, end in (exons or [(start, end) for start, end, _phase in cds])
    ]
    return TranscriptModel(
        gene_id=f"g_{candidate_id}", transcript_id=f"t_{candidate_id}",
        candidate_id=candidate_id, seqid="chr1", strand="+", source="test",
        exons=exon_features, cds=features, utr_refine_evidence=start_evidence,
    )


def metadata(*, active=False):
    return {"active_fusion": "1" if active else "0"}


class UtrGroupingTests(unittest.TestCase):
    def test_overlapping_models_from_same_predictor_gene_group(self):
        left = model("left", [(100, 150, 0), (201, 300, 0)])
        right = model("right", [(100, 150, 0), (201, 330, 0)])
        lineages = {
            "left": frozenset({("augustus", "g1")}),
            "right": frozenset({("augustus", "g1")}),
        }
        result = assess_pair(left, right, metadata(), metadata(), lineages)
        self.assertTrue(result.should_group)
        self.assertIn("SAME_PREDICTOR_GENE", result.reasons)

    def test_non_overlapping_models_do_not_group_from_lineage_label_alone(self):
        left = model("left", [(100, 180, 0)])
        right = model("right", [(300, 380, 0)])
        lineages = {
            "left": frozenset({("augustus", "g1")}),
            "right": frozenset({("augustus", "g1")}),
        }
        result = assess_pair(left, right, metadata(), metadata(), lineages)
        self.assertFalse(result.should_group)
        self.assertNotIn("SAME_PREDICTOR_GENE", result.reasons)

    def test_active_fusion_is_a_hard_grouping_boundary(self):
        left = model("left", [(100, 150, 0), (201, 300, 0)])
        right = model("right", [(100, 150, 0), (201, 330, 0)])
        result = assess_pair(
            left, right, metadata(active=True), metadata(),
            {"left": frozenset(), "right": frozenset()},
        )
        self.assertIn("SAME_CDS_INTRON_CHAIN", result.reasons)
        self.assertTrue(result.active_fusion_boundary)
        self.assertFalse(result.should_group)

    def test_same_frame_three_bp_overlap_is_postqc_signal_not_auto_merge(self):
        left = model("left", [(100, 180, 0)])
        right = model("right", [(178, 260, 0)])
        result = assess_pair(
            left, right, metadata(), metadata(),
            {"left": frozenset(), "right": frozenset()},
        )
        self.assertTrue(result.same_frame_cds_overlap_ge3)
        self.assertFalse(result.should_group)

    def test_shared_cds_intron_groups_alternative_models(self):
        left = model("left", [(100, 150, 0), (201, 260, 0), (301, 360, 0)])
        right = model("right", [(90, 150, 0), (201, 270, 0), (321, 390, 0)])
        result = assess_pair(
            left, right, metadata(), metadata(),
            {"left": frozenset(), "right": frozenset()},
        )
        self.assertIn("SHARED_CDS_INTRON", result.reasons)
        self.assertTrue(result.should_group)


if __name__ == "__main__":
    unittest.main()
