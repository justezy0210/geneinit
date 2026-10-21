import tempfile
from pathlib import Path
import unittest

from geneinit.models.candidate import (
    CdsPart,
    audit_orf,
    build_candidate_universe,
    parse_attributes,
    parse_predictor_gtf,
    read_fasta,
    structure_identifier,
)


class CandidateModelTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        sequence = list("N" * 160)
        sequence[0:4] = "ATGA"
        sequence[10:18] = "AACCCTAA"
        sequence[100:109] = "TTATTTCAT"
        self.genome_path = self.root / "genome.fa"
        self.genome_path.write_text(">chr1\n" + "".join(sequence) + "\n")
        self.genome = read_fasta(self.genome_path)

    def tearDown(self):
        self.temporary.cleanup()

    def test_attributes_accept_gtf_gff_and_augustus_bare_id(self):
        attrs, bare = parse_attributes('gene_id "g1"; transcript_id "g1.t1";')
        self.assertEqual(attrs["gene_id"], "g1")
        self.assertEqual(attrs["transcript_id"], "g1.t1")
        self.assertIsNone(bare)
        attrs, bare = parse_attributes("ID=tx1;Parent=g1")
        self.assertEqual(attrs, {"ID": "tx1", "Parent": "g1"})
        self.assertIsNone(bare)
        attrs, bare = parse_attributes("g1.t1")
        self.assertEqual(attrs, {})
        self.assertEqual(bare, "g1.t1")

    def test_augustus_and_genemark_same_cds_collapse_by_structure(self):
        augustus = self.root / "augustus.gtf"
        augustus.write_text(
            "chr1\tAUGUSTUS\tgene\t1\t18\t.\t+\t.\tg1\n"
            "chr1\tAUGUSTUS\ttranscript\t1\t18\t.\t+\t.\tg1.t1\n"
            'chr1\tAUGUSTUS\texon\t1\t4\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";\n'
            'chr1\tAUGUSTUS\tCDS\t1\t4\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
            'chr1\tAUGUSTUS\texon\t11\t18\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";\n'
            'chr1\tAUGUSTUS\tCDS\t11\t18\t.\t+\t2\ttranscript_id "g1.t1"; gene_id "g1";\n'
        )
        genemark = self.root / "genemark.gtf"
        genemark.write_text(
            'chr1\tGeneMark.hmm3\tCDS\t1\t4\t.\t+\t0\tgene_id "gm42"; transcript_id "gm42.t1";\n'
            'chr1\tGeneMark.hmm3\tCDS\t11\t18\t.\t+\t2\tgene_id "gm42"; transcript_id "gm42.t1";\n'
        )
        august_models = parse_predictor_gtf(augustus, "augustus", self.genome)
        gm_models = parse_predictor_gtf(genemark, "genemark_etp", self.genome)
        self.assertEqual(august_models[0].structure_id, gm_models[0].structure_id)
        self.assertNotEqual(august_models[0].source_model_id, gm_models[0].source_model_id)
        self.assertEqual(august_models[0].orf.status, "complete")
        self.assertEqual(august_models[0].orf.coding_sequence, "ATGAAACCCTAA")
        universe = build_candidate_universe(august_models + gm_models)
        self.assertEqual(len(universe), 1)
        self.assertEqual({model.predictor for model in universe[0].models}, {"augustus", "genemark_etp"})

    def test_direct_augustus_gff3_parent_ids_are_parsed(self):
        augustus = self.root / "augustus.gff3"
        augustus.write_text(
            "##gff-version 3\n"
            "chr1\tAUGUSTUS\tgene\t1\t18\t.\t+\t.\tID=g1\n"
            "chr1\tAUGUSTUS\tmRNA\t1\t18\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "chr1\tAUGUSTUS\texon\t1\t4\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1\n"
            "chr1\tAUGUSTUS\tCDS\t1\t4\t.\t+\t0\tID=g1.t1.c1;Parent=g1.t1\n"
            "chr1\tAUGUSTUS\texon\t11\t18\t.\t+\t.\tID=g1.t1.e2;Parent=g1.t1\n"
            "chr1\tAUGUSTUS\tCDS\t11\t18\t.\t+\t2\tID=g1.t1.c2;Parent=g1.t1\n"
        )
        model = parse_predictor_gtf(augustus, "augustus", self.genome)[0]
        self.assertEqual(model.source_gene_id, "g1")
        self.assertEqual(model.source_transcript_id, "g1.t1")
        self.assertEqual(model.orf.status, "complete")

    def test_minus_strand_cds_is_oriented_in_transcript_direction(self):
        path = self.root / "minus.gtf"
        path.write_text(
            'chr1\tGeneMark.hmm3\tCDS\t101\t109\t.\t-\t0\tgene_id "gm_minus"; transcript_id "gm_minus.t1";\n'
        )
        model = parse_predictor_gtf(path, "genemark_etp", self.genome)[0]
        self.assertEqual(model.orf.coding_sequence, "ATGAAATAA")
        self.assertEqual(model.orf.protein_sequence, "MK*")
        self.assertEqual(model.orf.status, "complete")

    def test_structure_id_ignores_input_order_but_preserves_phase(self):
        first = CdsPart(1, 4, 0)
        second = CdsPart(11, 18, 2)
        self.assertEqual(
            structure_identifier("chr1", "+", (first, second)),
            structure_identifier("chr1", "+", (second, first)),
        )
        self.assertNotEqual(
            structure_identifier("chr1", "+", (first, second)),
            structure_identifier("chr1", "+", (first, CdsPart(11, 18, 0))),
        )

    def test_orf_audit_flags_phase_inconsistency_without_dropping_model(self):
        audit = audit_orf(
            "chr1",
            "+",
            (CdsPart(1, 4, 0), CdsPart(11, 18, 0)),
            self.genome,
        )
        self.assertFalse(audit.phase_consistent)
        self.assertIn("phase_inconsistent", audit.status)
        self.assertEqual(audit.coding_sequence, "ATGAAACCCTAA")


if __name__ == "__main__":
    unittest.main()
