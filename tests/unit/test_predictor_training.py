from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

from geneinit.models.candidate import reverse_complement
from geneinit.predictor.common import PredictorError
from geneinit.predictor.training import (
    TrainingRecord, evaluation_metrics, nonredundant_records,
    read_training_records, require_clean_training, split_training,
    require_matching_selection,
    stratified_sample, training_errors,
)


def genbank(sequence, location, name="chr1", start=1):
    locus = f"{name}_{start}-{start + len(sequence) - 1}"
    return (f"LOCUS       {locus} {len(sequence)} bp DNA\n"
            f"FEATURES             Location/Qualifiers\n"
            f"     CDS             {location}\n"
            f'                     /gene="{locus}.gene"\n'
            f"ORIGIN\n        1 {sequence}\n//\n")


def record(index):
    locus = f"chr{index}_1-100"
    return TrainingRecord(locus, f"gene{index}", genbank("ATGAAATAA", "1..9", f"chr{index}"),
                          f"M{index}K", 300 * (index % 12 + 1), 100, index % 16 + 1,
                          0.3 + (index % 17) / 100, 10, 0, "", f"chr{index}", 1, 100)


class TrainingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.hits = self.root / "hits.tsv"
        self.hits.write_text("")

    def parse(self, text):
        path = self.root / "input.gb"
        path.write_text(text)
        return read_training_records(path)

    def test_reviewed_selection_gate_rejects_same_size_changed_models_or_roles(self):
        expected = self.root / "reviewed.tsv"
        observed = self.root / "regenerated.tsv"
        original = "locus\tgene\trole\nchr1\tgene1\tholdout\nchr2\tgene2\tonlytrain\n"
        expected.write_text(original)
        observed.write_text(original)
        require_matching_selection(observed, expected)
        for changed in (
            original.replace("gene1", "geneX"),
            original.replace("holdout", "SWAP").replace("onlytrain", "holdout").replace("SWAP", "onlytrain"),
        ):
            observed.write_text(changed)
            with self.assertRaisesRegex(PredictorError, "differs from the reviewed reference"):
                require_matching_selection(observed, expected)

    def test_normalize_only_dna_and_retain_gene_identifier(self):
        result = self.parse(genbank("cccATGaaaTAAccc", "4..12", "ChrX"))[0]
        self.assertEqual(result.protein, "MK")
        self.assertEqual(result.reason, "")
        self.assertIn("cccatgaaataaccc", result.raw)
        self.assertIn('gene="ChrX_1-15.gene"', result.raw)

    def test_wrapped_spliced_cds_both_strands(self):
        sequence = "CCCATGGTAAAGAAATAACCC"
        plus = genbank(sequence, "join(4..6,\n                     13..18)")
        minus = genbank(reverse_complement(sequence), "complement(join(4..9,16..18))", "chr2")
        for item in self.parse(plus + minus):
            self.assertEqual((item.protein, item.reason, item.exons), ("MK", "", 2))

    def test_bad_orf_and_splice_are_excluded(self):
        result = self.parse(genbank("ATGTAATAA", "1..9") +
                            genbank("ATGCCAAAGAAATAA", "join(1..3,10..15)", "chr2"))
        self.assertIn("internal_stop", result[0].reason)
        self.assertEqual(result[1].reason, "noncanonical_splice")

    def test_split_start_and_stop_codons_are_not_training_compatible(self):
        coding = "ATGAAATAA"
        for boundary, cut in (("start", 2), ("stop", 7)):
            sequence = "CCC" + coding[:cut] + "GTAAAG" + coding[cut:] + "CCC"
            coords = [(4, 3 + cut), (10 + cut, len(sequence) - 3)]
            for strand in ("+", "-"):
                with self.subTest(boundary=boundary, strand=strand):
                    oriented = sequence
                    parts = coords
                    if strand == "-":
                        oriented = reverse_complement(sequence)
                        parts = sorted((len(sequence) - b + 1, len(sequence) - a + 1)
                                       for a, b in coords)
                    location = "join(" + ",".join(f"{a}..{b}" for a, b in parts) + ")"
                    if strand == "-":
                        location = f"complement({location})"
                    result = self.parse(genbank(oriented, location))[0]
                    self.assertEqual(result.protein, "MK")
                    self.assertEqual(result.reason, f"unsupported_split_{boundary}_codon")

    def test_duplicate_or_incomplete_record_fails(self):
        text = genbank("ATGAAATAA", "1..9")
        with self.assertRaises(PredictorError):
            self.parse(text + text)
        with self.assertRaises(PredictorError):
            self.parse(text[:-3])

    def test_transitive_similarity_and_both_coverage_thresholds(self):
        items = [record(i) for i in range(4)]
        a, b, c, d = [item.locus for item in items]
        self.hits.write_text(f"{a}\t{b}\t90\t1\t90\t100\t1\t90\t100\n"
                             f"{b}\t{c}\t80\t1\t80\t100\t1\t80\t100\n"
                             f"{c}\t{d}\t99\t1\t90\t100\t1\t90\t200\n")
        representatives, membership = nonredundant_records(items, self.hits)
        self.assertEqual(len(representatives), 2)
        self.assertEqual(membership[a], membership[c])
        self.assertNotEqual(membership[a], membership[d])

    def test_exact_protein_and_overlapping_dna_grouping(self):
        a, b, c = [record(i) for i in range(3)]
        b = replace(b, protein=a.protein)
        c = replace(c, chromosome=b.chromosome, start=50, end=149)
        representatives, groups = nonredundant_records([a, b, c], self.hits)
        self.assertEqual(len(representatives), 1)
        self.assertEqual(len(set(groups.values())), 1)

    def test_stratified_selection_is_order_independent(self):
        items = [record(i) for i in range(200)]
        first, rest = stratified_sample(items, 53, 17, (0.35, 0.41))
        second, _ = stratified_sample(list(reversed(items)), 53, 17, (0.35, 0.41))
        self.assertEqual(first, second)
        self.assertEqual((len(first), len(rest)), (53, 147))
        self.assertFalse({r.locus for r in first} & {r.locus for r in rest})

    def test_split_disjointness_cap_and_exclusion_audit(self):
        items = [record(i) for i in range(30)]
        items[1] = replace(items[1], protein=items[0].protein)
        items[2] = replace(items[2], reason="noncanonical_splice")
        summary = split_training(items, self.hits, self.root, {items[3].locus: "bad structure"},
                                 max_loci=20, optimize_loci=4, holdout_loci=5)
        self.assertEqual(summary["nonredundant"]["loci"], 27)
        self.assertEqual([summary["splits"][k]["loci"] for k in ("holdout", "optimize", "onlytrain")],
                         [5, 4, 11])
        self.assertEqual(summary["cross_split_group_overlap"], 0)
        def loci(role):
            return {line.split()[1] for line in (self.root / f"{role}.gb").read_text().splitlines()
                    if line.startswith("LOCUS")}
        self.assertEqual(loci("train"), loci("optimize") | loci("onlytrain"))
        self.assertFalse(loci("train") & loci("holdout"))
        self.assertIn("noncanonical_splice", (self.root / "selection.tsv").read_text())

    def test_too_few_independent_records_fails(self):
        items = [replace(record(i), protein="MK") for i in range(30)]
        with self.assertRaises(PredictorError):
            split_training(items, self.hits, self.root, {}, max_loci=20, optimize_loci=4, holdout_loci=5)

    def test_training_exit_zero_does_not_hide_locus_errors(self):
        path = self.root / "training.log"
        path.write_text("exonmodel.cc: ERROR in sequence chr1_1-9: exon doesn't end in stop codon\n")
        self.assertIn("chr1_1-9", training_errors(path))
        with self.assertRaises(PredictorError):
            require_clean_training(path)

    def test_multiple_locus_warnings_preserve_stop_gate_and_all_reasons(self):
        path = self.root / "training.log"
        stop = "exon doesn't end in stop codon"
        splice = "intron has non-standard splice sites"
        path.write_text(
            f"exonmodel.cc: ERROR in sequence chr1_1-9: {stop}\n"
            f"intronmodel.cc: ERROR in sequence chr1_1-9: {splice}\n"
            f"exonmodel.cc: ERROR in sequence chr1_1-9: {stop}\n"
            "exonmodel.cc: ERROR in sequence chr2_1-9: missing start codon\n"
        )
        errors = training_errors(path)
        self.assertEqual(len(errors), 2)
        self.assertEqual(errors["chr1_1-9"], f"{stop}; {splice}")
        self.assertEqual(errors["chr2_1-9"], "missing start codon")
        self.assertTrue(any("end in stop codon" in reason for reason in errors.values()))
        with self.assertRaises(PredictorError):
            require_clean_training(path)

    def test_augustus_accuracy_table_and_incomplete_output(self):
        path = self.root / "evaluation.log"
        path.write_text("nucleotide level | 0.9 | 0.8 |\n"
                        "exon level | 100 | 90 | -- | 0.75 | 0.6 |\n"
                        "gene level | 50 | 45 | 30 | 0.6 | 0.5 |\n")
        metrics = evaluation_metrics(path)
        self.assertAlmostEqual(metrics["gene"]["f1"], 0.6 / 1.1)
        self.assertEqual(metrics["nucleotide"]["sensitivity"], 0.9)
        path.write_text("nucleotide level | 0.9 | 0.8 |\n")
        with self.assertRaises(PredictorError):
            evaluation_metrics(path)


if __name__ == "__main__":
    unittest.main()
