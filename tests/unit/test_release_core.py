import csv
import json
from pathlib import Path
import tempfile
import unittest

from geneinit.release.core import ReleaseError, build_structural_release


DECISION_COLUMNS = (
    "transcript_id", "candidate_id", "source_status", "input_gene_id",
    "postqc_gene_id", "canonical_transcript", "canonical_rank", "seqid",
    "strand", "start", "end", "selection_locus_id", "routes",
    "external_support", "operation", "reason",
)


class StructuralReleaseTests(unittest.TestCase):
    def _write_inputs(self, root: Path, *, complete_orf: bool = True) -> dict[str, Path]:
        sequence = list("C" * 180)
        sequence[9:99] = list(("ATG" if complete_orf else "CTG") + "GCC" * 28 + "TAA")
        genome = root / "genome.fa"
        genome.write_text(">chr1\n" + "".join(sequence) + "\n", encoding="utf-8")
        gff3 = root / "postqc.gff3"
        gff3.write_text(
            "##gff-version 3\n"
            "chr1\tGeneInIt\tgene\t7\t105\t.\t+\t.\tID=GPQC-1;\n"
            "chr1\tGeneInIt\tmRNA\t7\t105\t.\t+\t.\tID=tmp1;Parent=GPQC-1;candidate_id=C1;\n"
            "chr1\tGeneInIt\texon\t7\t105\t.\t+\t.\tID=e1;Parent=tmp1;\n"
            "chr1\tGeneInIt\tCDS\t10\t99\t.\t+\t0\tID=c1;Parent=tmp1;\n"
            "chr1\tGeneInIt\tfive_prime_UTR\t7\t9\t.\t+\t.\tID=u5;Parent=tmp1;\n"
            "chr1\tGeneInIt\tthree_prime_UTR\t100\t105\t.\t+\t.\tID=u3;Parent=tmp1;\n",
            encoding="utf-8",
        )
        decisions = root / "decisions.tsv"
        with decisions.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(DECISION_COLUMNS)
            writer.writerow((
                "tmp1", "C1", "RELEASED", "GGRP-1", "GPQC-1", 1, 1,
                "chr1", "+", 7, 105, "L1", "STRONG_RECURRENT_RNA", 1,
                "PRESERVE_GENE", "FIXED_POINT_COMPONENT",
            ))
        edges = root / "edges.tsv"
        edges.write_text(
            "left_transcript_id\tright_transcript_id\tseqid\tspan_overlap_bp\tdisposition\n",
            encoding="utf-8",
        )
        bridge_reviews = root / "bridge_reviews.tsv"
        bridge_reviews.write_text(
            "left_transcript_id\tright_transcript_id\tseqid\tintron_start\t"
            "intron_end\tintron_bp\tthreshold_class\tstatus\trequired_evidence\n",
            encoding="utf-8",
        )
        policy = root / "policy.json"
        policy.write_text(
            json.dumps({"policy_id": "fixture-policy", "status": "MIGRATION_PENDING"}),
            encoding="utf-8",
        )
        return {
            "genome": genome, "gff3": gff3, "decisions": decisions,
            "edges": edges, "bridge_reviews": bridge_reviews, "policy": policy,
        }

    def _run(self, root: Path, inputs: dict[str, Path]) -> dict[str, Path]:
        output = root / "release"
        paths = {
            name: output / filename for name, filename in {
                "annotation": "annotation.gff3", "transcripts": "transcripts.fna",
                "cds": "cds.fna", "proteins": "proteins.faa",
                "longest": "longest_proteins.faa", "gene_map": "gene_id_map.tsv",
                "transcript_map": "transcript_id_map.tsv", "lineage": "candidate_lineage.tsv",
                "decisions_out": "structural_decisions.tsv", "review": "review_flags.tsv",
                "qc_tsv": "structural_qc.tsv", "qc_json": "structural_qc.json",
                "status": "release_status.json", "manifest": "release_manifest.json",
                "checksums": "checksums.sha256", "policy_out": "policy.json",
                "complete": "release.complete.tsv",
            }.items()
        }
        build_structural_release(
            inputs["gff3"], inputs["decisions"], inputs["edges"],
            inputs["bridge_reviews"], inputs["genome"],
            inputs["policy"], ["chr1"], "TEST", "TESTG",
            paths["annotation"], paths["transcripts"], paths["cds"], paths["proteins"],
            paths["longest"], paths["gene_map"], paths["transcript_map"], paths["lineage"],
            paths["decisions_out"], paths["review"], paths["qc_tsv"], paths["qc_json"],
            paths["status"], paths["manifest"], paths["checksums"], paths["policy_out"],
            paths["complete"],
        )
        return paths

    def test_builds_deterministic_structural_release(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._run(root, self._write_inputs(root))
            self.assertIn("ID=TESTG000001;", paths["annotation"].read_text())
            self.assertIn("ID=TESTG000001.t1;Parent=TESTG000001;", paths["annotation"].read_text())
            self.assertIn(">TESTG000001.t1\n", paths["proteins"].read_text())
            self.assertIn("PASS_MIGRATION_PENDING", paths["complete"].read_text())
            status = json.loads(paths["status"].read_text())
            self.assertFalse(status["production_publication_authorized"])
            manifest = json.loads(paths["manifest"].read_text())
            self.assertEqual(manifest["gene_count"], 1)
            self.assertIn("release_status.json", manifest["files"])

    def test_rejects_incomplete_final_orf(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(ReleaseError, "non-complete final ORF"):
                self._run(root, self._write_inputs(root, complete_orf=False))

    def test_rejects_duplicate_overlapping_and_adjacent_exons(self):
        for kind, extra in [("duplicate", (7,105)), ("overlap", (30,70)), ("adjacent", None)]:
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                inputs = self._write_inputs(root)
                text = inputs["gff3"].read_text()
                if extra:
                    start, end = extra
                    text += f"chr1\tGeneInIt\texon\t{start}\t{end}\t.\t+\t.\tID=e2;Parent=tmp1;\n"
                else:
                    text = text.replace("exon\t7\t105", "exon\t7\t55")
                    text += "chr1\tGeneInIt\texon\t56\t105\t.\t+\t.\tID=e2;Parent=tmp1;\n"
                inputs["gff3"].write_text(text)
                with self.assertRaisesRegex(ReleaseError, "duplicate, overlapping or adjacent final exons"):
                    self._run(root, inputs)

    def test_rejects_final_utr_gene_hierarchy_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = self._write_inputs(root)
            inputs["gff3"].write_text(
                inputs["gff3"].read_text(encoding="utf-8").replace(
                    "GPQC-1", "GPQC-drift"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ReleaseError, "gene hierarchy differs from post-QC"
            ):
                self._run(root, inputs)

    def test_propagates_held_bridge_to_release_review_flags(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = self._write_inputs(root)
            with inputs["bridge_reviews"].open("a", encoding="utf-8") as handle:
                handle.write(
                    "tmp1\ttmp-review\tchr1\t106\t7105\t7000\t"
                    "ADDITIONAL_EVIDENCE\t"
                    "HELD_LONG_INTRON_ADDITIONAL_EVIDENCE_REQUIRED\t"
                    "SPANNING_PROTEIN_OR_FULL_LENGTH_TRANSCRIPT\n"
                )
            paths = self._run(root, inputs)
            review_text = paths["review"].read_text(encoding="utf-8")
            self.assertIn("HELD_LONG_INTRON_ADDITIONAL_EVIDENCE_REQUIRED", review_text)
            self.assertIn("TESTG000001.t1", review_text)


if __name__ == "__main__":
    unittest.main()
