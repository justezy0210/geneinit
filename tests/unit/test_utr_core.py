import csv
from pathlib import Path
import tempfile
import unittest

from geneinit.utr.core import (
    UtrError,
    aggregate_utr_shards,
    normalize_refinement,
    prepare_postqc_utr_shards,
    prepare_selected_coding_models,
    restore_postqc_gene_hierarchy,
    write_normalized_shard,
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def table(path: Path, columns, rows) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    return path


CODING_GTF = """chr1\tGeneInIt\tgene\t100\t300\t.\t+\t.\tgene_id \"GITMPG0000001\";
chr1\tGeneInIt\ttranscript\t100\t300\t.\t+\t.\tgene_id \"GITMPG0000001\"; transcript_id \"GITMPG0000001.t1\"; candidate_id \"C1\";
chr1\tGeneInIt\texon\t100\t150\t.\t+\t.\tgene_id \"GITMPG0000001\"; transcript_id \"GITMPG0000001.t1\"; candidate_id \"C1\";
chr1\tGeneInIt\texon\t201\t300\t.\t+\t.\tgene_id \"GITMPG0000001\"; transcript_id \"GITMPG0000001.t1\"; candidate_id \"C1\";
chr1\tGeneInIt\tCDS\t100\t150\t.\t+\t0\tgene_id \"GITMPG0000001\"; transcript_id \"GITMPG0000001.t1\"; candidate_id \"C1\";
chr1\tGeneInIt\tCDS\t201\t300\t.\t+\t0\tgene_id \"GITMPG0000001\"; transcript_id \"GITMPG0000001.t1\"; candidate_id \"C1\";
"""


RAW_GFF3 = """##gff-version 3
chr1\tGeneInIt\tgene\t90\t330\t.\t+\t.\tID=GITMPG0000001;
chr1\tGeneInIt\tmRNA\t90\t330\t.\t+\t.\tID=GITMPG0000001.t1;Parent=GITMPG0000001;candidate_id=C1;utr_refine_evidence=MSTRG.1.1;
chr1\tRNA\texon\t90\t150\t.\t+\t.\tParent=GITMPG0000001.t1;
chr1\tRNA\texon\t201\t330\t.\t+\t.\tParent=GITMPG0000001.t1;
chr1\tGeneInIt\tCDS\t100\t150\t.\t+\t0\tParent=GITMPG0000001.t1;
chr1\tGeneInIt\tCDS\t201\t300\t.\t+\t0\tParent=GITMPG0000001.t1;
chr1\tInGenAnnot\tfive_prime_UTR\t90\t99\t.\t+\t.\tParent=GITMPG0000001.t1;
chr1\tInGenAnnot\tthree_prime_UTR\t301\t330\t.\t+\t.\tParent=GITMPG0000001.t1;
"""


PREP_COLUMNS = (
    "candidate_id", "temporary_gene_id", "temporary_transcript_id", "seqid",
    "strand", "start", "end", "selection_locus_id", "predictors",
    "source_gene_ids", "source_transcript_ids", "routes", "active_fusion",
)


class UtrCoreTests(unittest.TestCase):
    def test_normalizes_repeated_shared_exon_per_transcript(self):
        second = CODING_GTF.replace("GITMPG0000001", "GITMPG0000002").replace('"C1"', '"C2"')
        raw_second = RAW_GFF3.replace("GITMPG0000001", "GITMPG0000002").replace("candidate_id=C1", "candidate_id=C2")
        parents = "Parent=GITMPG0000001.t1,GITMPG0000002.t1;"
        lines = []
        for line in (RAW_GFF3 + raw_second).splitlines():
            if "\texon\t" in line:
                line = line.replace("Parent=GITMPG0000001.t1;", parents).replace("Parent=GITMPG0000002.t1;", parents)
            lines.append(line)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, models = normalize_refinement(write(root/"coding.gtf", CODING_GTF+second), write(root/"raw.gff3", "\n".join(lines)+"\n"))
            self.assertEqual(len(models), 2)
            for model in models.values():
                self.assertEqual([(f.start, f.end) for f in model.exons], [(90,150),(201,330)])
                self.assertEqual(len(model.cds), 2)

    def test_normalizes_overlapping_exon_blocks_on_both_strands(self):
        line = "chr1\tRNA\texon\t90\t150\t.\t+\t.\tParent=GITMPG0000001.t1;"
        split = line.replace("90\t150", "90\t125")+"\n"+line.replace("90\t150", "120\t150")
        for strand in ("+", "-"):
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                coding = CODING_GTF.replace("\t+\t", f"\t{strand}\t")
                raw = RAW_GFF3.replace(line, split).replace("\t+\t", f"\t{strand}\t")
                _, models = normalize_refinement(write(root/"coding.gtf", coding), write(root/"raw.gff3", raw))
                model = models["GITMPG0000001.t1"]
                self.assertEqual([(f.start, f.end) for f in model.exons], [(90,150),(201,330)])
                self.assertEqual([(f.start, f.end, f.phase) for f in model.cds], [(100,150,"0"),(201,300,"0")])

    def test_repeated_cds_is_not_silently_deduplicated(self):
        duplicate = next(line for line in RAW_GFF3.splitlines() if "\tCDS\t" in line)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(UtrError, "changed CDS coordinates or phase"):
                normalize_refinement(write(root/"coding.gtf", CODING_GTF), write(root/"raw.gff3", RAW_GFF3+duplicate+"\n"))

    def test_selected_coding_preparation_writes_only_selected_candidates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidates = table(
                root / "candidates.tsv",
                (
                    "candidate_id", "seqid", "strand", "start", "end", "predictors",
                    "source_gene_ids", "source_transcript_ids", "cds_signature",
                ),
                (
                    ("C1", "chr1", "+", 100, 300, "augustus", "g1", "t1", "100-150:0,201-300:0"),
                    ("C2", "chr2", "+", 400, 500, "genemark", "g2", "t2", "400-500:0"),
                ),
            )
            audit = table(
                root / "audit.tsv",
                ("candidate_id", "locus_id", "selected", "routes", "active_fusion"),
                (("C1", "L1", 1, "RNA", 0), ("C2", "L2", 0, ".", 0)),
            )
            rows = prepare_selected_coding_models(
                candidates, audit, ["chr1", "chr2"],
                root / "coding.gff3", root / "mapping.tsv",
            )
            self.assertEqual([row["candidate_id"] for row in rows], ["C1"])
            self.assertIn("candidate_id=C1", (root / "coding.gff3").read_text())
            self.assertNotIn("candidate_id=C2", (root / "coding.gff3").read_text())

    def test_postqc_preparation_includes_promoted_isoform_and_writes_gene_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            postqc = write(
                root / "postqc.gff3",
                "##gff-version 3\n"
                "chr1\tGeneInIt-postQC\tgene\t100\t330\t.\t+\t.\tID=GPQC-1;\n"
                "chr1\tGeneInIt\tmRNA\t100\t300\t.\t+\t.\tID=t1;Parent=GPQC-1;candidate_id=C1;\n"
                "chr1\tGeneInIt\texon\t100\t300\t.\t+\t.\tParent=t1;\n"
                "chr1\tGeneInIt\tCDS\t100\t300\t.\t+\t0\tParent=t1;\n"
                "chr1\tGeneInIt-RNA-bridge\tmRNA\t120\t330\t.\t+\t.\tID=t2;Parent=GPQC-1;candidate_id=GRNA-1;\n"
                "chr1\tGeneInIt-RNA-bridge\texon\t120\t330\t.\t+\t.\tParent=t2;\n"
                "chr1\tGeneInIt-RNA-bridge\tCDS\t120\t330\t.\t+\t0\tParent=t2;\n",
            )
            decisions = table(
                root / "decisions.tsv",
                (
                    "transcript_id", "candidate_id", "postqc_gene_id",
                    "selection_locus_id", "routes", "source_status",
                ),
                (
                    ("t1", "C1", "GPQC-1", "L1", "RNA", "RELEASED"),
                    ("t2", "GRNA-1", "GPQC-1", ".", "POSTQC_STRONG_RNA_BRIDGE_COMPLETE_ORF", "PROMOTED_STRONG_RNA_BRIDGE"),
                ),
            )
            candidates = table(
                root / "candidates.tsv",
                ("candidate_id", "predictors", "source_gene_ids", "source_transcript_ids"),
                (("C1", "augustus", "g1", "a1"),),
            )
            audit = table(
                root / "audit.tsv", ("candidate_id", "active_fusion"),
                (("C1", 0),),
            )
            rna = write(
                root / "rna.gtf",
                'chr1\tStringTie\ttranscript\t90\t340\t.\t+\t.\tgene_id "M1"; transcript_id "M1.1";\n',
            )
            coding_paths = {
                "000001": root / "chr1.coding.gtf",
                "000002": root / "chr2.coding.gtf",
            }
            rna_paths = {
                "000001": root / "chr1.rna.gtf",
                "000002": root / "chr2.rna.gtf",
            }
            rows = prepare_postqc_utr_shards(
                postqc, decisions, candidates, audit, rna,
                {"chr1": "000001", "chr2": "000002"},
                coding_paths, rna_paths, root / "mapping.tsv", root / "tasks.tsv",
            )
            self.assertEqual({row["candidate_id"] for row in rows}, {"C1", "GRNA-1"})
            coding_text = coding_paths["000001"].read_text()
            self.assertEqual(coding_text.count("\tgene\t"), 2)
            self.assertIn('transcript_id "t1"', coding_text)
            self.assertIn('transcript_id "t2"', coding_text)
            order, models = normalize_refinement(coding_paths["000001"], None)
            prep = {row["candidate_id"]: row for row in rows}
            restore_postqc_gene_hierarchy(order, models, prep)
            write_normalized_shard(
                root / "normalized.gff3", root / "normalized.tsv", order, models, prep
            )
            self.assertEqual((root / "normalized.gff3").read_text().count("\tgene\t"), 1)

    def test_normalization_adds_boundary_utrs_without_changing_cds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            coding = write(root / "coding.gtf", CODING_GTF)
            raw = write(root / "raw.gff3", RAW_GFF3)
            order, models = normalize_refinement(coding, raw)
            model = models[order[0]]
            self.assertEqual([(item.start, item.end, item.phase) for item in model.cds], [(100, 150, "0"), (201, 300, "0")])
            self.assertEqual([(item.start, item.end) for item in model.five_prime_utrs], [(90, 99)])
            self.assertEqual([(item.start, item.end) for item in model.three_prime_utrs], [(301, 330)])

    def test_normalization_rejects_any_cds_phase_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            coding = write(root / "coding.gtf", CODING_GTF)
            raw = write(root / "raw.gff3", RAW_GFF3.replace("201\t300\t.\t+\t0", "201\t300\t.\t+\t1"))
            with self.assertRaisesRegex(UtrError, "changed CDS coordinates or phase"):
                normalize_refinement(coding, raw)

    def test_postqc_gene_restore_fails_closed_on_missing_assignment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            coding = write(root / "coding.gtf", CODING_GTF)
            order, models = normalize_refinement(coding, None)
            row = {column: "." for column in PREP_COLUMNS}
            row.update(
                {
                    "candidate_id": "C1",
                    "temporary_gene_id": "GITMPG0000001",
                    "temporary_transcript_id": "GITMPG0000001.t1",
                }
            )
            with self.assertRaisesRegex(UtrError, "missing final post-QC gene"):
                restore_postqc_gene_hierarchy(order, models, {"C1": row})

    def test_normalization_rejects_non_cds_exon_inside_coding_span(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            coding = write(root / "coding.gtf", CODING_GTF)
            raw_text = RAW_GFF3.replace(
                "chr1\tRNA\texon\t90\t150",
                "chr1\tRNA\texon\t90\t170",
            )
            raw = write(root / "raw.gff3", raw_text)
            with self.assertRaisesRegex(UtrError, "inside CDS span"):
                normalize_refinement(coding, raw)

    def test_shard_writer_and_aggregate_preserve_candidate_membership(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            coding = write(root / "coding.gtf", CODING_GTF)
            raw = write(root / "raw.gff3", RAW_GFF3)
            order, models = normalize_refinement(coding, raw)
            prep_row = {
                "candidate_id": "C1", "temporary_gene_id": "GITMPG0000001",
                "temporary_transcript_id": "GITMPG0000001.t1", "seqid": "chr1",
                "strand": "+", "start": "100", "end": "300", "selection_locus_id": "L1",
                "predictors": "augustus", "source_gene_ids": "g1",
                "source_transcript_ids": "t1", "routes": "RNA", "active_fusion": "0",
            }
            normalized = root / "normalized.gff3"
            mapping = root / "normalized.mapping.tsv"
            write_normalized_shard(normalized, mapping, order, models, {"C1": prep_row})
            preparation = table(root / "preparation.tsv", PREP_COLUMNS, (tuple(prep_row[column] for column in PREP_COLUMNS),))
            aggregate_utr_shards(
                [normalized], [mapping], preparation,
                root / "all.gff3", root / "all.mapping.tsv",
                root / "summary.tsv", root / "complete.tsv",
            )
            self.assertIn("five_prime_UTR", (root / "all.gff3").read_text())
            self.assertIn("PASS_MIGRATION_PENDING", (root / "complete.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
