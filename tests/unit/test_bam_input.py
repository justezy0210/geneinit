import csv
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest

import pysam

from geneinit.contracts.alignments import (
    AlignmentContractError, check_staging_paths, stage_alignment, validate_alignment,
)
from geneinit.contracts.preflight import PreflightError, run_input_preflight, run_preflight
from geneinit.contracts.samples import BAM_SAMPLE_COLUMNS, SampleContractError, load_samples


class BamInputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.lengths = {"ctg1": 1000, "ctg2": 100}
        self.genome = self.root / "genome.fa"
        self.reference = self.root / "mapping.fa"
        self.genome.write_text(">ctg1\n" + "acgt" * 250 + "\n>ctg2\n" + "ACGT" * 25 + "\n")
        self.reference.write_text(self.genome.read_text().upper().replace("CTG", "ctg"))
        self.bam = self.root / "S1.bam"
        self.csi = self.root / "S1.bam.csi"
        self.write_bam()
        self.samples = self.root / "samples.tsv"
        self.write_samples()
        self.repeats = self.root / "repeat.gff3"
        self.repeats.write_text("ctg1\tEDTA\trepeat_region\t1\t4\t.\t+\t.\tID=r1\n")
        self.proteins = self.root / "protein.fa"
        self.proteins.write_text(">P1\nMPEPTIDE\n")
        self.policy = self.root / "policy.json"
        self.policy.write_text('{"policy_id":"test"}\n')

    def tearDown(self):
        self.temporary.cleanup()

    def write_bam(self, lengths=None, paired=True, sort="coordinate", start=10, query="r1"):
        lengths = self.lengths if lengths is None else lengths
        with pysam.AlignmentFile(str(self.bam), "wb", header={
            "HD": {"VN": "1.6", "SO": sort},
            "SQ": [{"SN": key, "LN": value} for key, value in lengths.items()],
        }) as handle:
            read = pysam.AlignedSegment()
            read.query_name = query
            read.query_sequence = "ACGT" * 10
            read.flag = 99 if paired else 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = "20M50N20M"
            read.next_reference_id = 0 if paired else -1
            read.next_reference_start = start + 100 if paired else -1
            read.set_tag("XS", "+")
            handle.write(read)
        pysam.index("-c", str(self.bam), str(self.csi))

    def write_samples(self, layout="PE", strand="RF", duplicate=False):
        with self.samples.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(BAM_SAMPLE_COLUMNS)
            writer.writerow(["S1", "C", "1", layout, strand, self.bam, self.csi])
            if duplicate:
                writer.writerow(["S2", "C", "2", layout, strand, self.bam, self.csi])

    def record(self):
        return load_samples(self.samples, input_mode="bam")[0]

    def config(self):
        return dict(assembly_id="ASM1", organism="Genus species", assembly_scope="nuclear",
                    reference_representation="haploid", expected_protein_copies=1,
                    genetic_code=1, gene_prefix="ASM1G", genome=str(self.genome),
                    repeat_gff=str(self.repeats), protein_fasta=str(self.proteins),
                    samples=str(self.samples), policy=str(self.policy),
                    rna={"input_mode": "bam", "reference_fasta": str(self.reference)})

    def preflight(self):
        return run_preflight(self.config(), self.root / "preflight")

    def stage(self):
        return stage_alignment(self.record(), self.root / "preflight/input_manifest.tsv",
                               self.root / "staged/S1.bam", self.root / "staged/S1.bam.csi",
                               self.root / "staged/import.json")

    def test_input_receipt_preserves_bam_preflight_outputs(self):
        direct = run_preflight(self.config(), self.root / "direct")
        inputs = self.root / "inputs"
        run_input_preflight(self.config(), inputs)
        cached = run_preflight(self.config(), self.root / "cached",
                               input_validation=inputs / "input_validation.json")
        self.assertEqual(direct, cached)
        for name in ("normalized_config.json", "normalized_samples.tsv", "input_manifest.tsv"):
            self.assertEqual((self.root / "direct" / name).read_bytes(),
                             (self.root / "cached" / name).read_bytes())

    def test_input_receipt_rejects_retargeted_bam_or_csi(self):
        original_bam, original_csi = self.bam, self.csi
        bam_alias, csi_alias = self.root / "alias.bam", self.root / "alias.bam.csi"
        bam_alias.symlink_to(original_bam)
        csi_alias.symlink_to(original_csi)
        self.bam, self.csi = bam_alias, csi_alias
        self.write_samples()
        inputs = self.root / "inputs"
        run_input_preflight(self.config(), inputs)
        self.bam, self.csi = self.root / "different.bam", self.root / "different.bam.csi"
        self.write_bam(query="different")
        for alias, old, new in ((bam_alias, original_bam, self.bam),
                                (csi_alias, original_csi, self.csi)):
            with self.subTest(role=alias.name):
                alias.unlink()
                alias.symlink_to(new)
                with self.assertRaisesRegex(PreflightError, "sample input paths changed.*rerun input_preflight"):
                    run_preflight(self.config(), self.root / "cached",
                                  input_validation=inputs / "input_validation.json")
                alias.unlink()
                alias.symlink_to(old)

    def test_mask_only_change_accepts_bam_without_fastqs(self):
        result = self.preflight()
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["rna_input_mode"], "bam")
        self.assertEqual(result["alignment_checks"][0]["primary_reads_checked"], 1)
        self.assertIsNone(self.record().read1)
        normalized = load_samples(self.root / "preflight/normalized_samples.tsv", input_mode="bam")
        self.assertEqual(normalized, [self.record()])
        manifest = (self.root / "preflight/input_manifest.tsv").read_text()
        self.assertIn("rna_reference_fasta", manifest)
        self.assertNotIn("read1\t", manifest)

    def test_reference_wrapping_and_description_do_not_change_identity(self):
        first = self.preflight()
        self.reference.write_text(">ctg1 mapping reference\n" + ("ACGT" * 25 + "\n") * 10 + ">ctg2\n" + "ACGT" * 25 + "\n")
        self.assertEqual(first["semantic_config_sha256"], self.preflight()["semantic_config_sha256"])

    def test_same_length_reference_base_change_rejected(self):
        self.reference.write_text(self.reference.read_text().replace("ACGT", "TCGT", 1))
        with self.assertRaisesRegex(PreflightError, "differs from genome"):
            self.preflight()

    def test_reference_required(self):
        config = self.config()
        del config["rna"]["reference_fasta"]
        with self.assertRaisesRegex(PreflightError, "requires rna.reference_fasta"):
            run_preflight(config, self.root / "bad")

    def test_unknown_mode_rejected(self):
        config = self.config()
        config["rna"]["input_mode"] = "bma"
        with self.assertRaisesRegex(PreflightError, "fastq or bam"):
            run_preflight(config, self.root / "bad")

    def test_duplicate_bam_not_a_second_replicate(self):
        self.write_samples(duplicate=True)
        with self.assertRaisesRegex(SampleContractError, "duplicate BAM"):
            self.record()

    def test_missing_csi_rejected(self):
        self.csi.unlink()
        with self.assertRaisesRegex(SampleContractError, "readable file"):
            self.record()

    def test_bad_csi_rejected(self):
        self.csi.write_bytes(b"broken index")
        with self.assertRaisesRegex(PreflightError, "invalid BAM/CSI"):
            self.preflight()

    def test_bai_renamed_to_csi_rejected(self):
        pysam.index(str(self.bam))
        shutil.copyfile(str(self.bam) + ".bai", self.csi)
        with self.assertRaises(PreflightError):
            self.preflight()

    def test_csi_from_another_bam_rejected(self):
        empty = self.root / "empty.bam"
        with pysam.AlignmentFile(str(empty), "wb", header={
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": name, "LN": length} for name, length in self.lengths.items()],
        }):
            pass
        pysam.index("-c", str(empty))
        shutil.copyfile(str(empty) + ".csi", self.csi)
        with self.assertRaisesRegex(PreflightError, "CSI cannot retrieve"):
            self.preflight()

    def test_truncated_bam_rejected(self):
        self.bam.write_bytes(self.bam.read_bytes()[:-28])
        with self.assertRaisesRegex(PreflightError, "invalid BAM/CSI"):
            self.preflight()

    def test_header_dictionary_and_order_rejected(self):
        for lengths in ({"other": 1000, "ctg2": 100}, {"ctg1": 1001, "ctg2": 100},
                        {"ctg2": 100, "ctg1": 1000}):
            with self.subTest(lengths=lengths):
                self.write_bam(lengths=lengths)
                with self.assertRaisesRegex(PreflightError, "BAM @SQ"):
                    self.preflight()

    def test_wrong_sort_order_rejected(self):
        self.write_bam(sort="queryname")
        with self.assertRaisesRegex(PreflightError, "coordinate sort"):
            self.preflight()

    def test_layout_flags_checked(self):
        self.write_bam(paired=False)
        with self.assertRaisesRegex(PreflightError, "declared layout PE"):
            self.preflight()
        self.write_samples(layout="SE", strand="R")
        self.assertEqual(self.preflight()["layouts"], ["SE"])

    def test_large_coordinate_csi_query(self):
        lengths = {"ctg1": 2**29 + 100000}
        self.write_bam(lengths=lengths, start=2**29 + 10)
        result = validate_alignment(self.record(), lengths)
        self.assertEqual(result["csi_queries_beyond_512m"], 1)

    def test_import_is_idempotent_and_source_unchanged(self):
        before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (self.bam, self.csi)}
        self.preflight()
        self.stage()
        self.stage()
        for path, (content, mtime) in before.items():
            self.assertEqual(path.read_bytes(), content)
            self.assertEqual(path.stat().st_mtime_ns, mtime)
        self.assertTrue((self.root / "staged/S1.bam").is_symlink())
        self.assertEqual((self.root / "staged/S1.bam.csi").resolve(), self.csi)

    def test_changed_input_after_preflight_rejected(self):
        self.preflight()
        self.write_bam(query="different")
        with self.assertRaisesRegex(AlignmentContractError, "content differs"):
            self.stage()
        self.assertFalse((self.root / "staged").exists())

    def test_output_source_alias_rejected(self):
        with self.assertRaisesRegex(AlignmentContractError, "aliases an original"):
            check_staging_paths(self.record(), self.bam, self.csi)
        hardlink = self.root / "alias.bam"
        hardlink.hardlink_to(self.bam)
        with self.assertRaisesRegex(AlignmentContractError, "aliases an original"):
            check_staging_paths(self.record(), hardlink, self.root / "alias.csi")

    def test_audit_and_other_sample_sources_protected(self):
        with self.assertRaisesRegex(AlignmentContractError, "aliases an original"):
            check_staging_paths(self.record(), self.root / "new.bam", self.root / "new.csi", self.bam)
        other = self.root / "other.bam"
        other.write_bytes(b"other sample")
        with self.assertRaisesRegex(AlignmentContractError, "aliases an original"):
            check_staging_paths(self.record(), other, self.root / "new.csi", other_sources=[other])

    def test_existing_destination_not_overwritten(self):
        self.preflight()
        destination = self.root / "staged/S1.bam.csi"
        destination.parent.mkdir()
        destination.write_bytes(b"keep me")
        with self.assertRaisesRegex(AlignmentContractError, "already exists"):
            self.stage()
        self.assertEqual(destination.read_bytes(), b"keep me")
        self.assertFalse((self.root / "staged/S1.bam").exists())

    def test_relocation_not_semantic_but_bam_content_is(self):
        first = self.preflight()
        moved = self.root / "moved"
        moved.mkdir()
        for attr in ("bam", "csi", "reference"):
            old = getattr(self, attr)
            new = moved / old.name
            shutil.copyfile(old, new)
            setattr(self, attr, new)
        self.write_samples()
        self.assertEqual(first["semantic_config_sha256"], self.preflight()["semantic_config_sha256"])
        self.write_bam(query="new_content")
        self.assertNotEqual(first["semantic_config_sha256"], self.preflight()["semantic_config_sha256"])


if __name__ == "__main__":
    unittest.main()
