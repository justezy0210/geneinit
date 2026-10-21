import csv
import json
import hashlib
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from geneinit.contracts.preflight import (
    PreflightError, run_input_preflight, run_preflight, input_preflight_config, _sha256,
)


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.genome = self.root / "genome.fa"
        self.repeats = self.root / "repeats.gff3"
        self.proteins = self.root / "proteins.faa"
        self.read1 = self.root / "sample_R1.fastq"
        self.read2 = self.root / "sample_R2.fastq"
        self.samples = self.root / "samples.tsv"
        self.policy = self.root / "policy.json"
        self.genome.write_text(">ctg1\nAAAAccccNNNN\n>ctg2\nACGT\n", encoding="utf-8")
        self.repeats.write_text(
            "##gff-version 3\nctg1\tEDTA\trepeat_region\t5\t8\t.\t+\t.\tID=r1\n",
            encoding="utf-8",
        )
        self.proteins.write_text(">p1\nMPEPTIDE*\n", encoding="utf-8")
        self.read1.write_text("@r1\nACGT\n+\nFFFF\n", encoding="utf-8")
        self.read2.write_text("@r1\nTGCA\n+\nFFFF\n", encoding="utf-8")
        self.samples.write_text(
            "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
            f"S1\tC1\t1\tPE\tRF\t{self.read1}\t{self.read2}\n",
            encoding="utf-8",
        )
        self.policy.write_text('{"policy_id":"test"}\n', encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def config(self):
        return {
            "assembly_id": "ASM1",
            "organism": "Genus species",
            "assembly_scope": "nuclear",
            "reference_representation": "haploid",
            "expected_protein_copies": 1,
            "genetic_code": 1,
            "gene_prefix": "ASM1G",
            "genome": str(self.genome),
            "repeat_gff": str(self.repeats),
            "protein_fasta": str(self.proteins),
            "samples": str(self.samples),
            "policy": str(self.policy),
            "execution_id": "first_root",
            "execution_root": str(self.root / "runs"),
        }

    def test_valid_inputs_write_complete_contract(self):
        output = self.root / "preflight"
        result = run_preflight(self.config(), output)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["sample_count"], 1)
        self.assertTrue((output / "preflight.complete.tsv").is_file())
        manifest = (output / "input_manifest.tsv").read_text(encoding="utf-8")
        self.assertIn("sample_R1.fastq", manifest)
        self.assertIn("sample_R2.fastq", manifest)

    def test_input_receipt_preserves_full_semantics_without_rehashing_sources(self):
        config = self.config()
        direct = run_preflight(config, self.root / "direct")
        inputs = self.root / "inputs"
        run_input_preflight(config, inputs)
        with patch("geneinit.contracts.preflight._sha256", wraps=_sha256) as hashed:
            cached = run_preflight(config, self.root / "cached",
                                   input_validation=inputs / "input_validation.json")
        self.assertEqual(direct, cached)
        self.assertEqual([call.args[0] for call in hashed.call_args_list], [self.policy])
        self.assertEqual(
            json.loads((self.root / "direct/normalized_config.json").read_text()),
            json.loads((self.root / "cached/normalized_config.json").read_text()),
        )

    def test_input_receipt_binds_relative_paths_to_the_original_directory(self):
        config = self.config()
        for key in ("genome", "repeat_gff", "protein_fasta", "samples"):
            config[key] = Path(config[key]).name
        inputs = self.root / "inputs"
        with patch("os.getcwd", return_value=str(self.root)):
            run_input_preflight(config, inputs)
        with patch("os.getcwd", return_value=str(self.root / "elsewhere")):
            with self.assertRaisesRegex(PreflightError, "configuration differs"):
                run_preflight(config, self.root / "moved_cwd",
                              input_validation=inputs / "input_validation.json")

    def test_downstream_settings_do_not_change_input_scope(self):
        original = self.config()
        changed = {**original, "policy": "elsewhere.json", "gene_prefix": "OTHER",
                   "expected_protein_copies": 7, "te_evidence": "te.json",
                   "evidence": {"junction_min_mapq": 60},
                   "prediction": {"augustus_rounds": 8}, "deployment": {"mode": "other"}}
        self.assertEqual(input_preflight_config(original), input_preflight_config(changed))
        changed["genome"] = "different.fa"
        self.assertNotEqual(input_preflight_config(original), input_preflight_config(changed))

    def test_reused_input_receipt_still_validates_downstream_fields(self):
        config = self.config()
        inputs = self.root / "inputs"
        run_input_preflight(config, inputs)
        for update, message in (({"expected_protein_copies": 0}, "positive integer"),
                                ({"gene_prefix": "!"}, "public ID"),
                                ({"policy": "absent.json"}, "readable file")):
            with self.subTest(update=update), self.assertRaisesRegex(PreflightError, message):
                run_preflight({**config, **update}, self.root / "bad_run",
                              input_validation=inputs / "input_validation.json")

    def test_reused_input_receipt_rejects_changed_source_or_corruption(self):
        config = self.config()
        inputs = self.root / "inputs"
        run_input_preflight(config, inputs)
        path = inputs / "input_validation.json"
        original = path.read_text()
        data = json.loads(original)
        data["payload"]["semantic_samples"][0]["replicate"] = "different"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(PreflightError, "invalid input-preflight receipt"):
            run_preflight(config, self.root / "bad_receipt", input_validation=path)
        path.write_text(original)
        self.read1.write_text("@different\nACGT\n+\nFFFF\n")
        with self.assertRaisesRegex(PreflightError, "source changed"):
            run_preflight(config, self.root / "changed_input", input_validation=path)

    def test_receipt_rejects_retargeted_read_symlink(self):
        alias = self.root / "read_alias.fastq"
        alias.symlink_to(self.read1)
        self.samples.write_text(self.samples.read_text().replace(str(self.read1), str(alias)))
        inputs = self.root / "inputs"
        run_input_preflight(self.config(), inputs)
        different = self.root / "different.fastq"
        different.write_text("@different\nAAAA\n+\nFFFF\n")
        alias.unlink()
        alias.symlink_to(different)
        with self.assertRaisesRegex(PreflightError, "sample input paths changed.*rerun input_preflight"):
            run_preflight(self.config(), self.root / "cached",
                          input_validation=inputs / "input_validation.json")
        self.assertFalse((self.root / "cached/preflight.complete.tsv").exists())

    def test_receipt_rejects_relative_sample_paths_from_different_directory(self):
        self.samples.write_text(self.samples.read_text().replace(str(self.root) + "/", ""))
        other = self.root / "other"
        other.mkdir()
        (other / self.read1.name).write_bytes(self.read1.read_bytes())
        (other / self.read2.name).write_bytes(self.read2.read_bytes())
        inputs = self.root / "inputs"
        with patch("os.getcwd", return_value=str(self.root)):
            run_input_preflight(self.config(), inputs)
        with patch("os.getcwd", return_value=str(other)):
            with self.assertRaisesRegex(PreflightError, "sample input paths changed.*rerun input_preflight"):
                run_preflight(self.config(), self.root / "cached",
                              input_validation=inputs / "input_validation.json")

    def test_receipt_revalidates_duplicate_read_paths(self):
        alias = self.root / "read_alias.fastq"
        alias.symlink_to(self.read1)
        self.samples.write_text(self.samples.read_text().replace(str(self.read1), str(alias)))
        inputs = self.root / "inputs"
        run_input_preflight(self.config(), inputs)
        alias.unlink()
        alias.symlink_to(self.read2)
        with self.assertRaisesRegex(PreflightError, "read1 and read2 are identical.*rerun input_preflight"):
            run_preflight(self.config(), self.root / "cached",
                          input_validation=inputs / "input_validation.json")

    def test_runtime_fields_do_not_change_semantic_hash(self):
        first = run_preflight(self.config(), self.root / "first")
        second_config = self.config()
        second_config["execution_id"] = "second_root"
        second_config["execution_root"] = str(self.root / "elsewhere")
        second = run_preflight(second_config, self.root / "second")
        self.assertEqual(
            first["semantic_config_sha256"], second["semantic_config_sha256"]
        )

    def test_explicit_fastq_mode_preserves_semantic_hash(self):
        first = run_preflight(self.config(), self.root / "implicit")
        config = self.config()
        config["rna"] = {"input_mode": "fastq"}
        second = run_preflight(config, self.root / "explicit")
        self.assertEqual(first["semantic_config_sha256"], second["semantic_config_sha256"])

    def test_masking_recovery_changes_biological_fingerprint(self):
        original = run_preflight(self.config(), self.root / "original")
        config = self.config()
        config["prediction"] = {"genemark_proc": "protein_evidence", "masking_rescue": False}
        disabled = run_preflight(config, self.root / "disabled")
        self.assertEqual(original["semantic_config_sha256"], disabled["semantic_config_sha256"])
        config["prediction"]["masking_rescue"] = True
        enabled = run_preflight(config, self.root / "enabled")
        self.assertNotEqual(original["semantic_config_sha256"], enabled["semantic_config_sha256"])

    def test_effective_prediction_and_evidence_defaults_are_equivalent(self):
        original = run_preflight(self.config(), self.root / "implicit_defaults")
        config = self.config()
        config["prediction"] = {
            "genemark_proc": "protein_evidence", "augustus_species": "ASM1_etp",
            "masking_rescue": False, "augustus_rounds": 5, "augustus_kfold": 8,
            "augustus_training_flank": 1000,
            "augustus_training": {"max_loci": 6000, "optimize_loci": 500,
                "holdout_loci": 500, "seed": 17, "identity": 80.0, "coverage": 80.0},
        }
        config["evidence"] = {"junction_min_mapq": 10, "junction_min_anchor": 8}
        explicit = run_preflight(config, self.root / "explicit_defaults")
        self.assertEqual(original["semantic_config_sha256"], explicit["semantic_config_sha256"])

    def test_each_prediction_and_evidence_setting_changes_fingerprint(self):
        original = run_preflight(self.config(), self.root / "baseline")
        changes = [
            ("prediction", "genemark_proc", "other"),
            ("prediction", "augustus_species", "other"),
            ("prediction", "augustus_rounds", 9),
            ("prediction", "augustus_kfold", 4),
            ("prediction", "augustus_training_flank", 2000),
            ("evidence", "junction_min_mapq", 60),
            ("evidence", "junction_min_anchor", 30),
        ]
        for group, key, value in changes:
            with self.subTest(setting=key):
                config = self.config()
                config[group] = {key: value}
                changed = run_preflight(config, self.root / key)
                self.assertNotEqual(original["semantic_config_sha256"], changed["semantic_config_sha256"])
        for key, value in {"max_loci": 7000, "optimize_loci": 600, "holdout_loci": 600,
                           "seed": 99, "identity": 90, "coverage": 90}.items():
            with self.subTest(training_setting=key):
                config = self.config()
                config["prediction"] = {"augustus_training": {key: value}}
                changed = run_preflight(config, self.root / key)
                self.assertNotEqual(original["semantic_config_sha256"], changed["semantic_config_sha256"])

    def test_prediction_artifact_content_is_locked_without_path_dependence(self):
        artifact = self.root / "parameters.cfg"
        artifact.write_text("original parameters")
        other = self.root / "copy.cfg"
        other.write_bytes(artifact.read_bytes())
        config = self.config()
        config["prediction"] = {"augustus_optimization_start_parameters": str(artifact)}
        first = run_preflight(config, self.root / "artifact_first")
        config["prediction"]["augustus_optimization_start_parameters"] = str(other)
        copied = run_preflight(config, self.root / "artifact_copy")
        self.assertEqual(first["semantic_config_sha256"], copied["semantic_config_sha256"])
        other.write_text("modified parameters")
        changed = run_preflight(config, self.root / "artifact_changed")
        self.assertNotEqual(first["semantic_config_sha256"], changed["semantic_config_sha256"])

    def test_copied_fastq_cannot_supply_an_independent_replicate(self):
        copy1, copy2 = self.root / "copy_R1.fastq", self.root / "copy_R2.fastq"
        copy1.write_bytes(self.read1.read_bytes())
        copy2.write_bytes(self.read2.read_bytes())
        with self.samples.open("a") as handle:
            handle.write(f"S2\tC1\t2\tPE\tRF\t{copy1}\t{copy2}\n")
        with self.assertRaisesRegex(PreflightError, "duplicate FASTQ content"):
            run_preflight(self.config(), self.root / "duplicate_content")

    def test_identical_content_of_distinct_mates_in_one_sample_is_allowed(self):
        self.read2.write_bytes(self.read1.read_bytes())
        result = run_preflight(self.config(), self.root / "identical_mate_content")
        self.assertEqual(result["sample_count"], 1)

    def test_repeat_gff_may_omit_repeat_free_contig(self):
        result = run_preflight(self.config(), self.root / "subset")
        self.assertEqual(result["repeat_gff"]["seqids"], ["ctg1"])

    def test_repeat_gff_unknown_seqid_fails(self):
        self.repeats.write_text(
            "unknown\tEDTA\trepeat_region\t1\t2\t.\t+\t.\tID=r1\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(PreflightError, "absent from genome"):
            run_preflight(self.config(), self.root / "bad_gff")

    def test_single_end_contract(self):
        self.samples.write_text(
            "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
            f"S1\tC1\t1\tSE\tR\t{self.read1}\t\n",
            encoding="utf-8",
        )
        result = run_preflight(self.config(), self.root / "single")
        self.assertEqual(result["layouts"], ["SE"])
        self.assertEqual(result["strandedness"], ["R"])

    def test_paired_end_requires_second_read(self):
        self.samples.write_text(
            "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
            f"S1\tC1\t1\tPE\tRF\t{self.read1}\t\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(PreflightError, "requires read2"):
            run_preflight(self.config(), self.root / "missing_read2")

    def test_nonstandard_genetic_code_fails_closed(self):
        config = self.config()
        config["genetic_code"] = 4
        with self.assertRaisesRegex(PreflightError, "genetic_code: 1 only"):
            run_preflight(config, self.root / "bad_code")

    def test_deployment_path_is_runtime_but_sif_digest_is_semantic(self):
        first_sif = self.root / "one.sif"
        second_sif = self.root / "two.sif"
        first_sif.write_bytes(b"same locked image")
        second_sif.write_bytes(b"same locked image")
        digest = hashlib.sha256(b"same locked image").hexdigest()

        first_config = self.config()
        first_config["deployment"] = {
            "predictor_sif": str(first_sif),
            "predictor_sif_sha256": digest,
            "genemark_license_acknowledged": True,
            "apptainer_executable": sys.executable,
        }
        second_config = self.config()
        second_config["deployment"] = {
            **first_config["deployment"],
            "predictor_sif": str(second_sif),
        }
        first = run_preflight(first_config, self.root / "deployment_one")
        second = run_preflight(second_config, self.root / "deployment_two")
        self.assertEqual(
            first["semantic_config_sha256"], second["semantic_config_sha256"]
        )


if __name__ == "__main__":
    unittest.main()
