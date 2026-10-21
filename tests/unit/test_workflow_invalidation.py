"""Check the real DAG's reuse boundaries without running external predictors."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest


@pytest.mark.skipif(importlib.util.find_spec("snakemake") is None, reason="Snakemake is not installed")
def test_stage_specific_workflow_invalidation(tmp_path):
    original_source = Path(__file__).resolve().parents[2]
    source = tmp_path / "source"
    shutil.copytree(original_source, source, ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache"))
    fixture = source / "tests/fixtures/preflight"
    for name in ("genome.fa", "proteins.faa", "repeats.gff3", "sample_R1.fastq", "sample_R2.fastq"):
        shutil.copyfile(fixture / name, tmp_path / name)
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
        f"S1\tC\t1\tPE\tRF\t{tmp_path / 'sample_R1.fastq'}\t{tmp_path / 'sample_R2.fastq'}\n"
    )
    policy = tmp_path / "policy.json"
    shutil.copyfile(source / "policies/unified_v10_single_exon50_corroborated_multisample.json", policy)
    sif = tmp_path / "predictor.sif"
    sif.write_text("scheduler fixture only; no predictor execution")
    te = tmp_path / "te.json"
    te.write_text('{"schema_version": 1}\n')
    config = {
        "assembly_id": "TEST1", "organism": "Test species", "assembly_scope": "nuclear",
        "reference_representation": "haploid", "expected_protein_copies": 1,
        "genetic_code": 1, "gene_prefix": "TEST1G", "genome": str(tmp_path / "genome.fa"),
        "repeat_gff": str(tmp_path / "repeats.gff3"), "protein_fasta": str(tmp_path / "proteins.faa"),
        "samples": str(samples), "policy": str(policy),
        "execution_root": str(tmp_path / "runs"), "execution_id": "test",
        "deployment": {"mode": "in_container", "predictor_sif": str(sif),
                       "predictor_sif_sha256": hashlib.sha256(sif.read_bytes()).hexdigest(),
                       "genemark_license_acknowledged": True},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    env = dict(os.environ, PYTHONPATH=str(source) + os.pathsep + os.environ.get("PYTHONPATH", ""),
               PYTHONDONTWRITEBYTECODE="1")
    command = [sys.executable, "-m", "snakemake", "all", "--snakefile", str(source / "Snakefile"),
               "--cores", "1", "--runtime-source-cache-path", str(tmp_path / "cache"),
               "--configfile", str(config_path)]

    def run(label, *options, target="all", success=True):
        invocation = list(command)
        invocation[3] = target
        result = subprocess.run([*invocation, *options], cwd=tmp_path, env=env,
                                capture_output=True, text=True)
        output = result.stdout + result.stderr
        (tmp_path / f"{label}.log").write_text(output)
        assert (result.returncode == 0) == success, output
        return output

    # Run only the actual small validation rules, then materialize inert
    # analysis placeholders and let Snakemake record actual rule metadata.
    run("validate_initial_inputs", target="preflight")
    # Materialize inert placeholders, then let Snakemake record actual rule
    # metadata via --touch. No analysis command is executed.
    summary = run("initial_summary", "--summary")
    outputs = [line.split("\t", 1)[0] for line in summary.splitlines()
               if line.startswith(str(tmp_path / "runs") + "/")]
    assert outputs
    for name in outputs:
        output = Path(name)
        if output.name in {"work", "config", "trained_species"}:
            output.mkdir(parents=True, exist_ok=True)
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.touch()
    run("record_rule_metadata", "--touch", "--forceall")

    def planned(label):
        output = run(label, "--dry-run")
        return set(re.findall(r"(?:local)?rule ([A-Za-z0-9_]+):", output))

    assert planned("unchanged") == set()
    upstream = {"input_preflight", "fastp_pe", "hisat2_build", "hisat2_align_pe",
                "stringtie_sample", "stringtie_merge", "predictor_container_preflight",
                "genemark_etp", "build_full_etp_hints", "prepare_augustus_config", "augustus_direct"}
    changes = [
        ("copies", {"expected_protein_copies": 2}, {"preflight", "score_component_evidence"}),
        ("prefix", {"gene_prefix": "OTHER"}, {"preflight", "build_structural_release"}),
        ("te", {"te_evidence": str(te)}, {"preflight", "select_coding_models"}),
        ("junction_threshold", {"evidence": {"junction_min_mapq": 20}},
         {"preflight", "extract_sample_junctions", "aggregate_junctions"}),
    ]
    for label, update, expected in changes:
        config_path.write_text(json.dumps({**config, **update}))
        dirty = planned(label)
        assert expected <= dirty, (label, dirty)
        assert not dirty & upstream, (label, dirty & upstream)
        if label in {"prefix", "te"}:
            assert "score_component_evidence" not in dirty, (label, dirty)
    config_path.write_text(json.dumps(config))

    def changed_file(path, label):
        original, stat = path.read_bytes(), path.stat()
        path.write_bytes(original + b"\n")
        # File contents as well as mtime change, so this is not merely --touch.
        future = max(p.stat().st_mtime_ns for p in map(Path, outputs)) + 2_000_000_000
        os.utime(path, ns=(future, future))
        try:
            return planned(label)
        finally:
            path.write_bytes(original)
            os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    dirty = changed_file(policy, "policy_content")
    assert {"preflight", "select_coding_models", "run_structural_postqc"} <= dirty
    assert not dirty & upstream
    assert "score_component_evidence" not in dirty
    dirty = changed_file(tmp_path / "sample_R1.fastq", "read_content")
    assert {"input_preflight", "preflight", "fastp_pe", "hisat2_align_pe",
            "genemark_etp", "augustus_direct", "select_coding_models"} <= dirty
    config_path.write_text(json.dumps({**config, "prediction": {"genemark_proc": "protein_evidence", "augustus_rounds": 7}}))
    dirty = planned("training_setting")
    assert {"preflight", "augustus_direct"} <= dirty
    assert not dirty & {"input_preflight", "fastp_pe", "hisat2_build", "hisat2_align_pe", "genemark_etp"}
    config_path.write_text(json.dumps({**config, "deployment": {
        **config["deployment"], "predictor_sif_sha256": "0" * 64}}))
    dirty = planned("deployment_setting")
    assert {"preflight", "predictor_container_preflight", "genemark_etp", "augustus_direct"} <= dirty
    assert not dirty & {"input_preflight", "fastp_pe", "hisat2_build", "hisat2_align_pe"}

    config_path.write_text(json.dumps(config))
    code_cases = [
        ("rna_code", "geneinit/evidence/rna_cds.py", {"score_component_evidence", "run_structural_postqc"}),
        ("selector_code", "geneinit/selection/unified.py", {"select_coding_models", "run_structural_postqc"}),
        ("postqc_code", "geneinit/postqc/unified.py", {"run_structural_postqc"}),
        ("release_code", "geneinit/release/core.py", {"build_structural_release"}),
        ("predictor_code", "geneinit/predictor/commands.py", {"genemark_etp", "augustus_direct"}),
        ("training_code", "geneinit/predictor/training.py", {"augustus_direct"}),
        ("streaming_code", "geneinit/contracts/streaming.py", {"score_component_evidence", "run_structural_postqc"}),
    ]
    for label, relative, expected in code_cases:
        dirty = changed_file(source / relative, label)
        assert expected <= dirty, (label, dirty)
        assert not dirty & {"input_preflight", "fastp_pe", "hisat2_build", "hisat2_align_pe"}, (label, dirty)
        if label in {"selector_code", "postqc_code", "release_code"}:
            assert "score_component_evidence" not in dirty, (label, dirty)
        if label == "release_code":
            assert "select_coding_models" not in dirty and "run_structural_postqc" not in dirty

    # Fixture edits deliberately changed ctime too. Refresh its input receipt
    # before exercising a downstream-only validation failure.
    run("refresh_fixture_receipt", "--forcerun", "input_preflight", target="input_preflight")
    run("refresh_rule_metadata", "--touch", "--forceall")
    assert planned("restored_noop") == set()
    # A bad downstream contract must fail before either selection or release
    # can consume the otherwise complete upstream placeholders.
    config_path.write_text(json.dumps({**config, "te_evidence": str(te)}))
    protected = [tmp_path / "runs/test/07_selection/TEST1.selection_audit.tsv",
                 tmp_path / "runs/test/11_release.staging/TEST1/annotation.gff3"]
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in protected}
    output = run("failed_validation_blocks_consumers", success=False)
    assert "Error in rule preflight" in output
    assert "unsupported TE evidence schema" in output
    assert "localrule select_coding_models:" not in output
    assert "localrule build_structural_release:" not in output
    assert before == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in protected}
