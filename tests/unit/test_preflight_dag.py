"""Exercise config-only invalidation through the actual Snakemake scheduler."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.skipif(importlib.util.find_spec("snakemake") is None, reason="Snakemake is not installed")
def test_preflight_reruns_when_only_config_changes(tmp_path):
    source = Path(__file__).resolve().parents[2]
    fixture = source / "tests/fixtures/preflight"
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
        f"S1\tC\t1\tPE\tRF\t{fixture / 'sample_R1.fastq'}\t{fixture / 'sample_R2.fastq'}\n"
    )
    config = {
        "assembly_id": "TEST1", "organism": "Test species", "assembly_scope": "nuclear",
        "reference_representation": "haploid", "expected_protein_copies": 1,
        "genetic_code": 1, "gene_prefix": "TEST1G", "genome": str(fixture / "genome.fa"),
        "repeat_gff": str(fixture / "repeats.gff3"), "protein_fasta": str(fixture / "proteins.faa"),
        "samples": str(samples), "policy": str(source / "policies/production_v2.json"),
        "execution_root": str(tmp_path / "runs"), "execution_id": "test",
    }
    config_path = tmp_path / "config.json"
    env = dict(os.environ, PYTHONPATH=str(source) + os.pathsep + os.environ.get("PYTHONPATH", ""),
               PYTHONDONTWRITEBYTECODE="1")
    for copies in (1, 2):
        config["expected_protein_copies"] = copies
        config_path.write_text(json.dumps(config))
        run = subprocess.run(
            [sys.executable, "-m", "snakemake", "preflight", "--snakefile", str(source / "Snakefile"),
             "--cores", "1", "--runtime-source-cache-path", str(tmp_path / "cache"),
             "--configfile", str(config_path)],
            cwd=tmp_path, env=env, capture_output=True, text=True,
        )
        assert run.returncode == 0, run.stdout + run.stderr
        recorded = json.loads((tmp_path / "runs/test/00_preflight/normalized_config.json").read_text())
        assert recorded["expected_protein_copies"] == copies
