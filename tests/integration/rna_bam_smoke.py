"""Run real BAM-import/StringTie/junction smoke checks in the all-in-one image.

python tests/integration/rna_bam_smoke.py --workdir /new/writable/directory
Requires the image's snakemake, stringtie, and pysam on PATH/PYTHONPATH.
"""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[2]
    (root / "genome.fa").write_text(">ctg1\n" + "acgt" * 500 + "\n")
    (root / "mapping.fa").write_text(">ctg1\n" + "ACGT" * 500 + "\n")
    (root / "repeats.gff3").write_text("ctg1\tEDTA\trepeat_region\t1\t4\t.\t+\t.\tID=r1\n")
    (root / "proteins.fa").write_text(">P1\nMPEPTIDE\n")
    metadata = [("PE1", "PE", "RF"), ("SE1", "SE", "R")]
    for sample, layout, strand in metadata:
        records = []
        for index in range(200):
            for group, start, cigar in (("left", 100 + index % 80, "120M"),
                                        ("junction", 200, "100M200N100M"),
                                        ("right", 540 + index % 100, "140M")):
                read = pysam.AlignedSegment()
                read.query_name = f"{sample}_{group}_{index}"
                read.flag = 83 if layout == "PE" else 16
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigarstring = cigar
                read.query_sequence = "A" * sum(length for op, length in read.cigartuples if op in (0, 1, 4, 7, 8))
                read.set_tag("XS", "+")
                read.set_tag("NH", 1)
                if layout == "PE":
                    read.next_reference_id = 0
                    read.next_reference_start = 650
                    mate = pysam.AlignedSegment()
                    mate.query_name = read.query_name
                    mate.flag = 163
                    mate.reference_id = 0
                    mate.reference_start = 650
                    mate.next_reference_id = 0
                    mate.next_reference_start = start
                    mate.mapping_quality = 60
                    mate.cigarstring = "140M"
                    mate.query_sequence = "A" * 140
                    mate.set_tag("XS", "+")
                    mate.set_tag("NH", 1)
                    read.template_length = 790 - start
                    mate.template_length = -(790 - start)
                    records.append(mate)
                records.append(read)
        bam = root / f"{sample}.bam"
        with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                 "SQ": [{"SN": "ctg1", "LN": 2000}]}) as handle:
            for read in sorted(records, key=lambda read: (read.reference_id, read.reference_start)):
                handle.write(read)
        pysam.index("-c", str(bam))
    with (root / "bam_samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "condition", "replicate", "layout", "strandedness", "bam", "csi"])
        for index, (sample, layout, strand) in enumerate(metadata, 1):
            writer.writerow([sample, "C", index, layout, strand, root / f"{sample}.bam", root / f"{sample}.bam.csi"])
    config = dict(assembly_id="RNA_SMOKE", organism="Test species", assembly_scope="nuclear",
                  reference_representation="haploid", expected_protein_copies=1,
                  genetic_code=1, gene_prefix="TESTG", genome=str(root / "genome.fa"),
                  repeat_gff=str(root / "repeats.gff3"), protein_fasta=str(root / "proteins.fa"),
                  samples=str(root / "bam_samples.tsv"), policy=str(source / "policies/production_v2.json"),
                  execution_id="bam", execution_root=str(root / "runs"),
                  rna={"input_mode": "bam", "reference_fasta": str(root / "mapping.fa")})
    cfg = root / "bam.json"
    cfg.write_text(json.dumps(config, indent=2))
    def run(label, configfile, *targets, dry=False):
        command = ["snakemake", "--snakefile", str(source / "Snakefile"), "--cores", "2",
                   "--configfile", str(configfile), "--printshellcmds"]
        if dry:
            command.append("--dry-run")
        command += list(targets)
        result = subprocess.run(command, cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (root / f"{label}.log").write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(f"{label} failed: {result.stdout[-6000:]}")
        return result.stdout
    original = {p: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
                for p in root.glob("*.bam*")}
    dry = run("bam_dryrun", cfg, "rna_evidence", "junction_evidence", dry=True)
    for forbidden in ("rule fastp_", "rule hisat2_", "rule alignment_versions"):
        assert forbidden not in dry, forbidden
    assert "rule import_rna_bam" in dry
    run("bam_execute", cfg, "rna_evidence", "junction_evidence")
    preflight = json.loads((root / "runs/bam/00_preflight/preflight.json").read_text())
    assert preflight["rna_input_mode"] == "bam" and preflight["sample_count"] == 2
    merged = root / "runs/bam/04_transcripts/merged/RNA_SMOKE.stringtie.merged.gtf"
    assert "\ttranscript\t" in merged.read_text()
    junctions = root / "runs/bam/06_evidence/01_junctions/aggregate/RNA_SMOKE.hq_junctions.tsv.gz"
    with gzip.open(junctions, "rt") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    assert len(lines) >= 1, lines
    for p, value in original.items():
        assert value == (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
    assert len(list((root / "runs/bam/03_alignment").rglob("*.import.json"))) == 2
    resume = run("bam_resume", cfg, "rna_evidence", "junction_evidence")
    assert "Nothing to be done" in resume
    # A full prediction dry run confirms the imported BAM producer is shared.
    dummy = root / "dummy.image"
    dummy.write_bytes(b"dry-run only; never execute predictor")
    prediction = {**config, "execution_id": "bam_prediction", "deployment": {
        "mode": "in_container", "predictor_sif": str(dummy),
        "predictor_sif_sha256": hashlib.sha256(dummy.read_bytes()).hexdigest(),
        "genemark_license_acknowledged": True}}
    prediction_cfg = root / "prediction.json"
    prediction_cfg.write_text(json.dumps(prediction))
    pred = run("bam_prediction_dryrun", prediction_cfg, "prediction", dry=True)
    assert "rule import_rna_bam" in pred and "rule hisat2_" not in pred
    # Original FASTQ mode still resolves its original trimming and mapping rules.
    for mate in (1, 2):
        (root / f"read{mate}.fq").write_text("@r\nACGT\n+\nFFFF\n")
    (root / "fastq_samples.tsv").write_text("sample\tcondition\treplicate\tlayout\tstrandedness\tread1\tread2\n"
        + f"S1\tC\t1\tPE\tRF\t{root / 'read1.fq'}\t{root / 'read2.fq'}\n")
    raw = {key: value for key, value in config.items() if key != "rna"}
    raw.update(samples=str(root / "fastq_samples.tsv"), execution_id="fastq")
    raw_cfg = root / "fastq.json"
    raw_cfg.write_text(json.dumps(raw))
    raw_dry = run("fastq_dryrun", raw_cfg, "rna_evidence", dry=True)
    assert "rule fastp_pe" in raw_dry and "rule hisat2_align_pe" in raw_dry
    assert "rule import_rna_bam" not in raw_dry
    report = {"status": "PASS", "real_stringtie_and_junction_execution": True,
              "bam_sample_count": 2, "layouts": ["PE", "SE"], "source_files_unchanged": True,
              "restart_noop": True, "bam_prediction_dryrun": True, "fastq_dryrun": True,
              "full_prediction_executed": False}
    (root / "complete.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
