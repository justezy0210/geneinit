"""Fail-closed validation and content locking for GeneInIt inputs."""

from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

from .samples import BAM_SAMPLE_COLUMNS, SAMPLE_COLUMNS, SampleContractError, load_samples
from .alignments import AlignmentContractError, rna_input_config, validate_alignment
from .deployment import DeploymentError, load_predictor_deployment


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
SAFE_PREFIX = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,19}$")
DNA_ALPHABET = set("ACGTRYSWKMBDHVNUacgtryswkmbdhvnu")
PROTEIN_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-")
REPRESENTATIONS = {
    "haploid",
    "primary_haploidized",
    "phased_diploid",
    "unphased_diploid",
    "polyploid",
}
RUNTIME_FIELDS = {
    "execution_id",
    "execution_root",
    "deployment.predictor_sif",
    "deployment.apptainer_executable",
}


class PreflightError(ValueError):
    """Raised when a production input violates the GeneInIt contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _require_file(raw: Any, label: str) -> Path:
    if not isinstance(raw, str) or not raw.strip():
        raise PreflightError(f"{label} must be a non-empty path")
    path = Path(raw).expanduser().resolve()
    if not path.is_file():
        raise PreflightError(f"{label} is not a readable file: {path}")
    if not os.access(path, os.R_OK):
        raise PreflightError(f"{label} is not readable: {path}")
    if path.stat().st_size == 0:
        raise PreflightError(f"{label} is empty: {path}")
    return path


def _parse_fasta(path: Path, alphabet: set[str], label: str) -> dict[str, Any]:
    lengths: dict[str, int] = {}
    current: str | None = None
    current_length = 0
    lowercase_bases = 0
    sequence_digest = hashlib.sha256()

    def finish_record() -> None:
        nonlocal current, current_length
        if current is None:
            return
        if current_length == 0:
            raise PreflightError(f"{label} contains an empty record: {current}")
        lengths[current] = current_length
        sequence_digest.update(b"\0")

    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                identifiers = line[1:].split(maxsplit=1)
                identifier = identifiers[0] if identifiers else ""
                if not identifier or not SAFE_ID.fullmatch(identifier):
                    raise PreflightError(
                        f"{label} has an unsafe FASTA ID at line {line_number}: {identifier!r}"
                    )
                if identifier in lengths or identifier == current:
                    raise PreflightError(f"{label} has duplicate FASTA ID: {identifier}")
                sequence_digest.update(identifier.encode("utf-8") + b"\0")
                current = identifier
                current_length = 0
                continue
            if current is None:
                raise PreflightError(f"{label} sequence appears before its header")
            invalid = set(line) - alphabet
            if invalid:
                raise PreflightError(
                    f"{label} record {current} has invalid characters: {sorted(invalid)!r}"
                )
            sequence_digest.update(line.upper().encode("ascii"))
            lowercase_bases += sum(character.islower() for character in line)
            current_length += len(line)
    finish_record()
    if not lengths:
        raise PreflightError(f"{label} has no FASTA records: {path}")
    return {
        "record_count": len(lengths),
        "total_bases": sum(lengths.values()),
        "lowercase_bases": lowercase_bases,
        "lengths": lengths,
        "sequence_sha256": sequence_digest.hexdigest(),
    }


def _validate_repeat_gff(path: Path, lengths: Mapping[str, int]) -> dict[str, Any]:
    feature_count = 0
    observed_seqids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise PreflightError(
                    f"repeat_gff line {line_number} does not have 9 columns"
                )
            seqid = fields[0]
            if seqid not in lengths:
                raise PreflightError(
                    f"repeat_gff seqid is absent from genome: {seqid}"
                )
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as error:
                raise PreflightError(
                    f"repeat_gff has non-integer coordinates at line {line_number}"
                ) from error
            if start < 1 or end < start or end > lengths[seqid]:
                raise PreflightError(
                    f"repeat_gff coordinates are out of bounds at line {line_number}: "
                    f"{seqid}:{start}-{end}"
                )
            observed_seqids.add(seqid)
            feature_count += 1
    return {
        "feature_count": feature_count,
        "seqid_count": len(observed_seqids),
        "seqids": sorted(observed_seqids),
    }


def _write_tsv(path: Path, header: Iterable[str], rows: Iterable[Iterable[Any]]) -> None:
    import csv

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(path)


def _write_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


INPUT_CONFIG_FIELDS = (
    "assembly_id", "assembly_scope", "genetic_code", "genome", "repeat_gff",
    "protein_fasta", "samples", "rna",
)


def input_preflight_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Only settings consumed while validating source data, before selection."""
    scoped = {key: config[key] for key in INPUT_CONFIG_FIELDS if key in config}
    for key in ("genome", "repeat_gff", "protein_fasta", "samples"):
        if isinstance(scoped.get(key), str):
            scoped[key] = str(Path(scoped[key]).expanduser().resolve())
    if isinstance(scoped.get("rna"), dict):
        scoped["rna"] = dict(scoped["rna"])
        reference = scoped["rna"].get("reference_fasta")
        if isinstance(reference, str):
            scoped["rna"]["reference_fasta"] = str(Path(reference).expanduser().resolve())
    return scoped


def _file_state(path: Path) -> list[int]:
    stat = path.stat()
    return [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_dev, stat.st_ino]


def _validate_input_data(config: Mapping[str, Any]) -> dict[str, Any]:
    if not SAFE_ID.fullmatch(str(config.get("assembly_id", ""))):
        raise PreflightError("assembly_id must match [A-Za-z0-9_.-]+")
    if config.get("assembly_scope") != "nuclear":
        raise PreflightError("GeneInIt v1 accepts assembly_scope: nuclear only")
    if config.get("genetic_code") != 1:
        raise PreflightError("GeneInIt currently supports nuclear genetic_code: 1 only")
    genome = _require_file(config.get("genome"), "genome")
    repeat_gff = _require_file(config.get("repeat_gff"), "repeat_gff")
    protein_fasta = _require_file(config.get("protein_fasta"), "protein_fasta")
    samples_path = _require_file(config.get("samples"), "samples")
    try:
        input_mode, alignment_reference = rna_input_config(config)
    except AlignmentContractError as error:
        raise PreflightError(str(error)) from error
    genome_stats = _parse_fasta(genome, DNA_ALPHABET, "genome")
    if alignment_reference is not None:
        reference_stats = (
            genome_stats if alignment_reference == genome else
            _parse_fasta(alignment_reference, DNA_ALPHABET, "rna.reference_fasta")
        )
        if reference_stats["sequence_sha256"] != genome_stats["sequence_sha256"]:
            raise PreflightError(
                "rna.reference_fasta differs from genome in sequence IDs, order, or bases; "
                "only soft-mask case and FASTA formatting differences are allowed"
            )
    repeat_stats = _validate_repeat_gff(repeat_gff, genome_stats["lengths"])
    protein_stats = _parse_fasta(protein_fasta, PROTEIN_ALPHABET, "protein_fasta")
    try:
        samples = load_samples(samples_path, input_mode=input_mode)
    except SampleContractError as error:
        raise PreflightError(str(error)) from error
    if repeat_stats["feature_count"] and genome_stats["lowercase_bases"] == 0:
        raise PreflightError(
            "repeat_gff has features but genome has no lowercase soft-masked bases"
        )

    manifest_rows: list[dict[str, Any]] = []

    def lock(role: str, path: Path, sample: str = "") -> dict[str, Any]:
        row = {
            "role": role,
            "sample": sample,
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        manifest_rows.append(row)
        return row

    locked = {
        "genome": lock("genome", genome),
        "repeat_gff": lock("repeat_gff", repeat_gff),
        "protein_fasta": lock("protein_fasta", protein_fasta),
        "samples": lock("samples", samples_path),
    }
    if alignment_reference is not None:
        locked["rna_reference_fasta"] = lock("rna_reference_fasta", alignment_reference)
    alignment_checks = []
    semantic_samples = []
    read_content_samples: dict[str, str] = {}
    for sample in samples:
        semantic_sample = {
            "sample": sample.sample,
            "condition": sample.condition,
            "replicate": sample.replicate,
            "layout": sample.layout,
            "strandedness": sample.strandedness,
        }
        if input_mode == "bam":
            try:
                alignment_checks.append(validate_alignment(sample, genome_stats["lengths"]))
            except AlignmentContractError as error:
                raise PreflightError(str(error)) from error
            semantic_sample.update({
                "bam_sha256": lock("bam", sample.bam, sample.sample)["sha256"],
                "csi_sha256": lock("csi", sample.csi, sample.sample)["sha256"],
            })
        else:
            read1 = lock("read1", sample.read1, sample.sample)
            read2 = lock("read2", sample.read2, sample.sample) if sample.read2 else None
            # Compare across samples only: two distinct mates in a small library
            # can legitimately have identical bytes.
            for read in (read1, read2):
                if read is None:
                    continue
                previous_sample = read_content_samples.get(read["sha256"])
                if previous_sample is not None and previous_sample != sample.sample:
                    raise PreflightError(
                        f"duplicate FASTQ content assigned to multiple samples: "
                        f"{previous_sample}, {sample.sample}"
                    )
                read_content_samples[read["sha256"]] = sample.sample
            semantic_sample.update({
                "read1_sha256": read1["sha256"],
                "read2_sha256": read2["sha256"] if read2 else None,
            })
        semantic_samples.append(semantic_sample)

    return {
        "input_config": input_preflight_config(config),
        "locked": locked, "manifest_rows": manifest_rows,
        "file_states": {row["path"]: _file_state(Path(row["path"])) for row in manifest_rows},
        "genome_stats": genome_stats, "repeat_stats": repeat_stats,
        "protein_stats": protein_stats, "input_mode": input_mode,
        "alignment_checks": alignment_checks, "semantic_samples": semantic_samples,
    }


def _write_samples(path: Path, samples, input_mode: str) -> None:
    _write_tsv(
        path, BAM_SAMPLE_COLUMNS if input_mode == "bam" else SAMPLE_COLUMNS,
        ((row.sample, row.condition, row.replicate, row.layout, row.strandedness,
          str(row.bam) if input_mode == "bam" else str(row.read1),
          str(row.csi) if input_mode == "bam" else (str(row.read2) if row.read2 else ""))
         for row in samples),
    )


def _write_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    _write_tsv(path, ("role", "sample", "size_bytes", "sha256", "path"),
               ((row["role"], row["sample"], row["size_bytes"], row["sha256"], row["path"])
                for row in rows))


def run_input_preflight(config: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    """Validate and hash source data once, independently of downstream policy."""
    output_dir.mkdir(parents=True, exist_ok=True)
    data = _validate_input_data(input_preflight_config(config))
    digest = _canonical_sha256(data)
    _write_json(output_dir / "input_validation.json", {
        "schema_version": 1, "status": "PASS", "payload_sha256": digest, "payload": data,
    })
    samples = load_samples(Path(data["locked"]["samples"]["path"]), input_mode=data["input_mode"])
    _write_samples(output_dir / "normalized_samples.tsv", samples, data["input_mode"])
    _write_manifest(output_dir / "input_manifest.tsv", data["manifest_rows"])
    _write_tsv(output_dir / "preflight.complete.tsv", ("schema_version", "status", "input_sha256"),
               ((1, "PASS", digest),))
    return {"status": "PASS", "input_sha256": digest, "sample_count": len(samples)}


def _load_input_validation(config: Mapping[str, Any], path: Path) -> dict[str, Any]:
    receipt = json.loads(path.read_text(encoding="utf-8"))
    data = receipt.get("payload")
    if (receipt.get("schema_version") != 1 or receipt.get("status") != "PASS"
            or not isinstance(data, dict)
            or receipt.get("payload_sha256") != _canonical_sha256(data)):
        raise PreflightError("invalid input-preflight receipt")
    if data.get("input_config") != input_preflight_config(config):
        raise PreflightError("input-preflight configuration differs; rerun input_preflight")
    for row in data["manifest_rows"]:
        source = Path(row["path"])
        if not source.is_file() or _file_state(source) != data["file_states"].get(row["path"]):
            raise PreflightError(f"source changed since input preflight: {source}")
    try:
        samples = load_samples(Path(data["locked"]["samples"]["path"]), input_mode=data["input_mode"])
    except SampleContractError as error:
        raise PreflightError(
            f"sample inputs changed since input preflight: {error}; rerun input_preflight"
        ) from error
    roles = ("bam", "csi") if data["input_mode"] == "bam" else ("read1", "read2")
    expected_paths = sorted(
        (row["role"], row["sample"], row["path"])
        for row in data["manifest_rows"] if row["role"] in roles
    )
    current_paths = sorted(
        (role, sample.sample, str(getattr(sample, role)))
        for sample in samples for role in roles if getattr(sample, role) is not None
    )
    if current_paths != expected_paths:
        raise PreflightError(
            "sample input paths changed since input preflight; rerun input_preflight"
        )
    return data


def run_preflight(config: Mapping[str, Any], output_dir: Path, *,
                  input_validation: Path | None = None) -> dict[str, Any]:
    """Validate inputs and write the stage-00 content lock.

    Runtime fields are recorded but excluded from the biological semantic hash.
    A verified input-preflight receipt avoids rehashing unchanged source data.
    """

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    required_scalars = (
        "assembly_id",
        "organism",
        "assembly_scope",
        "reference_representation",
        "expected_protein_copies",
        "genetic_code",
        "gene_prefix",
    )
    missing = [key for key in required_scalars if key not in config]
    if missing:
        raise PreflightError(f"missing config fields: {', '.join(missing)}")
    if not SAFE_ID.fullmatch(str(config["assembly_id"])):
        raise PreflightError("assembly_id must match [A-Za-z0-9_.-]+")
    if not str(config["organism"]).strip():
        raise PreflightError("organism must be non-empty")
    if config["assembly_scope"] != "nuclear":
        raise PreflightError("GeneInIt v1 accepts assembly_scope: nuclear only")
    if config["reference_representation"] not in REPRESENTATIONS:
        raise PreflightError("invalid reference_representation")
    if not isinstance(config["expected_protein_copies"], int) or config[
        "expected_protein_copies"
    ] < 1:
        raise PreflightError("expected_protein_copies must be a positive integer")
    if config["genetic_code"] != 1:
        raise PreflightError(
            "GeneInIt currently supports nuclear genetic_code: 1 only"
        )
    if not SAFE_PREFIX.fullmatch(str(config["gene_prefix"])):
        raise PreflightError("gene_prefix violates the public ID contract")

    data = (_load_input_validation(config, input_validation) if input_validation is not None
            else _validate_input_data(config))
    locked = dict(data["locked"])
    manifest_rows = list(data["manifest_rows"])
    genome, repeat_gff, protein_fasta, samples_path = (
        Path(locked[key]["path"]) for key in ("genome", "repeat_gff", "protein_fasta", "samples")
    )
    genome_stats, repeat_stats, protein_stats = (
        data[key] for key in ("genome_stats", "repeat_stats", "protein_stats")
    )
    input_mode = data["input_mode"]
    alignment_reference = (Path(locked["rna_reference_fasta"]["path"])
                           if "rna_reference_fasta" in locked else None)
    alignment_checks, semantic_samples = data["alignment_checks"], data["semantic_samples"]
    samples = load_samples(samples_path, input_mode=input_mode)

    def lock(role: str, path: Path, sample: str = "") -> dict[str, Any]:
        row = {"role": role, "sample": sample, "path": str(path),
               "size_bytes": path.stat().st_size, "sha256": _sha256(path)}
        manifest_rows.append(row)
        return row

    policy = _require_file(config.get("policy"), "policy")
    locked["policy"] = lock("policy", policy)
    te_evidence = _require_file(config["te_evidence"], "te_evidence") if config.get("te_evidence") else None
    if te_evidence is not None:
        from geneinit.evidence.te_classification import TeClassifier
        TeClassifier(te_evidence, str(config["assembly_id"]))
        locked["te_evidence"] = lock("te_evidence", te_evidence)
    deployment = None
    if "deployment" in config:
        try:
            deployment = load_predictor_deployment(config)
        except DeploymentError as error:
            raise PreflightError(str(error)) from error
        locked["predictor_sif"] = lock("predictor_sif", deployment.sif)

    biological_config = {
        key: config[key]
        for key in required_scalars
    }
    biological_config.update(
        {
            "genome_sha256": locked["genome"]["sha256"],
            "repeat_gff_sha256": locked["repeat_gff"]["sha256"],
            "protein_fasta_sha256": locked["protein_fasta"]["sha256"],
            "policy_sha256": locked["policy"]["sha256"],
            "samples": semantic_samples,
        }
    )
    if te_evidence is not None:
        biological_config["te_evidence_sha256"] = locked["te_evidence"]["sha256"]
    if input_mode == "bam":
        biological_config["rna_input_mode"] = "bam"
        biological_config["rna_reference_sequence_sha256"] = genome_stats["sequence_sha256"]
    if deployment is not None:
        biological_config["predictor_container_sha256"] = deployment.sif_sha256
    prediction = config.get("prediction", {})
    training = prediction.get("augustus_training", {})
    effective_prediction = {
        "genemark_proc": str(prediction.get("genemark_proc", "protein_evidence")),
        "augustus_species": str(prediction.get("augustus_species", f"{config['assembly_id']}_etp")),
        "masking_rescue": bool(prediction.get("masking_rescue", False)),
        "augustus_rounds": int(prediction.get("augustus_rounds", 5)),
        "augustus_kfold": int(prediction.get("augustus_kfold", 8)),
        "augustus_training_flank": int(prediction.get("augustus_training_flank", 1000)),
        "augustus_training": {
            key: int(training.get(key, default))
            for key, default in {
                "max_loci": 6000, "optimize_loci": 500, "holdout_loci": 500, "seed": 17,
            }.items()
        } | {
            key: float(training.get(key, 80)) for key in ("identity", "coverage")
        },
    }
    prediction_files = {
        key: prediction[key]
        for key in ("augustus_prediction_resume", "augustus_optimization_start_parameters")
        if prediction.get(key)
    }
    prediction_files.update({
        f"augustus_training_snapshot_{key}": value
        for key, value in prediction.get("augustus_training_snapshot", {}).items()
    })
    for key, value in prediction_files.items():
        path = _require_file(value, f"prediction.{key}")
        effective_prediction[f"{key}_sha256"] = lock(key, path)["sha256"]
    evidence = config.get("evidence", {})
    biological_config["prediction"] = effective_prediction
    biological_config["evidence"] = {
        "junction_min_mapq": int(evidence.get("junction_min_mapq", 10)),
        "junction_min_anchor": int(evidence.get("junction_min_anchor", 8)),
    }
    biological_config["semantic_config_schema_version"] = 2
    semantic_hash = _canonical_sha256(biological_config)

    normalized_config = dict(config)
    normalized_config.update(
        {
            "genome": str(genome),
            "repeat_gff": str(repeat_gff),
            "protein_fasta": str(protein_fasta),
            "samples": str(samples_path),
            "policy": str(policy),
            "semantic_config_sha256": semantic_hash,
            "semantic_config_schema_version": 2,
            "runtime_fields_excluded_from_semantic_hash": sorted(RUNTIME_FIELDS),
        }
    )
    if deployment is not None:
        normalized_deployment = {
            **dict(config["deployment"]),
            "mode": deployment.mode,
            "predictor_sif": str(deployment.sif),
            "predictor_sif_sha256": deployment.sif_sha256,
        }
        if deployment.apptainer is not None:
            normalized_deployment["apptainer_executable"] = str(deployment.apptainer)
        else:
            normalized_deployment.pop("apptainer_executable", None)
        normalized_config["deployment"] = normalized_deployment
    normalized_config["rna"] = {"input_mode": input_mode}
    if alignment_reference is not None:
        normalized_config["rna"]["reference_fasta"] = str(alignment_reference)
    summary = {
        "schema_version": 1,
        "status": "PASS",
        "semantic_config_sha256": semantic_hash,
        "semantic_config_schema_version": 2,
        "genome": {key: value for key, value in genome_stats.items() if key != "lengths"},
        "repeat_gff": repeat_stats,
        "protein_fasta": {
            key: value for key, value in protein_stats.items() if key != "lengths"
        },
        "sample_count": len(samples),
        "rna_input_mode": input_mode,
        "alignment_checks": alignment_checks,
        "layouts": sorted({row.layout for row in samples}),
        "strandedness": sorted({row.strandedness for row in samples}),
        "predictor_deployment_locked": deployment is not None,
    }

    _write_json(output_dir / "normalized_config.json", normalized_config)
    _write_json(output_dir / "preflight.json", summary)
    _write_samples(output_dir / "normalized_samples.tsv", samples, input_mode)
    _write_manifest(output_dir / "input_manifest.tsv", manifest_rows)
    _write_tsv(
        output_dir / "preflight.complete.tsv",
        ("schema_version", "status", "semantic_config_sha256"),
        ((1, "PASS", semantic_hash),),
    )
    return summary
