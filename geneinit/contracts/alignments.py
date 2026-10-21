"""Validation and content-locked staging of precomputed RNA BAM + CSI inputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping

from .samples import SampleRecord


class AlignmentContractError(ValueError):
    """Raised when an imported alignment cannot be safely reused."""


def rna_input_config(config: Mapping[str, Any]) -> tuple[str, Path | None]:
    rna = config.get("rna", {})
    if not isinstance(rna, dict) or set(rna) - {"input_mode", "reference_fasta"}:
        raise AlignmentContractError("rna accepts only input_mode and reference_fasta")
    mode = rna.get("input_mode", "fastq")
    if mode not in ("fastq", "bam"):
        raise AlignmentContractError("rna.input_mode must be fastq or bam")
    if mode == "fastq":
        if "reference_fasta" in rna:
            raise AlignmentContractError("rna.reference_fasta is only valid in bam mode")
        return mode, None
    raw = rna.get("reference_fasta")
    if not isinstance(raw, str) or not raw.strip():
        raise AlignmentContractError("bam mode requires rna.reference_fasta used for alignment")
    reference = Path(raw).expanduser().resolve()
    if not reference.is_file() or not os.access(reference, os.R_OK) or not reference.stat().st_size:
        raise AlignmentContractError(f"rna.reference_fasta is not a readable nonempty file: {reference}")
    return mode, reference


def validate_alignment(sample: SampleRecord, lengths: Mapping[str, int]) -> dict[str, Any]:
    """Check header/EOF, CSI access, and a bounded sample of primary reads.

    This is not an exhaustive record scan or an inference of library strand.
    Library metadata and the alignment reference are supplied by the caller.
    """
    try:
        import pysam
    except ImportError as error:
        raise AlignmentContractError("BAM input validation requires pysam (envs/preflight.yaml)") from error
    if sample.bam is None or sample.csi is None:
        raise AlignmentContractError(f"sample {sample.sample} requires BAM and CSI")
    try:
        with gzip.open(sample.csi, "rb") as handle:
            if handle.read(4) != b"CSI\x01":
                raise AlignmentContractError("index must be CSI, not BAI")
        pysam.quickcheck(str(sample.bam))
        with pysam.AlignmentFile(str(sample.bam), "rb", index_filename=str(sample.csi)) as bam:
            if not bam.is_bam:
                raise AlignmentContractError("alignment must be BAM")
            if bam.header.to_dict().get("HD", {}).get("SO") != "coordinate":
                raise AlignmentContractError("BAM must declare coordinate sort order")
            if list(zip(bam.references, bam.lengths)) != list(lengths.items()):
                raise AlignmentContractError("BAM @SQ sequence IDs, order, or lengths differ from genome")
            if not bam.check_index():
                raise AlignmentContractError("BAM CSI index is unavailable")
            checked = 0
            anchor = None
            previous = (-1, -1)
            for index, read in enumerate(bam.fetch(until_eof=True)):
                if index >= 10000 or checked >= 256:
                    break
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                position = (read.reference_id, read.reference_start)
                if position < previous:
                    raise AlignmentContractError("sampled BAM records are not coordinate sorted")
                previous = position
                if read.is_paired != (sample.layout == "PE"):
                    raise AlignmentContractError(f"BAM read flags disagree with declared layout {sample.layout}")
                if anchor is None:
                    anchor = (read.reference_name, read.reference_start, read.query_name, read.flag)
                checked += 1
            if anchor is None:
                raise AlignmentContractError("no mapped primary reads in the first 10000 BAM records")
            name, start, query, flag = anchor
            if not any(read.query_name == query and read.flag == flag and read.reference_start == start
                       for read in bam.fetch(name, start, start + 1)):
                raise AlignmentContractError("CSI cannot retrieve a sampled BAM record")
            long_queries = 0
            for seqid, length in lengths.items():
                if length > 2**29:
                    # CSI must handle coordinates beyond the BAI limit. An empty
                    # interval is valid; coverage at that position is not required.
                    next(bam.fetch(seqid, 2**29, min(length, 2**29 + 10000)), None)
                    long_queries += 1
            return {
                "sample": sample.sample,
                "status": "PASS",
                "reference_count": len(lengths),
                "primary_reads_checked": checked,
                "csi_queries_beyond_512m": long_queries,
                "strandedness": sample.strandedness,
                "strandedness_source": "declared_library_metadata",
                "validation_scope": "header_eof_csi_and_sampled_records",
                "pysam_version": pysam.__version__,
                "alignment_programs": bam.header.to_dict().get("PG", []),
            }
    except AlignmentContractError as error:
        raise AlignmentContractError(f"sample {sample.sample}: {error}") from error
    except (OSError, ValueError, EOFError, pysam.SamtoolsError) as error:
        raise AlignmentContractError(f"sample {sample.sample}: invalid BAM/CSI: {error}") from error


def check_staging_paths(sample: SampleRecord, output_bam: Path, output_csi: Path,
                        audit: Path | None = None, *, other_sources=()) -> None:
    """Reject source/output aliasing before Snakemake can remove any output."""
    sources = {path.resolve() for path in (sample.bam, sample.csi, *other_sources) if path is not None}
    outputs = [output_bam, output_csi] + ([audit] if audit is not None else [])
    destinations = [path.parent.resolve() / path.name for path in outputs]
    if len(set(destinations)) != len(outputs):
        raise AlignmentContractError("RNA staging outputs must differ")
    for path in destinations:
        if path in sources or (path.exists() and not path.is_symlink()
                               and any(path.samefile(source) for source in sources)):
            raise AlignmentContractError(f"RNA staging output aliases an original input: {path}")


def stage_alignment(sample: SampleRecord, manifest: Path, output_bam: Path,
                    output_csi: Path, audit: Path) -> dict[str, Any]:
    """Verify preflight content locks, then symlink immutable source BAM + CSI."""
    check_staging_paths(sample, output_bam, output_csi, audit)
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["sample"] == sample.sample and row["role"] in {"bam", "csi"}]
    if len(rows) != 2 or {row["role"] for row in rows} != {"bam", "csi"}:
        raise AlignmentContractError(f"sample {sample.sample}: missing or duplicate preflight BAM/CSI locks")
    locked = {row["role"]: row for row in rows}
    for role, source in (("bam", sample.bam), ("csi", sample.csi)):
        if source is None or str(source.resolve()) != locked[role]["path"]:
            raise AlignmentContractError(f"sample {sample.sample}: {role} path differs from preflight lock")
        before = source.stat()
        digest = hashlib.sha256()
        with source.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
        after = source.stat()
        if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
            raise AlignmentContractError(f"sample {sample.sample}: {role} changed during import")
        if after.st_size != int(locked[role]["size_bytes"]) or digest.hexdigest() != locked[role]["sha256"]:
            raise AlignmentContractError(f"sample {sample.sample}: {role} content differs from preflight lock")
    # Validate all destinations before creating either symlink; never replace files.
    for source, output in ((sample.bam, output_bam), (sample.csi, output_csi)):
        if output.is_symlink() and output.resolve() == source.resolve():
            continue
        if output.exists() or output.is_symlink():
            raise AlignmentContractError(f"RNA staging output already exists: {output}")
    for source, output in ((sample.bam, output_bam), (sample.csi, output_csi)):
        output.parent.mkdir(parents=True, exist_ok=True)
        if not output.is_symlink():
            output.symlink_to(source.resolve())
    result = {
        "schema_version": 1, "status": "PASS", "input_mode": "bam",
        "sample": sample.sample, "condition": sample.condition, "replicate": sample.replicate,
        "layout": sample.layout, "strandedness": sample.strandedness,
        "staging": "symlink", "inputs": locked,
        "output_bam": str(output_bam.absolute()), "output_csi": str(output_csi.absolute()),
    }
    audit.parent.mkdir(parents=True, exist_ok=True)
    temporary = audit.with_suffix(audit.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(audit)
    return result
