"""RNA sample manifest and strandedness contracts."""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path


SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
SAMPLE_COLUMNS = (
    "sample",
    "condition",
    "replicate",
    "layout",
    "strandedness",
    "read1",
    "read2",
)
BAM_SAMPLE_COLUMNS = SAMPLE_COLUMNS[:5] + ("bam", "csi")
PROTOCOL_ARGUMENTS = {
    ("PE", "RF"): ("RF", "--rf"),
    ("PE", "FR"): ("FR", "--fr"),
    ("PE", "unstranded"): (None, None),
    ("SE", "R"): ("R", "--rf"),
    ("SE", "F"): ("F", "--fr"),
    ("SE", "unstranded"): (None, None),
}


class SampleContractError(ValueError):
    """Raised when the RNA manifest violates the declared protocol."""


@dataclass(frozen=True, slots=True)
class SampleRecord:
    sample: str
    condition: str
    replicate: str
    layout: str
    strandedness: str
    read1: Path | None = None
    read2: Path | None = None
    bam: Path | None = None
    csi: Path | None = None

    @property
    def hisat2_strandness(self) -> str | None:
        return PROTOCOL_ARGUMENTS[(self.layout, self.strandedness)][0]

    @property
    def stringtie_flag(self) -> str | None:
        return PROTOCOL_ARGUMENTS[(self.layout, self.strandedness)][1]

    @property
    def ingenannot_paired(self) -> bool:
        return self.layout == "PE"

    @property
    def ingenannot_stranded(self) -> bool:
        # InGenAnnot 0.0.20 implements reverse-stranded read filtering only.
        # Forward-stranded libraries therefore use its orientation-neutral
        # counter; strand-aware transcript assembly remains upstream.
        return self.strandedness in {"RF", "R"}

    @property
    def ingenannot_protocol_mode(self) -> str:
        if self.strandedness in {"RF", "R"}:
            return "REVERSE_STRANDED_EXACT"
        if self.strandedness == "unstranded":
            return "UNSTRANDED_EXACT"
        return "FORWARD_STRANDED_ORIENTATION_NEUTRAL_FALLBACK"


def _read_path(raw: str, label: str, require_files: bool) -> Path:
    path = Path(raw).expanduser().resolve()
    if require_files:
        if not path.is_file() or not os.access(path, os.R_OK):
            raise SampleContractError(f"{label} is not a readable file: {path}")
        if path.stat().st_size == 0:
            raise SampleContractError(f"{label} is empty: {path}")
    return path


def load_samples(path: Path, *, input_mode: str = "fastq", require_files: bool = True) -> list[SampleRecord]:
    """Load, validate, normalize, and sample-ID-sort one RNA manifest."""

    if input_mode not in {"fastq", "bam"}:
        raise SampleContractError("rna.input_mode must be fastq or bam")
    columns = BAM_SAMPLE_COLUMNS if input_mode == "bam" else SAMPLE_COLUMNS
    seen_bams: set[Path] = set()
    seen_reads: dict[Path, str] = {}
    path = path.expanduser().resolve()
    rows: list[SampleRecord] = []
    seen_samples: set[str] = set()
    seen_replicates: set[tuple[str, str]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(columns):
            raise SampleContractError(
                "samples header must be exactly: " + "\t".join(columns)
            )
        for line_number, source_row in enumerate(reader, 2):
            row = {key: (source_row.get(key) or "").strip() for key in columns}
            sample = row["sample"]
            if not SAFE_SAMPLE_ID.fullmatch(sample):
                raise SampleContractError(
                    f"samples line {line_number} has unsafe sample ID: {sample!r}"
                )
            if sample in seen_samples:
                raise SampleContractError(f"duplicate sample ID: {sample}")
            if not row["condition"] or not row["replicate"]:
                raise SampleContractError(
                    f"sample {sample} requires condition and replicate"
                )
            replicate_key = (row["condition"], row["replicate"])
            if replicate_key in seen_replicates:
                raise SampleContractError(
                    f"duplicate condition/replicate: {replicate_key[0]}/{replicate_key[1]}"
                )
            protocol = (row["layout"], row["strandedness"])
            if protocol not in PROTOCOL_ARGUMENTS:
                allowed = "PE: RF|FR|unstranded; SE: R|F|unstranded"
                raise SampleContractError(
                    f"sample {sample} has invalid protocol {protocol!r}; allowed {allowed}"
                )
            read1 = read2 = bam = csi = None
            if input_mode == "bam":
                for field in ("bam", "csi"):
                    if not row[field]:
                        raise SampleContractError(f"sample {sample} requires {field}")
                bam = _read_path(row["bam"], f"sample {sample} bam", require_files)
                csi = _read_path(row["csi"], f"sample {sample} csi", require_files)
                if bam == csi:
                    raise SampleContractError(f"sample {sample} bam and csi are identical")
                if bam in seen_bams or (require_files and any(bam.samefile(other) for other in seen_bams)):
                    raise SampleContractError(f"duplicate BAM assigned to multiple samples: {bam}")
                seen_bams.add(bam)
            else:
                if not row["read1"]:
                    raise SampleContractError(f"sample {sample} requires read1")
                if row["layout"] == "PE" and not row["read2"]:
                    raise SampleContractError(f"PE sample {sample} requires read2")
                if row["layout"] == "SE" and row["read2"]:
                    raise SampleContractError(f"SE sample {sample} must not define read2")
                read1 = _read_path(row["read1"], f"sample {sample} read1", require_files)
                read2 = (
                    _read_path(row["read2"], f"sample {sample} read2", require_files)
                    if row["read2"]
                    else None
                )
                if read2 is not None and (read1 == read2 or (
                    require_files and read1.samefile(read2)
                )):
                    raise SampleContractError(f"sample {sample} read1 and read2 are identical")
                for read in (read1, read2):
                    if read is None:
                        continue
                    for previous, previous_sample in seen_reads.items():
                        if read == previous or (require_files and read.samefile(previous)):
                            raise SampleContractError(
                                f"duplicate FASTQ assigned to multiple samples: "
                                f"{previous_sample}, {sample}: {read}"
                            )
                    seen_reads[read] = sample
            rows.append(
                SampleRecord(
                    sample=sample,
                    condition=row["condition"],
                    replicate=row["replicate"],
                    layout=row["layout"],
                    strandedness=row["strandedness"],
                    read1=read1,
                    read2=read2,
                    bam=bam,
                    csi=csi,
                )
            )
            seen_samples.add(sample)
            seen_replicates.add(replicate_key)
    if not rows:
        raise SampleContractError("samples table has no data rows")
    return sorted(rows, key=lambda row: row.sample)
