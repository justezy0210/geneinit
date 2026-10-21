"""Prepare, normalize, and aggregate fail-closed final InGenAnnot UTR shards.

InGenAnnot is used only to rank RNA isoforms and propose transcript boundaries.
This module remains the authority for model identity: every active post-QC CDS
segment and phase must survive the single final UTR refinement unchanged.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from dataclasses import dataclass, field
import hashlib
from itertools import zip_longest
import os
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO
from urllib.parse import quote, unquote

from geneinit.contracts.streaming import DiskBackedUniqueIds
from geneinit.models.candidate import CdsPart


GFF3_SAFE = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:^*$@!+_?-|"
EMPTY = {"", ".", "NA", "NONE"}
PREPARATION_COLUMNS = (
    "candidate_id", "temporary_gene_id", "temporary_transcript_id", "seqid",
    "strand", "start", "end", "selection_locus_id", "predictors",
    "source_gene_ids", "source_transcript_ids", "routes", "active_fusion",
)
POSTQC_PREPARATION_COLUMNS = PREPARATION_COLUMNS + (
    "postqc_gene_id", "source_status",
)
UTR_MAPPING_COLUMNS = PREPARATION_COLUMNS + (
    "utr_refine_evidence", "five_prime_utr_segments", "five_prime_utr_bp",
    "three_prime_utr_segments", "three_prime_utr_bp",
)


class UtrError(RuntimeError):
    """Raised when RNA refinement changes a selected coding structure."""


@dataclass(frozen=True, slots=True)
class Feature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: str
    attributes: Mapping[str, str]

    @property
    def signature(self) -> tuple[int, int, str]:
        return self.start, self.end, self.phase


@dataclass(slots=True)
class TranscriptModel:
    gene_id: str
    transcript_id: str
    candidate_id: str
    seqid: str
    strand: str
    source: str
    exons: list[Feature] = field(default_factory=list)
    cds: list[Feature] = field(default_factory=list)
    five_prime_utrs: list[Feature] = field(default_factory=list)
    three_prime_utrs: list[Feature] = field(default_factory=list)
    utr_refine_evidence: str = ""

    @property
    def start(self) -> int:
        return min(feature.start for feature in self.exons)

    @property
    def end(self) -> int:
        return max(feature.end for feature in self.exons)


@contextmanager
def atomic_text(path: Path, *, newline: str | None = None) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline=newline) as handle:
            yield handle
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for raw in text.strip().strip(";").split(";"):
        token = raw.strip()
        if not token:
            continue
        fields = token.split(None, 1)
        if len(fields) != 2:
            raise UtrError(f"malformed GTF attribute: {token!r}")
        key, value = fields
        value = value.strip().strip('"')
        if key in attributes:
            raise UtrError(f"duplicate GTF attribute: {key}")
        attributes[key] = value
    return attributes


def parse_gff3_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for token in text.strip().strip(";").split(";"):
        if not token:
            continue
        if "=" not in token:
            raise UtrError(f"malformed GFF3 attribute: {token!r}")
        key, value = token.split("=", 1)
        if not key or key in attributes:
            raise UtrError(f"empty/duplicate GFF3 attribute: {key!r}")
        attributes[key] = unquote(value)
    return attributes


def format_gff3_attributes(values: Iterable[tuple[str, str]]) -> str:
    return ";".join(f"{key}={quote(str(value), safe=GFF3_SAFE)}" for key, value in values) + ";"


def parse_cds_signature(text: str) -> tuple[CdsPart, ...]:
    parts: list[CdsPart] = []
    for token in text.split(","):
        try:
            coordinate, phase = token.rsplit(":", 1)
            start, end = coordinate.split("-", 1)
            part = CdsPart(int(start), int(end), int(phase))
        except (TypeError, ValueError) as error:
            raise UtrError(f"invalid CDS signature: {text!r}") from error
        if part.start < 1 or part.end < part.start:
            raise UtrError(f"invalid CDS interval: {token!r}")
        parts.append(part)
    if not parts:
        raise UtrError("empty CDS signature")
    return tuple(sorted(parts))


def _read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise UtrError(f"{path} lacks columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def _iter_table(
    path: Path, required: set[str]
) -> Iterator[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise UtrError(f"{path} lacks columns: {sorted(missing)}")
        for row in reader:
            yield {
                key: (value or "").strip() for key, value in row.items()
            }


def _retained_candidate_metadata(
    candidate_table: Path,
    selection_audit: Path,
    retained_ids: set[str] | None,
    *,
    require_structure: bool,
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Lockstep-stream candidate/audit tables and retain only needed rows."""

    candidate_required = {
        "candidate_id", "predictors", "source_gene_ids", "source_transcript_ids"
    }
    audit_required = {"candidate_id", "active_fusion"}
    if require_structure:
        candidate_required.update(
            {"seqid", "strand", "start", "end", "cds_signature"}
        )
        audit_required.update({"locus_id", "selected", "routes"})
    candidates = _iter_table(candidate_table, candidate_required)
    audits = _iter_table(selection_audit, audit_required)
    candidate_metadata: dict[str, dict[str, str]] = {}
    audit_metadata: dict[str, dict[str, str]] = {}
    sentinel = object()
    previous_key: tuple[object, ...] | None = None
    with DiskBackedUniqueIds() as unique_ids:
        for number, pair in enumerate(
            zip_longest(candidates, audits, fillvalue=sentinel), 2
        ):
            candidate, audit = pair
            if sentinel in pair:
                raise UtrError("candidate and selection-audit row counts differ")
            assert isinstance(candidate, dict)
            assert isinstance(audit, dict)
            candidate_id = candidate["candidate_id"]
            if candidate_id != audit["candidate_id"]:
                raise UtrError(
                    f"candidate/audit order differs at logical row {number}: "
                    f"{candidate_id}, {audit['candidate_id']}"
                )
            if require_structure:
                key = (
                    candidate["seqid"], int(candidate["start"]),
                    int(candidate["end"]), candidate["strand"], candidate_id,
                )
                if previous_key is not None and key < previous_key:
                    raise UtrError(
                        "candidate/audit tables are not in coordinate order"
                    )
                previous_key = key
            if not candidate_id or not unique_ids.add(candidate_id):
                raise UtrError(
                    f"duplicate/empty candidate ID: {candidate_id!r}"
                )
            if retained_ids is None:
                retain = audit["selected"] == "1"
            else:
                retain = candidate_id in retained_ids
            if retain:
                candidate_metadata[candidate_id] = candidate
                audit_metadata[candidate_id] = audit
    if retained_ids is not None:
        missing = retained_ids - set(candidate_metadata)
        if missing:
            raise UtrError(
                "required candidates are absent from candidate/audit inputs: "
                f"{sorted(missing)[:10]}"
            )
    return candidate_metadata, audit_metadata


def _write_table(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with atomic_text(path, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _gtf_attributes(values: Iterable[tuple[str, str]]) -> str:
    return " ".join(f'{key} "{value}";' for key, value in values)


def _parse_feature(raw: str, fmt: str, path: Path, number: int) -> Feature:
    fields = raw.rstrip("\n").split("\t")
    if len(fields) != 9:
        raise UtrError(f"{path}:{number} does not have 9 columns")
    try:
        start, end = int(fields[3]), int(fields[4])
    except ValueError as error:
        raise UtrError(f"{path}:{number} has non-integer coordinates") from error
    if start < 1 or end < start or fields[6] not in {"+", "-"}:
        raise UtrError(f"{path}:{number} has invalid coordinate/strand")
    attributes = parse_gtf_attributes(fields[8]) if fmt == "gtf" else parse_gff3_attributes(fields[8])
    return Feature(fields[0], fields[1], fields[2], start, end, fields[6], fields[7], attributes)


def _features(path: Path, fmt: str) -> list[Feature]:
    result: list[Feature] = []
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if raw.strip() and not raw.startswith("#"):
                result.append(_parse_feature(raw, fmt, path, number))
    return result


def load_coding_gtf(path: Path) -> tuple[list[str], dict[str, TranscriptModel]]:
    features = _features(path, "gtf")
    genes: dict[str, Feature] = {}
    order: list[str] = []
    models: dict[str, TranscriptModel] = {}
    for feature in features:
        attrs = feature.attributes
        if feature.feature_type == "gene":
            gene_id = attrs.get("gene_id", "")
            if not gene_id or gene_id in genes:
                raise UtrError(f"missing/duplicate coding gene ID: {gene_id!r}")
            genes[gene_id] = feature
        elif feature.feature_type in {"transcript", "mRNA"}:
            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            candidate_id = attrs.get("candidate_id", "")
            if gene_id not in genes or not transcript_id or not candidate_id or transcript_id in models:
                raise UtrError(f"invalid coding transcript hierarchy: {transcript_id!r}")
            models[transcript_id] = TranscriptModel(
                gene_id, transcript_id, candidate_id, feature.seqid, feature.strand, feature.source
            )
            order.append(transcript_id)
    for feature in features:
        if feature.feature_type not in {"exon", "CDS"}:
            continue
        transcript_id = feature.attributes.get("transcript_id", "")
        if transcript_id not in models:
            raise UtrError(f"coding child has unknown transcript: {transcript_id!r}")
        target = models[transcript_id].exons if feature.feature_type == "exon" else models[transcript_id].cds
        target.append(feature)
    for transcript_id in order:
        model = models[transcript_id]
        if not model.exons or not model.cds:
            raise UtrError(f"coding transcript lacks exon/CDS: {transcript_id}")
        if any(feature.phase not in {"0", "1", "2"} for feature in model.cds):
            raise UtrError(f"coding CDS has invalid phase: {transcript_id}")
    return order, models


def load_raw_gff3(path: Path) -> tuple[list[str], dict[str, TranscriptModel]]:
    features = _features(path, "gff3")
    order: list[str] = []
    models: dict[str, TranscriptModel] = {}
    for feature in features:
        if feature.feature_type not in {"transcript", "mRNA"}:
            continue
        attrs = feature.attributes
        transcript_id = attrs.get("ID") or attrs.get("transcript_id", "")
        gene_id = attrs.get("Parent") or attrs.get("gene_id", "")
        if not transcript_id or not gene_id or transcript_id in models:
            raise UtrError(f"invalid/duplicate refined transcript: {transcript_id!r}")
        models[transcript_id] = TranscriptModel(
            gene_id=gene_id,
            transcript_id=transcript_id,
            candidate_id=attrs.get("candidate_id") or attrs.get("candidate_key", ""),
            seqid=feature.seqid,
            strand=feature.strand,
            source=feature.source,
            utr_refine_evidence=attrs.get("utr_refine_evidence", ""),
        )
        order.append(transcript_id)
    for feature in features:
        if feature.feature_type not in {"exon", "CDS", "five_prime_UTR", "three_prime_UTR"}:
            continue
        parents = [value for value in feature.attributes.get("Parent", "").split(",") if value]
        if not parents:
            raise UtrError(f"refined child has no Parent at {feature.seqid}:{feature.start}")
        known = [parent for parent in parents if parent in models]
        if not known:
            raise UtrError(f"refined child has no known Parent at {feature.seqid}:{feature.start}")
        for parent in known:
            model = models[parent]
            target = {
                "exon": model.exons,
                "CDS": model.cds,
                "five_prime_UTR": model.five_prime_utrs,
                "three_prime_UTR": model.three_prime_utrs,
            }[feature.feature_type]
            target.append(feature)
    return order, models


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def _subtract(outer: Sequence[tuple[int, int]], inner: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    remaining: list[tuple[int, int]] = []
    for outer_start, outer_end in outer:
        cursor = outer_start
        for inner_start, inner_end in inner:
            if inner_end < cursor:
                continue
            if inner_start > outer_end:
                break
            if inner_start > cursor:
                remaining.append((cursor, min(outer_end, inner_start - 1)))
            cursor = max(cursor, inner_end + 1)
            if cursor > outer_end:
                break
        if cursor <= outer_end:
            remaining.append((cursor, outer_end))
    return remaining


def _new_feature(model: TranscriptModel, feature_type: str, start: int, end: int) -> Feature:
    return Feature(model.seqid, "InGenAnnot", feature_type, start, end, model.strand, ".", {})


def normalize_refinement(
    coding_gtf: Path,
    raw_gff3: Path | None,
) -> tuple[list[str], dict[str, TranscriptModel]]:
    expected_order, expected = load_coding_gtf(coding_gtf)
    if raw_gff3 is None:
        return expected_order, expected
    observed_order, observed = load_raw_gff3(raw_gff3)
    if observed_order != expected_order:
        raise UtrError("InGenAnnot changed selected transcript order or membership")
    for transcript_id in expected_order:
        source = expected[transcript_id]
        result = observed[transcript_id]
        if (result.gene_id, result.seqid, result.strand) != (source.gene_id, source.seqid, source.strand):
            raise UtrError(f"InGenAnnot changed hierarchy/axis: {transcript_id}")
        if result.candidate_id in EMPTY:
            result.candidate_id = source.candidate_id
        if result.candidate_id != source.candidate_id:
            raise UtrError(f"InGenAnnot changed candidate identity: {transcript_id}")
        expected_cds = sorted(feature.signature for feature in source.cds)
        observed_cds = sorted(feature.signature for feature in result.cds)
        if expected_cds != observed_cds:
            raise UtrError(f"InGenAnnot changed CDS coordinates or phase: {transcript_id}")
        if not result.exons:
            raise UtrError(f"refined transcript lacks exon: {transcript_id}")
        if any(feature.seqid != result.seqid or feature.strand != result.strand for feature in result.exons + result.cds):
            raise UtrError(f"refined child changes sequence/strand: {transcript_id}")
        exon_intervals = _merge_intervals((feature.start, feature.end) for feature in result.exons)
        # InGenAnnot may repeat a shared, multi-parent exon under each isoform.
        # Keep one genomic block per transcript for both GFF and sequence export.
        originals = {(feature.start, feature.end): feature for feature in reversed(result.exons)}
        result.exons = [
            originals.get(interval) or _new_feature(result, "exon", *interval)
            for interval in exon_intervals
        ]
        cds_intervals = _merge_intervals((feature.start, feature.end) for feature in result.cds)
        residual = _subtract(exon_intervals, cds_intervals)
        min_cds, max_cds = cds_intervals[0][0], cds_intervals[-1][1]
        lower: list[tuple[int, int]] = []
        upper: list[tuple[int, int]] = []
        for interval in residual:
            if interval[1] < min_cds:
                lower.append(interval)
            elif interval[0] > max_cds:
                upper.append(interval)
            else:
                raise UtrError(f"non-CDS exon lies inside CDS span: {transcript_id}")
        if (lower or upper) and result.utr_refine_evidence in EMPTY:
            raise UtrError(f"UTR proposal lacks RNA evidence ID: {transcript_id}")
        five, three = (lower, upper) if result.strand == "+" else (upper, lower)
        result.five_prime_utrs = [_new_feature(result, "five_prime_UTR", *interval) for interval in five]
        result.three_prime_utrs = [_new_feature(result, "three_prime_UTR", *interval) for interval in three]
    return expected_order, observed


def restore_postqc_gene_hierarchy(
    order: Sequence[str],
    models: Mapping[str, TranscriptModel],
    preparation_rows: Mapping[str, Mapping[str, str]],
) -> None:
    """Replace per-transcript refinement containers with final post-QC genes."""

    for transcript_id in order:
        model = models[transcript_id]
        if model.candidate_id not in preparation_rows:
            raise UtrError(f"missing preparation mapping: {model.candidate_id}")
        preparation = preparation_rows[model.candidate_id]
        postqc_gene_id = preparation.get("postqc_gene_id", "").strip()
        if postqc_gene_id in EMPTY:
            raise UtrError(
                f"missing final post-QC gene assignment: {model.candidate_id}"
            )
        if model.gene_id != preparation["temporary_gene_id"]:
            raise UtrError(f"UTR refinement gene lineage mismatch: {transcript_id}")
        model.gene_id = postqc_gene_id


def _write_gff_row(
    handle: TextIO,
    model: TranscriptModel,
    feature_type: str,
    start: int,
    end: int,
    phase: str,
    attributes: Iterable[tuple[str, str]],
    source: str,
) -> None:
    handle.write(
        f"{model.seqid}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{model.strand}\t{phase}\t"
        f"{format_gff3_attributes(attributes)}\n"
    )


def write_normalized_shard(
    path: Path,
    mapping_path: Path,
    order: Sequence[str],
    models: Mapping[str, TranscriptModel],
    preparation_rows: Mapping[str, Mapping[str, str]],
) -> list[dict[str, object]]:
    mapping_rows: list[dict[str, object]] = []
    genes: defaultdict[str, list[str]] = defaultdict(list)
    gene_order: list[str] = []
    for transcript_id in order:
        model = models[transcript_id]
        if model.gene_id not in genes:
            gene_order.append(model.gene_id)
        genes[model.gene_id].append(transcript_id)
    with atomic_text(path) as handle:
        handle.write("##gff-version 3\n")
        for gene_id in gene_order:
            transcript_ids = genes[gene_id]
            gene_models = [models[value] for value in transcript_ids]
            axes = {(model.seqid, model.strand) for model in gene_models}
            if len(axes) != 1:
                raise UtrError(f"gene crosses sequence/strand axes: {gene_id}")
            representative = gene_models[0]
            _write_gff_row(
                handle, representative, "gene",
                min(model.start for model in gene_models),
                max(model.end for model in gene_models), ".",
                (("ID", gene_id), ("transcript_count", str(len(gene_models)))),
                "GeneInIt",
            )
            for transcript_id in transcript_ids:
                model = models[transcript_id]
                if model.candidate_id not in preparation_rows:
                    raise UtrError(f"missing preparation mapping: {model.candidate_id}")
                tx_attrs = [("ID", model.transcript_id), ("Parent", model.gene_id), ("candidate_id", model.candidate_id)]
                if model.utr_refine_evidence not in EMPTY:
                    tx_attrs.append(("utr_refine_evidence", model.utr_refine_evidence))
                _write_gff_row(handle, model, "mRNA", model.start, model.end, ".", tx_attrs, "GeneInIt")
                feature_groups = (
                    ("exon", model.exons), ("CDS", model.cds),
                    ("five_prime_UTR", model.five_prime_utrs),
                    ("three_prime_UTR", model.three_prime_utrs),
                )
                for feature_type, features in feature_groups:
                    for index, feature in enumerate(sorted(features, key=lambda item: (item.start, item.end, item.phase)), 1):
                        phase = feature.phase if feature_type == "CDS" else "."
                        _write_gff_row(
                            handle, model, feature_type, feature.start, feature.end, phase,
                            (("ID", f"{feature_type}-{model.transcript_id}-{index:03d}"), ("Parent", model.transcript_id)),
                            feature.source or "GeneInIt",
                        )
                source_row = preparation_rows[model.candidate_id]
                mapping_rows.append(
                    {
                        **{column: source_row[column] for column in PREPARATION_COLUMNS},
                        "utr_refine_evidence": model.utr_refine_evidence or ".",
                        "five_prime_utr_segments": len(model.five_prime_utrs),
                        "five_prime_utr_bp": sum(item.end - item.start + 1 for item in model.five_prime_utrs),
                        "three_prime_utr_segments": len(model.three_prime_utrs),
                        "three_prime_utr_bp": sum(item.end - item.start + 1 for item in model.three_prime_utrs),
                    }
                )
    _write_table(mapping_path, UTR_MAPPING_COLUMNS, mapping_rows)
    return mapping_rows


def prepare_selected_coding_models(
    candidate_table: Path,
    selection_audit: Path,
    reference_sequence_ids: Sequence[str],
    output_gff3: Path,
    output_mapping: Path,
) -> list[dict[str, object]]:
    """Materialize the selected, coding-only transcript set before post-QC."""

    by_id, selected = _retained_candidate_metadata(
        candidate_table, selection_audit, None, require_structure=True
    )
    if not selected:
        raise UtrError("selection audit contains no selected models")
    missing = set(selected) - set(by_id)
    if missing:
        raise UtrError(f"selection audit references unknown candidates: {len(missing)}")
    sequence_rank = {seqid: index for index, seqid in enumerate(reference_sequence_ids)}
    if len(sequence_rank) != len(reference_sequence_ids):
        raise UtrError("reference sequence IDs contain duplicates")
    unknown = {by_id[candidate_id]["seqid"] for candidate_id in selected} - set(sequence_rank)
    if unknown:
        raise UtrError(f"selected models use unknown reference seqids: {sorted(unknown)}")
    ordered_rows = sorted(
        (by_id[candidate_id] for candidate_id in selected),
        key=lambda row: (
            sequence_rank[row["seqid"]], int(row["start"]), int(row["end"]),
            row["strand"], row["candidate_id"],
        ),
    )
    order: list[str] = []
    models: dict[str, TranscriptModel] = {}
    preparation: dict[str, dict[str, object]] = {}
    for index, row in enumerate(ordered_rows, 1):
        candidate_id = row["candidate_id"]
        gene_id = f"GITMPG{index:07d}"
        transcript_id = f"{gene_id}.t1"
        parts = parse_cds_signature(row["cds_signature"])
        model = TranscriptModel(
            gene_id=gene_id,
            transcript_id=transcript_id,
            candidate_id=candidate_id,
            seqid=row["seqid"],
            strand=row["strand"],
            source="GeneInIt",
            exons=[
                Feature(row["seqid"], "GeneInIt", "exon", part.start, part.end, row["strand"], ".", {})
                for part in parts
            ],
            cds=[
                Feature(
                    row["seqid"], "GeneInIt", "CDS", part.start, part.end,
                    row["strand"], str(part.phase), {},
                )
                for part in parts
            ],
        )
        order.append(transcript_id)
        models[transcript_id] = model
        audit = selected[candidate_id]
        preparation[candidate_id] = {
            "candidate_id": candidate_id,
            "temporary_gene_id": gene_id,
            "temporary_transcript_id": transcript_id,
            "seqid": row["seqid"],
            "strand": row["strand"],
            "start": parts[0].start,
            "end": parts[-1].end,
            "selection_locus_id": audit["locus_id"],
            "predictors": row["predictors"],
            "source_gene_ids": row["source_gene_ids"],
            "source_transcript_ids": row["source_transcript_ids"],
            "routes": audit["routes"],
            "active_fusion": audit["active_fusion"],
        }
    return write_normalized_shard(
        output_gff3, output_mapping, order, models, preparation
    )


def _write_rna_gtf_shards(
    merged_rna_gtf: Path,
    seqid_to_shard: Mapping[str, str],
    rna_paths: Mapping[str, Path],
) -> Counter[str]:
    if set(rna_paths) != set(seqid_to_shard.values()):
        raise UtrError("RNA output shards differ from the seqid/shard contract")
    handles: dict[str, TextIO] = {}
    temporary_paths: dict[str, Path] = {}
    transcript_counts: Counter[str] = Counter()
    try:
        for shard, output in rna_paths.items():
            output.parent.mkdir(parents=True, exist_ok=True)
            temporary = output.with_name(f".{output.name}.partial.{os.getpid()}")
            temporary_paths[shard] = temporary
            handles[shard] = temporary.open("w", encoding="utf-8")
        with merged_rna_gtf.open(encoding="utf-8") as source:
            for number, raw in enumerate(source, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise UtrError(f"{merged_rna_gtf}:{number} does not have 9 columns")
                seqid = fields[0]
                if seqid not in seqid_to_shard:
                    raise UtrError(f"RNA transcript uses unknown reference seqid: {seqid}")
                handles[seqid_to_shard[seqid]].write(
                    raw if raw.endswith("\n") else raw + "\n"
                )
                if fields[2] in {"transcript", "mRNA"}:
                    transcript_counts[seqid] += 1
    finally:
        for handle in handles.values():
            handle.close()
    for shard, output in rna_paths.items():
        os.replace(temporary_paths[shard], output)
    return transcript_counts


def prepare_postqc_utr_shards(
    postqc_gff3: Path,
    transcript_decisions: Path,
    candidate_table: Path,
    selection_audit: Path,
    merged_rna_gtf: Path,
    seqid_to_shard: Mapping[str, str],
    coding_paths: Mapping[str, Path],
    rna_paths: Mapping[str, Path],
    preparation_mapping: Path,
    task_manifest: Path,
) -> list[dict[str, object]]:
    """Prepare every active post-QC transcript for the single final UTR pass."""

    order, models = load_raw_gff3(postqc_gff3)
    if not order:
        raise UtrError("post-QC GFF3 contains no active transcripts")
    decision_rows = _read_table(
        transcript_decisions,
        {
            "transcript_id", "candidate_id", "postqc_gene_id",
            "selection_locus_id", "routes", "source_status",
        },
    )
    decision_ids = [row["transcript_id"] for row in decision_rows]
    if len(decision_ids) != len(set(decision_ids)):
        raise UtrError("post-QC transcript decisions contain duplicate transcript IDs")
    decisions = {
        row["transcript_id"]: row
        for row in decision_rows
        if row["postqc_gene_id"] not in EMPTY
    }
    if set(decisions) != set(models):
        raise UtrError("post-QC GFF3 and active transcript decisions differ")
    released_candidate_ids = {
        row["candidate_id"]
        for row in decisions.values()
        if row["source_status"] == "RELEASED"
    }
    candidate_metadata, audit_metadata = _retained_candidate_metadata(
        candidate_table, selection_audit, released_candidate_ids,
        require_structure=False,
    )
    if any(model.candidate_id in EMPTY for model in models.values()):
        raise UtrError("post-QC transcript lacks candidate identity")
    if len({model.candidate_id for model in models.values()}) != len(models):
        raise UtrError("post-QC transcripts do not have one-to-one candidate identities")
    if set(coding_paths) != set(seqid_to_shard.values()):
        raise UtrError("coding output shards differ from the seqid/shard contract")

    gene_members: defaultdict[str, list[str]] = defaultdict(list)
    gene_order: list[str] = []
    for transcript_id in order:
        model = models[transcript_id]
        decision = decisions[transcript_id]
        if model.candidate_id != decision["candidate_id"]:
            raise UtrError(f"post-QC candidate lineage mismatch: {transcript_id}")
        if model.gene_id != decision["postqc_gene_id"]:
            raise UtrError(f"post-QC gene assignment mismatch: {transcript_id}")
        if model.five_prime_utrs or model.three_prime_utrs:
            raise UtrError(
                "post-QC input already contains UTR features; final UTR reconstruction "
                f"must start from coding-only models: {transcript_id}"
            )
        if _merge_intervals((feature.start, feature.end) for feature in model.exons) != _merge_intervals(
            (feature.start, feature.end) for feature in model.cds
        ):
            raise UtrError(
                "post-QC input exons extend beyond CDS before final UTR reconstruction: "
                f"{transcript_id}"
            )
        if model.gene_id not in gene_members:
            gene_order.append(model.gene_id)
        gene_members[model.gene_id].append(transcript_id)
    for gene_id in gene_order:
        axes = {
            (models[transcript_id].seqid, models[transcript_id].strand)
            for transcript_id in gene_members[gene_id]
        }
        if len(axes) != 1:
            raise UtrError(f"post-QC gene crosses sequence/strand axes: {gene_id}")

    refinement_gene_ids = {
        transcript_id: f"GIUTRG{index:07d}"
        for index, transcript_id in enumerate(order, 1)
    }
    temporary_paths: dict[str, Path] = {}
    handles: dict[str, TextIO] = {}
    try:
        for shard, output in coding_paths.items():
            output.parent.mkdir(parents=True, exist_ok=True)
            temporary = output.with_name(f".{output.name}.partial.{os.getpid()}")
            temporary_paths[shard] = temporary
            handles[shard] = temporary.open("w", encoding="utf-8")
        for transcript_id in order:
            model = models[transcript_id]
            seqid, strand = model.seqid, model.strand
            if seqid not in seqid_to_shard:
                raise UtrError(f"post-QC model uses unknown reference seqid: {seqid}")
            handle = handles[seqid_to_shard[seqid]]
            refinement_gene_id = refinement_gene_ids[transcript_id]
            gene_attrs = _gtf_attributes((("gene_id", refinement_gene_id),))
            handle.write(
                f"{seqid}\tGeneInIt-postQC\tgene\t"
                f"{model.start}\t{model.end}\t.\t{strand}\t.\t{gene_attrs}\n"
            )
            tx_attrs = _gtf_attributes(
                (("gene_id", refinement_gene_id), ("transcript_id", transcript_id),
                 ("candidate_id", model.candidate_id))
            )
            handle.write(
                f"{seqid}\tGeneInIt-postQC\ttranscript\t{model.start}\t{model.end}\t.\t"
                f"{strand}\t.\t{tx_attrs}\n"
            )
            for feature_type, features in (("exon", model.exons), ("CDS", model.cds)):
                for feature in sorted(
                    features, key=lambda item: (item.start, item.end, item.phase)
                ):
                    phase = feature.phase if feature_type == "CDS" else "."
                    handle.write(
                        f"{seqid}\t{feature.source or 'GeneInIt-postQC'}\t{feature_type}\t"
                        f"{feature.start}\t{feature.end}\t.\t{strand}\t{phase}\t{tx_attrs}\n"
                    )
    finally:
        for handle in handles.values():
            handle.close()
    for shard, output in coding_paths.items():
        os.replace(temporary_paths[shard], output)

    preparation_rows: list[dict[str, object]] = []
    coding_counts: Counter[str] = Counter()
    for transcript_id in order:
        model = models[transcript_id]
        decision = decisions[transcript_id]
        candidate = candidate_metadata.get(model.candidate_id, {})
        audit = audit_metadata.get(model.candidate_id, {})
        preparation_rows.append(
            {
                "candidate_id": model.candidate_id,
                "temporary_gene_id": refinement_gene_ids[transcript_id],
                "temporary_transcript_id": model.transcript_id,
                "seqid": model.seqid,
                "strand": model.strand,
                "start": model.start,
                "end": model.end,
                "selection_locus_id": decision["selection_locus_id"] or ".",
                "predictors": candidate.get("predictors", decision["source_status"]),
                "source_gene_ids": candidate.get("source_gene_ids", "."),
                "source_transcript_ids": candidate.get("source_transcript_ids", "."),
                "routes": decision["routes"] or ".",
                "active_fusion": audit.get("active_fusion", "0"),
                "postqc_gene_id": model.gene_id,
                "source_status": decision["source_status"],
            }
        )
        coding_counts[model.seqid] += 1
    _write_table(
        preparation_mapping, POSTQC_PREPARATION_COLUMNS, preparation_rows
    )
    rna_counts = _write_rna_gtf_shards(merged_rna_gtf, seqid_to_shard, rna_paths)
    task_rows = [
        {
            "index": index,
            "seqid": seqid,
            "shard": shard,
            "coding_models": coding_counts[seqid],
            "rna_transcripts": rna_counts[seqid],
        }
        for index, (seqid, shard) in enumerate(seqid_to_shard.items(), 1)
    ]
    _write_table(
        task_manifest,
        ("index", "seqid", "shard", "coding_models", "rna_transcripts"),
        task_rows,
    )
    return preparation_rows


def aggregate_utr_shards(
    normalized_paths: Sequence[Path],
    mapping_paths: Sequence[Path],
    preparation_mapping: Path,
    output_gff3: Path,
    output_mapping: Path,
    summary_path: Path,
    complete_path: Path,
) -> None:
    if len(normalized_paths) != len(mapping_paths):
        raise UtrError("normalized GFF3 and mapping shard counts differ")
    expected_rows = _read_table(preparation_mapping, set(PREPARATION_COLUMNS))
    expected = {row["candidate_id"] for row in expected_rows}
    observed_rows: list[dict[str, str]] = []
    with atomic_text(output_gff3) as handle:
        handle.write("##gff-version 3\n")
        for path in normalized_paths:
            with path.open(encoding="utf-8") as source:
                for raw in source:
                    if raw.strip() and not raw.startswith("#"):
                        handle.write(raw if raw.endswith("\n") else raw + "\n")
    for path in mapping_paths:
        observed_rows.extend(_read_table(path, set(UTR_MAPPING_COLUMNS)))
    observed = [row["candidate_id"] for row in observed_rows]
    if len(observed) != len(set(observed)):
        raise UtrError("UTR shards contain duplicate candidate IDs")
    if set(observed) != expected:
        raise UtrError(
            f"UTR candidate membership differs: missing={len(expected-set(observed))}, extra={len(set(observed)-expected)}"
        )
    by_id = {row["candidate_id"]: row for row in observed_rows}
    ordered_rows = [by_id[row["candidate_id"]] for row in expected_rows]
    _write_table(output_mapping, UTR_MAPPING_COLUMNS, ordered_rows)
    updated = [row for row in ordered_rows if row["utr_refine_evidence"] not in EMPTY]
    metrics = (
        ("status", "PASS_MIGRATION_PENDING"),
        ("selected_model_count", len(ordered_rows)),
        ("utr_updated_model_count", len(updated)),
        ("five_prime_utr_bp", sum(int(row["five_prime_utr_bp"]) for row in ordered_rows)),
        ("three_prime_utr_bp", sum(int(row["three_prime_utr_bp"]) for row in ordered_rows)),
    )
    _write_table(summary_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in metrics))
    completion = (
        ("schema_version", 1), ("status", "PASS_MIGRATION_PENDING"),
        ("selected_model_count", len(ordered_rows)), ("utr_updated_model_count", len(updated)),
        ("preparation_mapping_sha256", sha256_file(preparation_mapping)),
        ("utr_gff3_sha256", sha256_file(output_gff3)),
        ("utr_mapping_sha256", sha256_file(output_mapping)),
        ("summary_sha256", sha256_file(summary_path)),
    )
    _write_table(complete_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in completion))
