"""Build deterministic structural release files from the final UTR GFF3."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO

from geneinit.evidence.te_classification import TeClassifier, confirmed_only, model_identity, protein_sha256
from geneinit.models.candidate import CdsPart, audit_orf, read_fasta, reverse_complement
from geneinit.utr.core import (
    EMPTY,
    Feature,
    TranscriptModel,
    _write_table,
    atomic_text,
    format_gff3_attributes,
    load_raw_gff3,
    sha256_file,
)


GENE_MAP_COLUMNS = (
    "postqc_gene_id", "final_gene_id", "seqid", "strand", "start", "end",
    "transcript_count",
)
TRANSCRIPT_MAP_COLUMNS = (
    "temporary_transcript_id", "candidate_id", "postqc_gene_id", "final_gene_id",
    "final_transcript_id", "canonical_transcript", "transcript_rank",
)
LINEAGE_COLUMNS = (
    "candidate_id", "temporary_transcript_id", "final_gene_id",
    "final_transcript_id", "source_status", "selection_locus_id", "routes",
    "postqc_operation", "postqc_reason",
)
REVIEW_COLUMNS = (
    "flag_type", "severity", "subject_a", "subject_b", "seqid", "start", "end",
    "reason", "release_blocking",
)


class ReleaseError(RuntimeError):
    """Raised when a canonical structural release fails validation."""


def _read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ReleaseError(f"{path} lacks columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def _feature_intervals(features: Iterable[Feature]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((feature.start, feature.end) for feature in features))


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return tuple(merged)


def _splice_sequence(
    features: Iterable[Feature], genome: Mapping[str, str], seqid: str, strand: str
) -> str:
    if seqid not in genome:
        raise ReleaseError(f"annotation seqid absent from genome: {seqid}")
    sequence = genome[seqid]
    ordered = sorted(features, key=lambda item: item.start, reverse=strand == "-")
    parts: list[str] = []
    for feature in ordered:
        if feature.start < 1 or feature.end > len(sequence):
            raise ReleaseError(f"feature exceeds genome: {seqid}:{feature.start}-{feature.end}")
        part = sequence[feature.start - 1 : feature.end]
        parts.append(reverse_complement(part) if strand == "-" else part)
    return "".join(parts).upper()


def _write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> None:
    with atomic_text(path) as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for offset in range(0, len(sequence), width):
                handle.write(sequence[offset : offset + width] + "\n")


def _write_gff_row(
    handle: TextIO, seqid: str, source: str, feature_type: str, start: int, end: int,
    strand: str, phase: str, attributes: Iterable[tuple[str, str]],
) -> None:
    handle.write(
        f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{strand}\t{phase}\t"
        f"{format_gff3_attributes(attributes)}\n"
    )


def _structural_validate(model: TranscriptModel, genome: Mapping[str, str]) -> tuple[str, str, str]:
    if not model.exons or not model.cds:
        raise ReleaseError(f"transcript lacks exon/CDS: {model.transcript_id}")
    if any(feature.strand != model.strand or feature.seqid != model.seqid for feature in model.exons + model.cds):
        raise ReleaseError(f"mixed child axis: {model.transcript_id}")
    ordered_exons = sorted(model.exons, key=lambda feature: (feature.start, feature.end))
    if any(right.start <= left.end + 1 for left, right in zip(ordered_exons, ordered_exons[1:])):
        raise ReleaseError(f"duplicate, overlapping or adjacent final exons: {model.transcript_id}")
    if _merge_intervals(_feature_intervals(model.exons)) != _merge_intervals(
        _feature_intervals(model.cds + model.five_prime_utrs + model.three_prime_utrs)
    ):
        raise ReleaseError(f"exons differ from CDS+UTR: {model.transcript_id}")
    parts = tuple(
        sorted(
            CdsPart(feature.start, feature.end, int(feature.phase))
            for feature in model.cds
        )
    )
    audit = audit_orf(model.seqid, model.strand, parts, genome)
    if audit.status != "complete":
        raise ReleaseError(f"non-complete final ORF {model.transcript_id}: {audit.status}")
    transcript_sequence = _splice_sequence(model.exons, genome, model.seqid, model.strand)
    cds_sequence = _splice_sequence(model.cds, genome, model.seqid, model.strand)
    if cds_sequence != audit.coding_sequence:
        raise ReleaseError(f"CDS extraction differs from ORF audit: {model.transcript_id}")
    protein = audit.protein_sequence[:-1] if audit.protein_sequence.endswith("*") else audit.protein_sequence
    if not protein or "*" in protein:
        raise ReleaseError(f"final protein is empty or has internal stop: {model.transcript_id}")
    return transcript_sequence, cds_sequence, protein


def build_structural_release(
    final_utr_gff3: Path,
    transcript_decisions: Path,
    edge_audit: Path,
    bridge_reviews: Path,
    genome_path: Path,
    policy_path: Path,
    reference_sequence_ids: Sequence[str],
    assembly_id: str,
    gene_prefix: str,
    annotation_gff3: Path,
    transcripts_fna: Path,
    cds_fna: Path,
    proteins_faa: Path,
    longest_proteins_faa: Path,
    gene_id_map: Path,
    transcript_id_map: Path,
    candidate_lineage: Path,
    structural_decisions: Path,
    review_flags: Path,
    qc_tsv: Path,
    qc_json: Path,
    release_status: Path,
    manifest_path: Path,
    checksums_path: Path,
    frozen_policy_path: Path,
    complete_path: Path,
    te_evidence_path: Path | None = None,
    candidates_path: Path | None = None,
    selection_audit_path: Path | None = None,
    generated_te_path: Path | None = None,
) -> None:
    genome = read_fasta(genome_path)
    _order, models_by_transcript = load_raw_gff3(final_utr_gff3)
    decision_rows = _read_table(
        transcript_decisions,
        {
            "transcript_id", "candidate_id", "source_status", "input_gene_id",
            "postqc_gene_id", "canonical_transcript", "canonical_rank",
            "selection_locus_id", "routes", "operation", "reason",
        },
    )
    if not decision_rows:
        raise ReleaseError("post-QC transcript decisions are empty")
    decision_ids = [row["transcript_id"] for row in decision_rows]
    if len(decision_ids) != len(set(decision_ids)):
        raise ReleaseError("post-QC transcript decisions contain duplicate transcript IDs")
    decisions = {
        row["transcript_id"]: row for row in decision_rows if row["postqc_gene_id"] not in EMPTY
    }
    if set(decisions) != set(models_by_transcript):
        raise ReleaseError("post-QC GFF3 and active transcript decisions differ")
    with policy_path.open(encoding="utf-8") as handle:
        policy = json.load(handle)
    if not policy.get("policy_id"):
        raise ReleaseError("policy lacks policy_id")
    policy_status = str(policy.get("status", "MIGRATION_PENDING"))
    classifier = TeClassifier(te_evidence_path, assembly_id) if confirmed_only(policy.get("selector_core", {})) else None
    te_annotations = {}
    selection_te_rows = {}
    if classifier is not None:
        if candidates_path is None or selection_audit_path is None or generated_te_path is None:
            raise ReleaseError("confirmed-only release requires raw and post-QC classification audits")
        selection_te_rows = {r["candidate_id"]: r for r in _read_table(selection_audit_path, {"candidate_id", "te_policy_id"})}
        for model in models_by_transcript.values():
            identity = {**selection_te_rows.get(model.candidate_id, {}), **model_identity(model)}
            result = classifier.classify(identity, protein_sha256(_structural_validate(model, genome)[2]))
            if result["candidate_pool"] != "GENERAL":
                raise ReleaseError(f"confirmed TE in general release: {model.candidate_id}")
            te_annotations[model.candidate_id] = result

    genes: defaultdict[str, list[TranscriptModel]] = defaultdict(list)
    for model in models_by_transcript.values():
        genes[model.gene_id].append(model)
        decision = decisions[model.transcript_id]
        if model.candidate_id != decision["candidate_id"]:
            raise ReleaseError(f"candidate lineage mismatch: {model.transcript_id}")
        if model.gene_id != decision["postqc_gene_id"]:
            raise ReleaseError(
                "final UTR gene hierarchy differs from post-QC decision: "
                f"{model.transcript_id} ({model.gene_id} != "
                f"{decision['postqc_gene_id']})"
            )
    if len(reference_sequence_ids) != len(set(reference_sequence_ids)):
        raise ReleaseError("reference_sequence_ids contains duplicates")
    sequence_rank = {seqid: index for index, seqid in enumerate(reference_sequence_ids)}
    missing_seqids = {model.seqid for model in models_by_transcript.values()} - set(sequence_rank)
    if missing_seqids:
        raise ReleaseError(f"annotation uses unregistered reference seqids: {sorted(missing_seqids)}")
    for gene_id, gene_models in genes.items():
        axes = {(model.seqid, model.strand) for model in gene_models}
        if len(axes) != 1:
            raise ReleaseError(f"post-QC gene crosses sequence/strand axes: {gene_id}")
    ordered_genes = sorted(
        genes,
        key=lambda gene_id: (
            sequence_rank[genes[gene_id][0].seqid],
            min(model.start for model in genes[gene_id]),
            max(model.end for model in genes[gene_id]),
            genes[gene_id][0].strand,
            gene_id,
        ),
    )
    gene_map: dict[str, str] = {
        gene_id: f"{gene_prefix}{index:06d}" for index, gene_id in enumerate(ordered_genes, 1)
    }
    transcript_map: dict[str, str] = {}
    ranked_by_gene: dict[str, list[TranscriptModel]] = {}
    for gene_id in ordered_genes:
        ranked = sorted(
            genes[gene_id],
            key=lambda model: (int(decisions[model.transcript_id]["canonical_rank"]), model.candidate_id),
        )
        ranks = [int(decisions[model.transcript_id]["canonical_rank"]) for model in ranked]
        if ranks != list(range(1, len(ranked) + 1)):
            raise ReleaseError(f"post-QC canonical ranks are not contiguous: {gene_id}")
        canonical_flags = [int(decisions[model.transcript_id]["canonical_transcript"]) for model in ranked]
        if canonical_flags != [1] + [0] * (len(ranked) - 1):
            raise ReleaseError(f"post-QC canonical flags disagree with ranks: {gene_id}")
        ranked_by_gene[gene_id] = ranked
        final_gene = gene_map[gene_id]
        for rank, model in enumerate(ranked, 1):
            transcript_map[model.transcript_id] = f"{final_gene}.t{rank}"

    sequence_records: dict[str, tuple[str, str, str]] = {}
    with atomic_text(annotation_gff3) as handle:
        handle.write("##gff-version 3\n")
        for postqc_gene in ordered_genes:
            final_gene = gene_map[postqc_gene]
            ranked = ranked_by_gene[postqc_gene]
            gene_start = min(model.start for model in ranked)
            gene_end = max(model.end for model in ranked)
            _write_gff_row(
                handle, ranked[0].seqid, "GeneInIt", "gene", gene_start, gene_end,
                ranked[0].strand, ".",
                (("ID", final_gene), ("transcript_count", str(len(ranked))),
                 ("canonical_transcript", transcript_map[ranked[0].transcript_id]),
                 ("postqc_gene_id", postqc_gene)),
            )
            for rank, model in enumerate(ranked, 1):
                final_transcript = transcript_map[model.transcript_id]
                decision = decisions[model.transcript_id]
                sequence_records[final_transcript] = _structural_validate(model, genome)
                attrs = [
                    ("ID", final_transcript), ("Parent", final_gene),
                    ("candidate_id", model.candidate_id),
                    ("canonical", "true" if rank == 1 else "false"),
                    ("source_status", decision["source_status"]),
                ]
                if classifier is not None:
                    te = te_annotations[model.candidate_id]
                    attrs.extend([("candidate_pool", "GENERAL"), ("te_class", te["effective_class"]),
                                  ("te_review_flags", te["te_review_flags"]), ("host_state", te["host_state"])])
                if model.utr_refine_evidence not in EMPTY:
                    attrs.append(("utr_refine_evidence", model.utr_refine_evidence))
                _write_gff_row(
                    handle, model.seqid, model.source, "mRNA", model.start, model.end,
                    model.strand, ".", attrs,
                )
                for feature_type, features in (
                    ("exon", model.exons), ("CDS", model.cds),
                    ("five_prime_UTR", model.five_prime_utrs),
                    ("three_prime_UTR", model.three_prime_utrs),
                ):
                    for index, feature in enumerate(sorted(features, key=lambda item: (item.start, item.end, item.phase)), 1):
                        _write_gff_row(
                            handle, model.seqid, feature.source, feature_type,
                            feature.start, feature.end, model.strand,
                            feature.phase if feature_type == "CDS" else ".",
                            (("ID", f"{feature_type}-{final_transcript}-{index:03d}"),
                             ("Parent", final_transcript)),
                        )

    ordered_transcripts = [
        transcript_map[model.transcript_id]
        for gene_id in ordered_genes for model in ranked_by_gene[gene_id]
    ]
    _write_fasta(transcripts_fna, ((identifier, sequence_records[identifier][0]) for identifier in ordered_transcripts))
    _write_fasta(cds_fna, ((identifier, sequence_records[identifier][1]) for identifier in ordered_transcripts))
    _write_fasta(proteins_faa, ((identifier, sequence_records[identifier][2]) for identifier in ordered_transcripts))
    longest_records: list[tuple[str, str]] = []
    for postqc_gene in ordered_genes:
        identifiers = [transcript_map[model.transcript_id] for model in ranked_by_gene[postqc_gene]]
        longest = max(identifiers, key=lambda value: (len(sequence_records[value][2]), -int(value.rsplit("t", 1)[1])))
        longest_records.append((longest, sequence_records[longest][2]))
    _write_fasta(longest_proteins_faa, longest_records)

    gene_rows = [
        {
            "postqc_gene_id": gene_id, "final_gene_id": gene_map[gene_id],
            "seqid": genes[gene_id][0].seqid, "strand": genes[gene_id][0].strand,
            "start": min(model.start for model in genes[gene_id]),
            "end": max(model.end for model in genes[gene_id]),
            "transcript_count": len(genes[gene_id]),
        }
        for gene_id in ordered_genes
    ]
    _write_table(gene_id_map, GENE_MAP_COLUMNS, gene_rows)
    transcript_rows = []
    lineage_rows = []
    structural_rows = []
    for row in decision_rows:
        active = row["postqc_gene_id"] not in EMPTY
        final_gene = gene_map[row["postqc_gene_id"]] if active else "."
        final_transcript = transcript_map[row["transcript_id"]] if active else "."
        if active:
            transcript_rows.append(
                {
                    "temporary_transcript_id": row["transcript_id"],
                    "candidate_id": row["candidate_id"],
                    "postqc_gene_id": row["postqc_gene_id"],
                    "final_gene_id": final_gene,
                    "final_transcript_id": final_transcript,
                    "canonical_transcript": row["canonical_transcript"],
                    "transcript_rank": row["canonical_rank"],
                }
            )
        lineage_rows.append(
            {
                "candidate_id": row["candidate_id"],
                "temporary_transcript_id": row["transcript_id"],
                "final_gene_id": final_gene,
                "final_transcript_id": final_transcript,
                "source_status": row["source_status"],
                "selection_locus_id": row["selection_locus_id"],
                "routes": row["routes"],
                "postqc_operation": row["operation"],
                "postqc_reason": row["reason"],
            }
        )
        structural_rows.append({**row, "final_gene_id": final_gene, "final_transcript_id": final_transcript})
    _write_table(transcript_id_map, TRANSCRIPT_MAP_COLUMNS, transcript_rows)
    _write_table(candidate_lineage, LINEAGE_COLUMNS, lineage_rows)
    structural_columns = tuple(decision_rows[0]) + ("final_gene_id", "final_transcript_id")
    _write_table(structural_decisions, structural_columns, structural_rows)

    edge_rows = _read_table(
        edge_audit,
        {"left_transcript_id", "right_transcript_id", "seqid", "span_overlap_bp", "disposition"},
    )
    flags = [
        {
            "flag_type": "OVERLAPPING_FRAME_OR_STRUCTURE_KEEP_SEPARATE",
            "severity": "INFO",
            "subject_a": row["left_transcript_id"],
            "subject_b": row["right_transcript_id"],
            "seqid": row["seqid"], "start": ".", "end": ".",
            "reason": row["disposition"], "release_blocking": 0,
        }
        for row in edge_rows
        if row["disposition"] in {
            "KEEP_SEPARATE_FRAME_OR_STRUCTURE",
            "KEEP_SEPARATE_ACTIVE_FUSION",
        }
    ]
    bridge_review_rows = _read_table(
        bridge_reviews,
        {
            "left_transcript_id", "right_transcript_id", "seqid",
            "intron_start", "intron_end", "intron_bp", "threshold_class",
            "status", "required_evidence",
        },
    )
    flags.extend(
        {
            "flag_type": row["status"],
            "severity": "REVIEW" if row["threshold_class"] == "ADDITIONAL_EVIDENCE" else "INFO",
            "subject_a": transcript_map.get(row["left_transcript_id"], row["left_transcript_id"]),
            "subject_b": transcript_map.get(row["right_transcript_id"], row["right_transcript_id"]),
            "seqid": row["seqid"],
            "start": row["intron_start"], "end": row["intron_end"],
            "reason": f"intron_bp={row['intron_bp']};required={row['required_evidence']}",
            "release_blocking": 0,
        }
        for row in bridge_review_rows
    )
    te_files = {}
    if classifier is not None:
        from geneinit.release.te import write_te_release
        flags.extend(write_te_release(annotation_gff3.parent, candidates_path, selection_audit_path,
            generated_te_path, classifier, genome, models_by_transcript, gene_map, transcript_map,
            lineage_rows, decisions, gene_prefix))
        for model in models_by_transcript.values():
            for flag in filter(None, te_annotations[model.candidate_id]["te_review_flags"].split(",")):
                flags.append(dict(flag_type=flag, severity="INFO", subject_a=transcript_map[model.transcript_id],
                    subject_b=".", seqid=model.seqid, start=model.start, end=model.end,
                    reason="TE annotation only; general candidate retained", release_blocking=0))
        te_files = {name: annotation_gff3.parent/name for name in
            ("te_annotation.gff3", "te_cds.fna", "te_proteins.faa", "te_candidate_audit.tsv", "te_summary.json")}
    else:
        for name in ("te_annotation.gff3", "te_cds.fna", "te_proteins.faa", "te_candidate_audit.tsv", "te_summary.json"):
            with atomic_text(annotation_gff3.parent/name) as handle:
                handle.write('{"status": "LEGACY_POLICY"}\n' if name.endswith('.json') else '##gff-version 3\n' if name.endswith('.gff3') else '')
    _write_table(review_flags, REVIEW_COLUMNS, flags)

    qc_values = [
        ("status", "PASS"),
        ("gene_count", len(ordered_genes)),
        ("transcript_count", len(ordered_transcripts)),
        ("canonical_transcript_count", len(ordered_genes)),
        ("cds_count", len(ordered_transcripts)),
        ("protein_count", len(ordered_transcripts)),
        ("longest_protein_count", len(longest_records)),
        ("retired_transcript_count", sum(row["postqc_gene_id"] in EMPTY for row in decision_rows)),
        ("nonblocking_review_flag_count", len(flags)),
        ("release_blocking_flag_count", 0),
        ("orf_failure_count", 0),
        ("hierarchy_failure_count", 0),
        ("sequence_mismatch_count", 0),
    ]
    _write_table(qc_tsv, ("metric", "value"), ({"metric": key, "value": value} for key, value in qc_values))
    qc_document = {key: value for key, value in qc_values}
    with atomic_text(qc_json) as handle:
        json.dump(qc_document, handle, indent=2, sort_keys=True)
        handle.write("\n")

    with atomic_text(frozen_policy_path) as handle:
        handle.write(policy_path.read_text(encoding="utf-8"))
    release_gate = "PASS" if policy_status == "PRODUCTION" else "PASS_MIGRATION_PENDING"
    status_document = {
        "schema_version": 1,
        "status": release_gate,
        "assembly_id": assembly_id,
        "policy_id": policy["policy_id"],
        "policy_status": policy_status,
        "structural_qc": "PASS",
        "production_publication_authorized": policy_status == "PRODUCTION",
    }
    with atomic_text(release_status) as handle:
        json.dump(status_document, handle, indent=2, sort_keys=True)
        handle.write("\n")
    release_files = {
        "annotation.gff3": annotation_gff3,
        "transcripts.fna": transcripts_fna,
        "cds.fna": cds_fna,
        "proteins.faa": proteins_faa,
        "longest_proteins.faa": longest_proteins_faa,
        "gene_id_map.tsv": gene_id_map,
        "transcript_id_map.tsv": transcript_id_map,
        "candidate_lineage.tsv": candidate_lineage,
        "structural_decisions.tsv": structural_decisions,
        "review_flags.tsv": review_flags,
        "structural_qc.tsv": qc_tsv,
        "structural_qc.json": qc_json,
        "release_status.json": release_status,
        "policy.json": frozen_policy_path,
    }
    release_files.update(te_files)
    checksums = {name: sha256_file(path) for name, path in release_files.items()}
    with atomic_text(checksums_path) as handle:
        for name in sorted(checksums):
            handle.write(f"{checksums[name]}  {name}\n")
    manifest = {
        **status_document,
        "gene_count": len(ordered_genes),
        "transcript_count": len(ordered_transcripts),
        "input_sha256": {
            "final_utr_gff3": sha256_file(final_utr_gff3),
            "transcript_decisions": sha256_file(transcript_decisions),
            "bridge_reviews": sha256_file(bridge_reviews),
            "genome": sha256_file(genome_path),
            "policy": sha256_file(policy_path),
        },
        "files": {name: {"sha256": checksums[name], "bytes": path.stat().st_size} for name, path in release_files.items()},
    }
    with atomic_text(manifest_path) as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    completion = (
        ("schema_version", 1), ("status", release_gate),
        ("assembly_id", assembly_id), ("gene_count", len(ordered_genes)),
        ("transcript_count", len(ordered_transcripts)),
        ("structural_qc_status", "PASS"),
        ("release_manifest_sha256", sha256_file(manifest_path)),
        ("checksums_sha256", sha256_file(checksums_path)),
    )
    _write_table(complete_path, ("metric", "value"), ({"metric": key, "value": value} for key, value in completion))
