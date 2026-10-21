"""Validate cached ranking and restore only source-proven exon strands."""

from collections import Counter, defaultdict
import csv
from pathlib import Path

from geneinit.utr.core import Feature, UtrError, _features, atomic_text, parse_gtf_attributes


def rna_features(path: Path) -> list[Feature]:
    # RNA may legitimately have unknown strand; coding-model parsing stays strict.
    result = []
    for number, line in enumerate(path.open(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise UtrError(f"{path}:{number} does not have 9 columns")
        start, end = int(fields[3]), int(fields[4])
        if start < 1 or end < start or fields[6] not in {"+", "-", "."}:
            raise UtrError(f"{path}:{number} has invalid RNA coordinates/strand")
        result.append(Feature(fields[0], fields[1], fields[2], start, end,
                              fields[6], fields[7], parse_gtf_attributes(fields[8])))
    return result


def validate_ranked_rna(ranked: Path, rna: Path, output: Path, audit: Path) -> dict[str, int]:
    """Retain ranking, IDs and coordinates; reject every change except the
    known InGenAnnot 0.0.20 exon-strand serialization defect. The original RNA
    must independently establish the transcript axis and complete exon chain.
    """
    source = rna_features(rna)
    transcripts = {}
    exons = defaultdict(Counter)
    for feature in source:
        if feature.feature_type in {"transcript", "mRNA"}:
            tid = feature.attributes.get("transcript_id")
            if not tid or tid in transcripts:
                raise UtrError("duplicate or missing source RNA transcript ID")
            transcripts[tid] = feature
    for feature in source:
        if feature.feature_type != "exon":
            continue
        tid = feature.attributes.get("transcript_id")
        parent = transcripts.get(tid)
        if parent is None or (feature.seqid, feature.strand, feature.attributes.get("gene_id")) != (
            parent.seqid, parent.strand, parent.attributes.get("gene_id")
        ):
            raise UtrError(f"source RNA exon hierarchy/axis mismatch: {tid}")
        exons[tid][(feature.start, feature.end)] += 1
    observed = _features(ranked, "gtf")
    seen = set()
    observed_exons = defaultdict(Counter)
    repairs = []
    excluded = set()
    exclusions = []
    for feature in observed:
        tid = feature.attributes.get("transcript_id")
        parent = transcripts.get(tid)
        if parent is None:
            raise UtrError(f"ranked RNA transcript absent from source: {tid}")
        if (feature.seqid, feature.attributes.get("gene_id")) != (
            parent.seqid, parent.attributes.get("gene_id")
        ):
            raise UtrError(f"ranked RNA identity/sequence changed: {tid}")
        if feature.feature_type in {"transcript", "mRNA"}:
            if tid in seen or (feature.start, feature.end) != (parent.start, parent.end) or (
                parent.strand != "." and feature.strand != parent.strand
            ):
                raise UtrError(f"ranked RNA transcript coordinates/strand changed: {tid}")
            seen.add(tid)
            if parent.strand == ".":
                # InGenAnnot exports unknown strand as +. Never use this invented
                # orientation to attach RNA evidence to a directional CDS model.
                excluded.add(tid)
                exclusions.append((tid, feature.seqid, feature.start, feature.end,
                                   feature.strand, ".", "EXCLUDE_UNKNOWN_SOURCE_STRAND"))
        elif feature.feature_type == "exon":
            coordinates = (feature.start, feature.end)
            observed_exons[tid][coordinates] += 1
            if coordinates not in exons[tid]:
                raise UtrError(f"ranked RNA exon absent from source: {tid}")
            if parent.strand != "." and feature.strand != parent.strand:
                if feature.strand not in {"+", "-"}:
                    raise UtrError(f"unknown ranked RNA exon strand: {tid}")
                repairs.append((tid, feature.seqid, feature.start, feature.end,
                                feature.strand, parent.strand, "RESTORE_SOURCE_EXON_STRAND"))
        else:
            raise UtrError(f"unexpected ranked RNA feature: {feature.feature_type}")
    if set(observed_exons) != seen or any(observed_exons[tid] != exons[tid] for tid in seen):
        raise UtrError("ranked RNA exon chain/membership changed")
    with atomic_text(output) as destination, ranked.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                tid = parse_gtf_attributes(fields[8])["transcript_id"]
                if tid in excluded:
                    continue
                if fields[2] == "exon":
                    if fields[6] != transcripts[tid].strand:
                        fields[6] = transcripts[tid].strand
                        line = "\t".join(fields) + "\n"
            destination.write(line)
    with atomic_text(audit, newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("transcript_id", "seqid", "start", "end", "ranked_strand", "source_strand", "action"))
        writer.writerows(repairs + exclusions)
    return {"exon_strands_restored": len(repairs), "unstranded_transcripts_excluded": len(excluded)}
