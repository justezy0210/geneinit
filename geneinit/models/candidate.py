"""Normalize predictor GTF records into a deterministic candidate universe.

The stable structure identifier intentionally excludes predictor names and the
predictors' own gene/transcript identifiers.  Those values remain attached as
lineage metadata, while identical CDS structures proposed by multiple sources
collapse into one candidate.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
import csv
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Iterable, Iterator, Mapping, Sequence, TextIO


VALID_STRANDS = {"+", "-"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
DNA_COMPLEMENT = str.maketrans(
    "ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
    "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn",
)

SOURCE_COLUMNS = (
    "source_model_id",
    "structure_id",
    "predictor",
    "source_gene_id",
    "source_transcript_id",
    "seqid",
    "strand",
    "start",
    "end",
    "cds_part_count",
    "exon_part_count",
    "raw_cds_length",
    "translated_cds_length",
    "first_phase",
    "phase_consistent",
    "start_codon",
    "terminal_stop",
    "internal_stop_count",
    "ambiguous_codon_count",
    "orf_status",
    "cds_signature",
)

UNIVERSE_COLUMNS = (
    "candidate_id",
    "structure_id",
    "seqid",
    "strand",
    "start",
    "end",
    "cds_part_count",
    "exon_part_count_max",
    "raw_cds_length",
    "translated_cds_length",
    "first_phase",
    "phase_consistent",
    "start_codon",
    "terminal_stop",
    "internal_stop_count",
    "ambiguous_codon_count",
    "orf_status",
    "predictor_count",
    "predictors",
    "source_model_count",
    "source_model_ids",
    "source_gene_ids",
    "source_transcript_ids",
    "cds_signature",
)


class CandidateError(RuntimeError):
    """Raised when predictor models cannot be normalized safely."""


@dataclass(frozen=True, order=True, slots=True)
class CdsPart:
    start: int
    end: int
    phase: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True, slots=True)
class OrfAudit:
    raw_cds_length: int
    translated_cds_length: int
    first_phase: int
    phase_consistent: bool
    start_codon: str
    terminal_stop: str
    internal_stop_count: int
    ambiguous_codon_count: int
    status: str
    coding_sequence: str
    protein_sequence: str


@dataclass(frozen=True, slots=True)
class CandidateModel:
    predictor: str
    source_gene_id: str
    source_transcript_id: str
    seqid: str
    strand: str
    cds_parts: tuple[CdsPart, ...]
    exon_parts: tuple[tuple[int, int], ...]
    structure_id: str
    source_model_id: str
    orf: OrfAudit

    @property
    def start(self) -> int:
        return min(part.start for part in self.cds_parts)

    @property
    def end(self) -> int:
        return max(part.end for part in self.cds_parts)

    @property
    def cds_signature(self) -> str:
        return cds_signature(self.cds_parts)


@dataclass(frozen=True, slots=True)
class CandidateUniverse:
    structure_id: str
    candidate_id: str
    models: tuple[CandidateModel, ...]

    @property
    def representative(self) -> CandidateModel:
        return self.models[0]


@dataclass(slots=True)
class _TranscriptBuilder:
    source_gene_id: str
    source_transcript_id: str
    seqid: str
    strand: str
    cds_parts: list[CdsPart] = field(default_factory=list)
    exon_parts: list[tuple[int, int]] = field(default_factory=list)


def parse_attributes(text: str) -> tuple[dict[str, str], str | None]:
    """Parse GTF/GFF-like attributes and retain a bare AUGUSTUS identifier."""

    attributes: dict[str, str] = {}
    bare: str | None = None
    for raw_token in text.strip().strip(";").split(";"):
        token = raw_token.strip()
        if not token:
            continue
        if "=" in token and not re.search(r"\s", token.split("=", 1)[0]):
            key, value = token.split("=", 1)
            attributes[key.strip()] = value.strip().strip('"')
            continue
        match = re.fullmatch(r"([^\s]+)\s+\"([^\"]*)\"", token)
        if match:
            attributes[match.group(1)] = match.group(2)
            continue
        match = re.fullmatch(r"([^\s]+)\s+(.+)", token)
        if match:
            attributes[match.group(1)] = match.group(2).strip().strip('"')
            continue
        if bare is not None and bare != token:
            raise CandidateError(f"multiple bare identifiers in attributes: {text!r}")
        bare = token.strip('"')
    return attributes, bare


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split(maxsplit=1)[0]
                if not current:
                    raise CandidateError(f"{path}:{number} has an empty FASTA ID")
                if current in sequences:
                    raise CandidateError(f"duplicate FASTA ID: {current}")
                sequences[current] = []
                continue
            if current is None:
                raise CandidateError(f"{path}:{number} sequence precedes FASTA header")
            sequence = re.sub(r"\s+", "", line).upper()
            if not re.fullmatch(r"[A-Z*.-]+", sequence):
                raise CandidateError(f"{path}:{number} contains invalid FASTA characters")
            sequences[current].append(sequence)
    if not sequences:
        raise CandidateError(f"no FASTA records: {path}")
    return {seqid: "".join(parts) for seqid, parts in sequences.items()}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def transcript_order(parts: Sequence[CdsPart], strand: str) -> tuple[CdsPart, ...]:
    return tuple(sorted(parts, key=lambda part: part.start, reverse=strand == "-"))


def cds_signature(parts: Sequence[CdsPart]) -> str:
    return ",".join(f"{part.start}-{part.end}:{part.phase}" for part in parts)


def structure_identifier(seqid: str, strand: str, parts: Sequence[CdsPart]) -> str:
    ordered_genomic = tuple(sorted(parts, key=lambda part: (part.start, part.end, part.phase)))
    payload = {
        "schema": 1,
        "seqid": seqid,
        "strand": strand,
        "cds": [[part.start, part.end, part.phase] for part in ordered_genomic],
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return f"GSTR-{digest[:24]}"


def source_model_identifier(predictor: str, structure_id: str, transcript_id: str) -> str:
    payload = f"1\0{predictor}\0{structure_id}\0{transcript_id}".encode()
    return f"GSRC-{hashlib.sha256(payload).hexdigest()[:24]}"


def _expected_next_phase(cumulative_length: int) -> int:
    return (3 - (cumulative_length % 3)) % 3


def audit_orf(
    seqid: str,
    strand: str,
    parts: Sequence[CdsPart],
    genome: Mapping[str, str],
) -> OrfAudit:
    if seqid not in genome:
        raise CandidateError(f"model seqid absent from genome: {seqid}")
    sequence = genome[seqid]
    ordered = transcript_order(parts, strand)
    fragments: list[str] = []
    phase_consistent = True
    cumulative = 0
    for index, part in enumerate(ordered):
        if part.start < 1 or part.end < part.start or part.end > len(sequence):
            raise CandidateError(
                f"CDS coordinate out of bounds: {seqid}:{part.start}-{part.end}"
            )
        fragment = sequence[part.start - 1 : part.end]
        if strand == "-":
            fragment = reverse_complement(fragment)
        fragments.append(fragment)
        if index > 0 and part.phase != _expected_next_phase(cumulative):
            phase_consistent = False
        cumulative += part.length

    raw = "".join(fragments)
    first_phase = ordered[0].phase
    coding = raw[first_phase:]
    complete_codons = len(coding) // 3
    codons = [coding[index : index + 3] for index in range(0, complete_codons * 3, 3)]
    start = codons[0] if codons else ""
    terminal = codons[-1] if codons else ""
    internal = sum(codon in STOP_CODONS for codon in codons[:-1])
    ambiguous = sum(not re.fullmatch(r"[ACGT]{3}", codon) for codon in codons)
    translation_table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }
    protein = "".join(translation_table.get(codon, "X") for codon in codons)

    issues: list[str] = []
    if first_phase != 0:
        issues.append("partial_5prime")
    if len(coding) % 3:
        issues.append("length_not_mod3")
    if start != "ATG" and first_phase == 0:
        issues.append("missing_start")
    if terminal not in STOP_CODONS:
        issues.append("missing_stop")
    if internal:
        issues.append("internal_stop")
    if ambiguous:
        issues.append("ambiguous_codon")
    if not phase_consistent:
        issues.append("phase_inconsistent")
    status = "complete" if not issues else ";".join(issues)
    return OrfAudit(
        raw_cds_length=len(raw),
        translated_cds_length=len(coding),
        first_phase=first_phase,
        phase_consistent=phase_consistent,
        start_codon=start,
        terminal_stop=terminal,
        internal_stop_count=internal,
        ambiguous_codon_count=ambiguous,
        status=status,
        coding_sequence=coding,
        protein_sequence=protein,
    )


def _resolve_transcript_id(
    feature: str, attributes: Mapping[str, str], bare: str | None
) -> str | None:
    resolved = (
        attributes.get("transcript_id")
        or attributes.get("transcript")
        or (attributes.get("ID") if feature in {"transcript", "mRNA"} else None)
        or (attributes.get("Parent") if feature in {"CDS", "exon"} else None)
        or (bare if feature in {"transcript", "mRNA"} else None)
    )
    if resolved and "," in resolved:
        raise CandidateError(
            f"multi-parent {feature} records are not supported: {resolved}"
        )
    return resolved


def _resolve_gene_id(
    feature: str,
    attributes: Mapping[str, str],
    bare: str | None,
    transcript_id: str | None,
) -> str:
    return (
        attributes.get("gene_id")
        or attributes.get("gene")
        or attributes.get("Parent")
        or (attributes.get("ID") if feature == "gene" else None)
        or (bare if feature == "gene" else None)
        or (transcript_id.rsplit(".", 1)[0] if transcript_id and "." in transcript_id else "")
        or transcript_id
        or "unknown_gene"
    )


def parse_predictor_gtf(
    path: Path,
    predictor: str,
    genome: Mapping[str, str],
) -> list[CandidateModel]:
    """Parse AUGUSTUS-like or GeneMark-like GTF into audited CDS models."""

    builders: dict[str, _TranscriptBuilder] = {}
    parent_genes: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise CandidateError(f"{path}:{number} does not have 9 columns")
            seqid, _source, feature, start_text, end_text, _score, strand, phase_text, attr_text = fields
            if strand not in VALID_STRANDS:
                if feature in {"gene", "transcript", "mRNA", "CDS", "exon"}:
                    raise CandidateError(f"{path}:{number} has unsupported strand {strand!r}")
                continue
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise CandidateError(f"{path}:{number} has non-integer coordinates") from error
            if start < 1 or end < start:
                raise CandidateError(f"{path}:{number} has invalid coordinates")
            if seqid not in genome or end > len(genome[seqid]):
                raise CandidateError(f"{path}:{number} coordinate is outside genome")
            attributes, bare = parse_attributes(attr_text)
            transcript_id = _resolve_transcript_id(feature, attributes, bare)
            gene_id = _resolve_gene_id(feature, attributes, bare, transcript_id)

            if feature == "gene":
                gene_key = attributes.get("ID") or bare or gene_id
                if gene_key:
                    parent_genes[gene_key] = gene_id
                continue
            if feature not in {"transcript", "mRNA", "CDS", "exon"}:
                continue
            if transcript_id is None:
                raise CandidateError(f"{path}:{number} {feature} lacks transcript identifier")
            parent = attributes.get("Parent")
            if parent and parent in parent_genes:
                gene_id = parent_genes[parent]
            existing = builders.get(transcript_id)
            if existing is None:
                existing = _TranscriptBuilder(gene_id, transcript_id, seqid, strand)
                builders[transcript_id] = existing
            elif existing.seqid != seqid or existing.strand != strand:
                raise CandidateError(
                    f"{path}:{number} transcript ID reused across seqid/strand: {transcript_id}"
                )
            elif existing.source_gene_id in {"unknown_gene", transcript_id} and gene_id:
                existing.source_gene_id = gene_id
            if feature == "CDS":
                if phase_text not in {"0", "1", "2"}:
                    raise CandidateError(f"{path}:{number} CDS has invalid phase {phase_text!r}")
                existing.cds_parts.append(CdsPart(start, end, int(phase_text)))
            elif feature == "exon":
                existing.exon_parts.append((start, end))

    models: list[CandidateModel] = []
    seen_source_ids: set[str] = set()
    for transcript_id in sorted(builders):
        builder = builders[transcript_id]
        if not builder.cds_parts:
            continue
        genomic_parts = tuple(sorted(builder.cds_parts))
        if len(set(genomic_parts)) != len(genomic_parts):
            raise CandidateError(f"duplicate CDS part in transcript {transcript_id}")
        for left, right in zip(genomic_parts, genomic_parts[1:]):
            if right.start <= left.end:
                raise CandidateError(f"overlapping CDS parts in transcript {transcript_id}")
        exons = tuple(sorted(set(builder.exon_parts)))
        structure_id = structure_identifier(builder.seqid, builder.strand, genomic_parts)
        source_id = source_model_identifier(predictor, structure_id, transcript_id)
        if source_id in seen_source_ids:
            raise CandidateError(f"source model identifier collision: {source_id}")
        seen_source_ids.add(source_id)
        models.append(
            CandidateModel(
                predictor=predictor,
                source_gene_id=builder.source_gene_id,
                source_transcript_id=transcript_id,
                seqid=builder.seqid,
                strand=builder.strand,
                cds_parts=genomic_parts,
                exon_parts=exons,
                structure_id=structure_id,
                source_model_id=source_id,
                orf=audit_orf(builder.seqid, builder.strand, genomic_parts, genome),
            )
        )
    if not models:
        raise CandidateError(f"no CDS-bearing transcript models found: {path}")
    return sorted(models, key=model_sort_key)


def model_sort_key(model: CandidateModel) -> tuple[object, ...]:
    return (
        model.seqid,
        model.start,
        model.end,
        model.strand,
        model.structure_id,
        model.predictor,
        model.source_transcript_id,
    )


def build_candidate_universe(models: Iterable[CandidateModel]) -> list[CandidateUniverse]:
    grouped: dict[str, list[CandidateModel]] = defaultdict(list)
    signatures: dict[str, tuple[str, str, tuple[CdsPart, ...]]] = {}
    for model in models:
        signature = (model.seqid, model.strand, model.cds_parts)
        previous = signatures.setdefault(model.structure_id, signature)
        if previous != signature:
            raise CandidateError(f"structure identifier collision: {model.structure_id}")
        grouped[model.structure_id].append(model)
    universe: list[CandidateUniverse] = []
    for structure_id, members in grouped.items():
        ordered = tuple(sorted(members, key=model_sort_key))
        universe.append(
            CandidateUniverse(
                structure_id=structure_id,
                candidate_id=f"GCAN-{structure_id.removeprefix('GSTR-')}",
                models=ordered,
            )
        )
    return sorted(
        universe,
        key=lambda item: (
            item.representative.seqid,
            item.representative.start,
            item.representative.end,
            item.representative.strand,
            item.structure_id,
        ),
    )


def _atomic_text(path: Path) -> tuple[Path, TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    return temporary, temporary.open("w", encoding="utf-8", newline="")


def _source_row(model: CandidateModel) -> tuple[object, ...]:
    audit = model.orf
    return (
        model.source_model_id, model.structure_id, model.predictor,
        model.source_gene_id, model.source_transcript_id, model.seqid, model.strand,
        model.start, model.end, len(model.cds_parts), len(model.exon_parts),
        audit.raw_cds_length, audit.translated_cds_length, audit.first_phase,
        int(audit.phase_consistent), audit.start_codon, audit.terminal_stop,
        audit.internal_stop_count, audit.ambiguous_codon_count, audit.status,
        model.cds_signature,
    )


def _universe_row(candidate: CandidateUniverse) -> tuple[object, ...]:
    representative = candidate.representative
    audit = representative.orf
    models = candidate.models
    predictors = sorted({model.predictor for model in models})
    return (
        candidate.candidate_id, candidate.structure_id, representative.seqid,
        representative.strand, representative.start, representative.end,
        len(representative.cds_parts), max(len(model.exon_parts) for model in models),
        audit.raw_cds_length, audit.translated_cds_length, audit.first_phase,
        int(audit.phase_consistent), audit.start_codon, audit.terminal_stop,
        audit.internal_stop_count, audit.ambiguous_codon_count, audit.status,
        len(predictors), ",".join(predictors), len(models),
        ",".join(model.source_model_id for model in models),
        ",".join(sorted({model.source_gene_id for model in models})),
        ",".join(model.source_transcript_id for model in models),
        representative.cds_signature,
    )


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    temporary, handle = _atomic_text(path)
    try:
        with handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_source_models(path: Path, models: Sequence[CandidateModel]) -> None:
    write_tsv(path, SOURCE_COLUMNS, (_source_row(model) for model in models))


def write_universe(path: Path, universe: Sequence[CandidateUniverse]) -> None:
    write_tsv(path, UNIVERSE_COLUMNS, (_universe_row(candidate) for candidate in universe))


def fasta_records(universe: Sequence[CandidateUniverse], kind: str) -> Iterator[tuple[str, str]]:
    if kind not in {"cds", "protein"}:
        raise ValueError(f"unsupported FASTA kind: {kind}")
    for candidate in universe:
        audit = candidate.representative.orf
        sequence = audit.coding_sequence if kind == "cds" else audit.protein_sequence
        yield candidate.candidate_id, sequence


def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> None:
    temporary, handle = _atomic_text(path)
    try:
        with handle:
            for identifier, sequence in records:
                handle.write(f">{identifier}\n")
                for offset in range(0, len(sequence), width):
                    handle.write(sequence[offset : offset + width] + "\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_candidate_gff3(path: Path, universe: Sequence[CandidateUniverse]) -> None:
    temporary, handle = _atomic_text(path)
    try:
        with handle:
            handle.write("##gff-version 3\n")
            for candidate in universe:
                model = candidate.representative
                candidate_id = candidate.candidate_id
                predictors = ",".join(sorted({item.predictor for item in candidate.models}))
                handle.write(
                    f"{model.seqid}\tGeneInIt\tmRNA\t{model.start}\t{model.end}\t.\t"
                    f"{model.strand}\t.\tID={candidate_id};predictors={predictors};"
                    f"structure_id={candidate.structure_id}\n"
                )
                for index, part in enumerate(transcript_order(model.cds_parts, model.strand), 1):
                    handle.write(
                        f"{model.seqid}\tGeneInIt\tCDS\t{part.start}\t{part.end}\t.\t"
                        f"{model.strand}\t{part.phase}\tID={candidate_id}.cds{index};"
                        f"Parent={candidate_id}\n"
                    )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_summary(
    path: Path,
    models: Sequence[CandidateModel],
    universe: Sequence[CandidateUniverse],
) -> None:
    status_counts: dict[str, int] = defaultdict(int)
    for candidate in universe:
        status_counts[candidate.representative.orf.status] += 1
    rows: list[tuple[str, object]] = [
        ("source_model_count", len(models)),
        ("candidate_structure_count", len(universe)),
        ("exact_structure_duplicates_collapsed", len(models) - len(universe)),
        ("multi_predictor_candidate_count", sum(len({m.predictor for m in c.models}) > 1 for c in universe)),
        ("complete_orf_candidate_count", status_counts.get("complete", 0)),
    ]
    rows.extend((f"orf_status::{status}", count) for status, count in sorted(status_counts.items()))
    write_tsv(path, ("metric", "value"), rows)
