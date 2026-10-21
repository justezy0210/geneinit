"""Prepare bounded, nonredundant AUGUSTUS training and evaluation sets."""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import json
from pathlib import Path
import random
import re

from geneinit.models.candidate import CdsPart, audit_orf, reverse_complement
from geneinit.predictor.common import PredictorError


@dataclass(frozen=True)
class TrainingRecord:
    locus: str
    gene: str
    raw: str
    protein: str
    cds_bp: int
    genomic_bp: int
    exons: int
    gc: float
    flank: int
    ambiguous_bp: int
    reason: str
    chromosome: str
    start: int
    end: int


def read_training_records(path: Path) -> list[TrainingRecord]:
    """Read the complete, single-CDS GenBank records emitted by gff2gbSmallDNA."""
    records = []
    chunk = []
    seen_loci, seen_genes = set(), set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            chunk.append(line)
            if line.strip() != "//":
                continue
            raw = "".join(chunk)
            chunk = []
            locus_match = re.search(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp", raw)
            if locus_match is None or "ORIGIN" not in raw:
                raise PredictorError("malformed training GenBank record")
            locus, expected_length = locus_match.group(1), int(locus_match.group(2))
            head, origin = raw.split("ORIGIN", 1)
            sequence = re.sub(r"[^A-Za-z]", "", origin).upper()
            if len(sequence) != expected_length:
                raise PredictorError(f"GenBank sequence length mismatch: {locus}")
            genes = re.findall(r'/gene="([^"]+)"', head)
            if len(genes) != 1 or len(re.findall(r"^     CDS\s", head, re.M)) != 1:
                raise PredictorError(f"training record requires one gene/CDS: {locus}")
            gene = genes[0]
            if locus in seen_loci or gene in seen_genes:
                raise PredictorError(f"duplicate training locus or gene: {locus}, {gene}")
            seen_loci.add(locus)
            seen_genes.add(gene)
            location = re.sub(r"\s+", "", head.split("     CDS", 1)[1].split("/gene=", 1)[0])
            coords = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", location)]
            if not coords or "<" in location or ">" in location:
                raise PredictorError(f"unsupported partial/complex training CDS: {locus}")
            coords.sort()
            if any(a < 1 or b < a or b > len(sequence) for a, b in coords):
                raise PredictorError(f"training CDS out of bounds: {locus}")
            if any(b >= c for (_, b), (c, _) in zip(coords, coords[1:])):
                raise PredictorError(f"overlapping training CDS parts: {locus}")
            strand = "-" if "complement" in location else "+"
            parts, cumulative = [], 0
            for start, end in sorted(coords, reverse=strand == "-"):
                parts.append(CdsPart(start, end, (3 - cumulative % 3) % 3))
                cumulative += end - start + 1
            audit = audit_orf(locus, strand, parts, {locus: sequence})
            reason = "" if audit.status == "complete" else f"orf:{audit.status}"
            for (_, end), (start, _) in zip(coords, coords[1:]):
                intron = sequence[end:start - 1]
                if strand == "-":
                    intron = reverse_complement(intron)
                if len(intron) < 4 or (intron[:2], intron[-2:]) not in {
                    ("GT", "AG"), ("GC", "AG"), ("AT", "AC")
                }:
                    reason = reason or "noncanonical_splice"
            # AUGUSTUS etraining requires start/stop codons to fit inside their
            # respective terminal CDS exons, even when the spliced ORF is valid.
            boundary_issues = []
            if parts[0].length < 3:
                boundary_issues.append("unsupported_split_start_codon")
            if parts[-1].length < 3:
                boundary_issues.append("unsupported_split_stop_codon")
            reason = reason or ";".join(boundary_issues)
            genomic = re.fullmatch(r"(.+)_(\d+)-(\d+)", locus)
            if genomic is None:
                raise PredictorError(f"missing source coordinates in training locus: {locus}")
            chromosome, start, end = genomic.group(1), int(genomic.group(2)), int(genomic.group(3))
            if end - start + 1 != len(sequence):
                raise PredictorError(f"source interval length mismatch: {locus}")
            acgt = sum(sequence.count(x) for x in "ACGT")
            records.append(TrainingRecord(
                locus, gene, head + "ORIGIN" + origin.lower(),
                audit.protein_sequence.rstrip("*"), cumulative, len(sequence), len(coords),
                (sequence.count("G") + sequence.count("C")) / max(1, acgt),
                min(coords[0][0] - 1, len(sequence) - coords[-1][1]),
                len(sequence) - acgt, reason, chromosome, start, end,
            ))
    if any(line.strip() for line in chunk):
        raise PredictorError("unterminated training GenBank record")
    if not records:
        raise PredictorError("empty training GenBank")
    return records


def write_genbank(path: Path, records: list[TrainingRecord]) -> None:
    path.write_text("".join(record.raw for record in records), encoding="utf-8")


def write_proteins(path: Path, records: list[TrainingRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in sorted(records, key=lambda record: record.locus):
            handle.write(f">{record.locus}\n{record.protein}\n")


def training_errors(log: Path) -> dict[str, str]:
    """Preserve every distinct error per locus, including earlier fatal errors."""
    errors = defaultdict(list)
    for line in log.read_text(encoding="utf-8").splitlines():
        match = re.search(r"in sequence (\S+):\s*(.*)", line)
        if match:
            locus, message = match.groups()
            if message not in errors[locus]:
                errors[locus].append(message)
    return {locus: "; ".join(messages) for locus, messages in errors.items()}


def require_clean_training(log: Path) -> None:
    errors = training_errors(log)
    if errors:
        raise PredictorError(f"etraining reported {len(errors)} inconsistent loci; inspect {log}")


def require_matching_selection(observed: Path, expected: Path) -> None:
    """Do not train if regenerated assignments differ from the reviewed set."""
    if observed.read_bytes() != expected.read_bytes():
        raise PredictorError(
            f"training selection differs from the reviewed reference: {expected}; "
            f"inspect {observed} before training"
        )


def nonredundant_records(records: list[TrainingRecord], hits: Path,
                         identity: float = 80, coverage: float = 80):
    """One representative per connected protein-similarity / overlapping-DNA group."""
    by_id = {record.locus: record for record in records}
    parent = {locus: locus for locus in by_id}

    def find(locus):
        while parent[locus] != locus:
            parent[locus] = parent[parent[locus]]
            locus = parent[locus]
        return locus

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[max(left, right)] = min(left, right)

    exact = {}
    for record in records:
        if record.protein in exact:
            union(record.locus, exact[record.protein])
        exact[record.protein] = record.locus
    with hits.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) != 9:
                raise PredictorError(f"invalid DIAMOND training hit at line {number}")
            query, subject = fields[:2]
            if query not in by_id or subject not in by_id:
                raise PredictorError("DIAMOND training hit references unknown locus")
            pid, qs, qe, qlen, ss, se, slen = map(float, fields[2:])
            if qlen <= 0 or slen <= 0:
                raise PredictorError("nonpositive protein length in DIAMOND hit")
            if (pid >= identity and 100 * (abs(qe - qs) + 1) / qlen >= coverage
                    and 100 * (abs(se - ss) + 1) / slen >= coverage):
                union(query, subject)
    # Prevent overlapping genomic flanks from entering separate evaluation sets.
    previous = None
    for record in sorted(records, key=lambda r: (r.chromosome, r.start, r.end, r.locus)):
        if previous and record.chromosome == previous.chromosome and record.start <= previous.end:
            union(record.locus, previous.locus)
            if record.end > previous.end:
                previous = record
        else:
            previous = record
    groups = defaultdict(list)
    for record in records:
        groups[find(record.locus)].append(record)
    representatives = [min(group, key=lambda r: (
        r.ambiguous_bp / r.genomic_bp, -r.flank, -len(r.protein), r.locus
    )) for _, group in sorted(groups.items())]
    membership = {record.locus: group for group, members in groups.items() for record in members}
    return representatives, membership


def stratified_sample(records: list[TrainingRecord], count: int, seed: int,
                      gc_bounds: tuple[float, float]):
    """Seeded proportional quotas preserve GC, CDS length and exon-count strata."""
    if not 0 <= count <= len(records):
        raise PredictorError(f"cannot select {count} loci from {len(records)}")
    groups = defaultdict(list)
    for record in records:
        key = (bisect_right(gc_bounds, record.gc),
               bisect_right((600, 1500, 3000), record.cds_bp),
               bisect_right((2, 6, 11), record.exons))
        groups[key].append(record)
    total = len(records)
    quotas = {key: count * len(group) // max(1, total) for key, group in groups.items()}
    remaining = count - sum(quotas.values())
    ordered = sorted(groups, key=lambda key: (-(count * len(groups[key]) % max(1, total)), key))
    for key in ordered[:remaining]:
        quotas[key] += 1
    selected, rest = [], []
    for key, group in sorted(groups.items()):
        group.sort(key=lambda record: record.locus)
        digest = hashlib.sha256(f"{seed}:{key}".encode()).digest()
        random.Random(int.from_bytes(digest[:8], "big")).shuffle(group)
        selected.extend(group[:quotas[key]])
        rest.extend(group[quotas[key]:])
    return sorted(selected, key=lambda r: r.locus), sorted(rest, key=lambda r: r.locus)


def split_training(records: list[TrainingRecord], hits: Path, output: Path,
                   excluded: dict[str, str], *, max_loci: int = 6000,
                   optimize_loci: int = 500, holdout_loci: int = 500,
                   seed: int = 17, identity: float = 80, coverage: float = 80) -> dict:
    if not 0 < identity <= 100 or not 0 < coverage <= 100:
        raise PredictorError("training similarity identity/coverage must be in (0, 100]")
    if min(optimize_loci, holdout_loci) < 1 or max_loci <= optimize_loci + holdout_loci:
        raise PredictorError("training sizes must leave nonempty evaluation and onlytrain sets")
    eligible = [record for record in records if not record.reason and record.locus not in excluded]
    representatives, membership = nonredundant_records(eligible, hits, identity, coverage)
    selected_count = min(max_loci, len(representatives))
    if selected_count <= optimize_loci + holdout_loci:
        raise PredictorError(f"only {selected_count} independent loci survive; reduce requested split sizes")
    gcs = sorted(record.gc for record in representatives)
    bounds = (gcs[len(gcs) // 3], gcs[2 * len(gcs) // 3])
    holdout, rest = stratified_sample(representatives, holdout_loci, seed, bounds)
    optimize, rest = stratified_sample(rest, optimize_loci, seed + 1, bounds)
    onlytrain, unused = stratified_sample(rest, selected_count - optimize_loci - holdout_loci, seed + 2, bounds)
    splits = {"holdout": holdout, "optimize": optimize, "onlytrain": onlytrain}
    roles = {record.locus: role for role, group in splits.items() for record in group}
    selected_groups = [membership[locus] for locus in roles]
    if len(set(selected_groups)) != len(selected_groups):
        raise PredictorError("training/evaluation similarity group leakage")
    output.mkdir(parents=True, exist_ok=True)
    for role, group in splits.items():
        write_genbank(output / f"{role}.gb", group)
    write_genbank(output / "train.gb", sorted(onlytrain + optimize, key=lambda record: record.locus))
    with (output / "selection.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("locus", "gene", "group", "role", "reason", "cds_bp", "genomic_bp", "exons", "gc"))
        representative_ids = {record.locus for record in representatives}
        for record in sorted(records, key=lambda record: record.locus):
            reason = record.reason or excluded.get(record.locus, "")
            role = roles.get(record.locus, "excluded" if reason else
                             "unused" if record.locus in representative_ids else "redundant")
            writer.writerow((record.locus, record.gene, membership.get(record.locus, ""), role,
                             reason, record.cds_bp, record.genomic_bp, record.exons, f"{record.gc:.6f}"))
    def describe(group):
        return {"loci": len(group), "genomic_bp": sum(r.genomic_bp for r in group),
                "cds_bp": sum(r.cds_bp for r in group),
                "single_exon_loci": sum(r.exons == 1 for r in group),
                "gc_mean": sum(r.gc for r in group) / max(1, len(group))}
    summary = {
        "schema_version": 1, "method": "quality-protein-spatial-stratified-v1",
        "seed": seed, "max_loci": max_loci, "identity": identity, "coverage": coverage,
        "gc_tertile_boundaries": bounds,
        "input": describe(records), "eligible": describe(eligible),
        "nonredundant": describe(representatives), "splits": {k: describe(v) for k, v in splits.items()},
        "unused_nonredundant_loci": len(unused),
        "orf_splice_exclusions": dict(Counter(r.reason for r in records if r.reason)),
        "etraining_excluded_loci": len(excluded), "cross_split_group_overlap": 0,
    }
    (output / "selection.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def evaluation_metrics(path: Path) -> dict:
    metrics = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*(nucleotide|exon|gene) level\s*\|", line)
        if match:
            columns = [column.strip() for column in line.split("|") if column.strip()]
            try:
                sensitivity, specificity = map(float, columns[-2:])
            except ValueError as error:
                raise PredictorError(f"invalid held-out evaluation row: {line}") from error
            if not (0 <= sensitivity <= 1 and 0 <= specificity <= 1):
                raise PredictorError(f"invalid held-out evaluation accuracy: {line}")
            metrics[match.group(1)] = {"sensitivity": sensitivity, "specificity": specificity,
                "f1": 2 * sensitivity * specificity / (sensitivity + specificity) if sensitivity + specificity else 0.0}
    if set(metrics) != {"nucleotide", "exon", "gene"}:
        raise PredictorError(f"incomplete held-out AUGUSTUS evaluation: {path}")
    return metrics
