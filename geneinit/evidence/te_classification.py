"""Confirmed-only TE decisions from versioned, sequence-bound evidence.

Search programs produce the evidence bundle described in docs/te-evidence.md.
This module never treats a missing search as a negative result. Protein results
are cached by SHA256; intact structures and reviews remain candidate-specific.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path

POLICY_ID = "geneinit_te_confirmed_only_v2"
TE_COLUMNS = (
    "te_policy_id", "protein_sha256", "te_evidence_sha256", "te_signature_version",
    "te_exception_version", "host_state", "te_direct_status", "te_host_status",
    "te_direct_pass", "te_intact_pass", "te_architecture_pass",
    "te_separation_eligible", "automatic_class", "review_result", "effective_class",
    "candidate_pool", "te_review_flags", "te_selected",
)


class TeEvidenceError(ValueError):
    pass


def confirmed_only(policy):
    mode = policy.get("te_policy_id")
    if mode not in (None, POLICY_ID):
        raise TeEvidenceError(f"unsupported te_policy_id: {mode}")
    return mode == POLICY_ID


def repeat_neutral(record):
    return record.row.get("te_policy_id") == POLICY_ID


def protein_sha256(sequence):
    return hashlib.sha256(sequence.upper().removesuffix("*").encode()).hexdigest()


def model_identity(model):
    parts = sorted(model.cds, key=lambda p: (p.start, p.end))
    return dict(candidate_id=model.candidate_id, seqid=model.seqid,
                strand=model.strand, start=str(parts[0].start), end=str(parts[-1].end),
                cds_signature=",".join(f"{p.start}-{p.end}:{p.phase}" for p in parts))


def partition_pools(loci, edges):
    """Retain the raw locus IDs in audit, isolate TE before any graph operation."""
    from collections import defaultdict
    grouped = defaultdict(list)
    membership = {}
    for locus, members in loci.items():
        for record in members:
            pool = record.row.get("candidate_pool", "GENERAL")
            key = (locus, pool)
            grouped[key].append(record)
            membership[record.key] = key
    pool_edges = [replace(edge, locus_id=membership[edge.left]) for edge in edges
                  if membership[edge.left] == membership[edge.right]]
    return grouped, pool_edges


def _number(value):
    try:
        result = float(value)
    except (ValueError, TypeError) as error:
        raise TeEvidenceError(f"invalid numeric TE evidence: {value!r}") from error
    if not math.isfinite(result):
        raise TeEvidenceError("non-finite TE evidence")
    return result


def review_snapshot(bundle, context):
    payload = {"protein": bundle.get("proteins", {}).get(context.get("protein_sha256")),
               "signatures": bundle.get("signatures"),
               "signature_version": bundle.get("signature_version"),
               "exception_version": bundle.get("exception_version"),
               "context": {key: value for key, value in context.items() if key != "review"}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


class TeClassifier:
    def __init__(self, path: Path | None = None, assembly_id: str | None = None):
        self.bundle = {}
        self.digest = "NOT_RUN"
        self.cache = {}
        if path is not None:
            raw = path.read_bytes()  # explicitly configured missing file is an error
            self.digest = hashlib.sha256(raw).hexdigest()
            self.bundle = json.loads(raw)
            if self.bundle.get("schema") != "geneinit_te_evidence_v2":
                raise TeEvidenceError("unsupported TE evidence schema")
            if not assembly_id or self.bundle.get("assembly_id") != assembly_id:
                raise TeEvidenceError("TE evidence assembly mismatch")
            for field in ("signature_version", "exception_version"):
                if not self.bundle.get(field):
                    raise TeEvidenceError(f"TE evidence lacks {field}")
            # Validate every configured search, including unused query records.
            for digest, result in self.bundle.get("proteins", {}).items():
                if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                    raise TeEvidenceError("invalid protein SHA256")
                self._search(result.get("direct", {}), digest)
                self._search(result.get("host", {}), digest)

    def _search(self, search, digest):
        state = search.get("status", "NOT_RUN")
        if state == "FAILED":
            raise TeEvidenceError("BLOCKED_INPUT: configured TE/host search failed")
        if state not in {"NOT_RUN", "COMPLETE_HIT", "COMPLETE_NO_HIT"}:
            raise TeEvidenceError(f"invalid search status: {state}")
        hits = search.get("hits", [])
        if state == "NOT_RUN":
            if hits:
                raise TeEvidenceError("uncompleted search contains hits")
            return state, []
        if search.get("query_sha256") != digest:
            raise TeEvidenceError("search query SHA256 mismatch")
        if not all(search.get(key) for key in ("tool_version", "database_sha256", "options", "receipt")):
            raise TeEvidenceError("completed search lacks provenance/receipt")
        if (state == "COMPLETE_HIT") != bool(hits):
            raise TeEvidenceError("search status and hits disagree")
        return state, hits

    def _protein(self, digest):
        if digest in self.cache:
            return self.cache[digest]
        result = self.bundle.get("proteins", {}).get(digest, {})
        ds, direct = self._search(result.get("direct", {}), digest)
        hs, host = self._search(result.get("host", {}), digest)
        supported = any(
            (h.get("kind") == "BUSCO" and h.get("status") in {"Complete", "Duplicated"})
            or (h.get("kind") == "CELLULAR_SWISSPROT" and h.get("accession")
                and _number(h["identity"]) >= 40 and _number(h["query_coverage"]) >= .8
                and _number(h["subject_coverage"]) >= .7 and _number(h["alignment_aa"]) >= 150
                and _number(h["evalue"]) <= 1e-20)
            or h.get("kind") == "CURATED_DOMESTICATED_GENE"
            for h in host
        )
        hstate = "SUPPORTED" if supported else "NOT_ASSESSED" if hs == "NOT_RUN" else "NO_STRONG_MATCH"
        matches = {}
        for family, signature in self.bundle.get("signatures", {}).items():
            if signature.get("validated") is not True:
                continue
            required = signature.get("domains", [])
            if not required or not signature.get("intact_features") or not signature.get("validation_reference"):
                raise TeEvidenceError(f"incomplete validated signature: {family}")
            positions = []
            for domain in required:
                passing = [hit for hit in direct
                           if hit.get("kind") == "DIRECT_DOMAIN" and hit.get("family") == family
                           and hit.get("domain") == domain["id"]
                           and _number(hit["score"]) >= _number(domain["min_score"])
                           and _number(hit["coverage"]) >= _number(domain["min_coverage"])]
                if not passing:
                    break
                positions.append(sorted({int(hit["start_aa"]) for hit in passing}))
            if len(positions) == len(required):
                matches[family] = True
                if signature.get("ordered", True):
                    # The earliest usable occurrence leaves every later hit
                    # available, including repeated domains in the signature.
                    previous = None
                    for starts in positions:
                        position = next((start for start in starts
                                         if previous is None or start > previous), None)
                        if position is None:
                            matches[family] = False
                            break
                        previous = position
        suspect = bool(result.get("eggnog_te_terms")) or any(
            hit.get("kind") == "DIRECT_DOMAIN" or
            (hit.get("kind") == "TE_PROTEIN" and _number(hit["evalue"]) <= 1e-10
             and _number(hit["identity"]) >= 30 and _number(hit["query_coverage"]) >= .5
             and _number(hit["alignment_aa"]) >= 80) for hit in direct
        )
        self.cache[digest] = ds, hs, hstate, matches, suspect
        return self.cache[digest]

    def classify(self, row, digest):
        ds, hs, host, matches, suspect = self._protein(digest)
        direct_pass = bool(matches)
        architecture = any(matches.values())
        intact = False
        context = self.bundle.get("candidates", {}).get(row["candidate_id"])
        review = "NOT_REVIEWED"
        flags = []
        if context is not None:
            for key in ("seqid", "strand", "start", "end", "cds_signature"):
                if str(context.get(key)) != str(row[key]):
                    raise TeEvidenceError(f"candidate TE coordinate/CDS mismatch: {row['candidate_id']} {key}")
            if context.get("protein_sha256") != digest:
                raise TeEvidenceError(f"candidate TE protein mismatch: {row['candidate_id']}")
            for structure in context.get("structures", []):
                family = structure.get("family")
                sig = self.bundle.get("signatures", {}).get(family, {})
                if not matches.get(family):
                    continue
                intact |= (
                    structure.get("status") == "INTACT"
                    and structure.get("seqid") == row["seqid"]
                    and structure.get("strand") == row["strand"]
                    and int(structure["start"]) <= int(row["start"])
                    and int(structure["end"]) >= int(row["end"])
                    and bool(structure.get("source_sha256"))
                    and bool(structure.get("structure_id"))
                    and set(sig["intact_features"]) <= set(structure.get("features", []))
                )
            reviewed = context.get("review", {})
            if reviewed:
                # The context already binds the review to the exact sequence/locus.
                snapshot = review_snapshot(self.bundle, context)
                if reviewed.get("evidence_snapshot_sha256") != snapshot:
                    flags.append("STALE_REVIEW_IGNORED")
                else:
                    review = reviewed.get("result", "UNRESOLVED")
                    if review not in {"HOST", "TE", "UNRESOLVED"} or not reviewed.get("reason"):
                        raise TeEvidenceError("invalid TE review result/reason")
        eligible = direct_pass and architecture and intact and host == "NO_STRONG_MATCH"
        if eligible:
            automatic = "TE_HIGH_CONFIDENCE"
        elif (suspect or direct_pass) and host == "SUPPORTED":
            automatic = "HOST_TE_CONFLICT"
        elif direct_pass and host == "NOT_ASSESSED":
            automatic = "HOST_SCREEN_PENDING"
        elif direct_pass:
            automatic = "TE_INCOMPLETE_EVIDENCE"
        elif suspect:
            automatic = "TE_SUSPECT"
        else:
            automatic = "HOST_SUPPORTED" if host == "SUPPORTED" else "FUNCTION_UNRESOLVED"
        # HOST/UNRESOLVED reviews retain GENERAL. A TE opinion never bypasses evidence.
        if review in {"HOST", "UNRESOLVED"}:
            eligible = False
        effective = "REVIEWED_TE" if eligible and review == "TE" else automatic
        if not eligible and automatic == "TE_HIGH_CONFIDENCE":
            effective = "REVIEWED_HOST" if review == "HOST" else "FUNCTION_UNRESOLVED"
        if effective not in {"TE_HIGH_CONFIDENCE", "REVIEWED_TE", "HOST_SUPPORTED"}:
            flags.append(effective)
        if float(row.get("repeat_cds_fraction", 0) or 0) > 0:
            flags.append("CLASSIFIED_TE_OVERLAP")
        return dict(te_policy_id=POLICY_ID, protein_sha256=digest,
                    te_evidence_sha256=self.digest,
                    te_signature_version=self.bundle.get("signature_version", "NOT_ASSESSED"),
                    te_exception_version=self.bundle.get("exception_version", "NOT_ASSESSED"),
                    host_state=host, te_direct_status=ds, te_host_status=hs,
                    te_direct_pass=str(int(direct_pass)), te_intact_pass=str(int(intact)),
                    te_architecture_pass=str(int(architecture)), te_separation_eligible=str(int(eligible)),
                    automatic_class=automatic, review_result=review, effective_class=effective,
                    candidate_pool="TE" if eligible else "GENERAL", te_review_flags=",".join(flags),
                    te_selected="0")
