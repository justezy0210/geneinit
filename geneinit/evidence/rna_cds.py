"""CDS/RNA compatibility, including introns crossing coding boundaries.

These are structural observations, not proof of translation. A complete
compatible RNA preserves the candidate's entire CDS chain; ORF/phase validity
remains the responsibility of the existing genomic ORF audit.
"""

from itertools import chain, islice
import json

from geneinit.evidence.metrics import interval_bp, overlap_bp


RNA_CDS_BASE_COLUMNS = (
    "rna_cds_boundary_version", "rna_cds_status",
    "rna_cds_support_transcripts_json", "rna_cds_conflict_transcripts_json",
    "rna_cds_recurrent_conflicts_json", "rna_cds_conflict_bp",
    "rna_cds_start_conflict", "rna_cds_stop_conflict",
)
RNA_CDS_COLUMNS = RNA_CDS_BASE_COLUMNS + (
    "rna_cds_chain_support_transcripts_json", "rna_cds_self_strong_chain",
    "rna_cds_countered_conflicts_json", "rna_cds_unopposed_conflicts_json",
    "rna_cds_unopposed_conflict_fraction",
)
RNA_CDS_STATUSES = {"SUPPORTED", "MIXED", "CONFLICT", "UNRESOLVED"}
RNA_CDS_SAMPLE_COLUMNS = (
    "rna_cds_sample_provenance_version", "rna_cds_support_sample_sources_json",
    "rna_cds_support_sample_count",
)


def full_cds_sample_metrics(sample_sources):
    """Count distinct canonical sample-GTF sources, not biological replicates."""
    sources = sorted(set(sample_sources))
    return {
        "rna_cds_sample_provenance_version": 1,
        "rna_cds_support_sample_sources_json": json.dumps(sources, separators=(",", ":")),
        "rna_cds_support_sample_count": len(sources),
    }


def terminal_positions(cds, strand):
    """Return codon coordinates in transcription order, including split codons."""
    ordered = sorted(cds, key=lambda part: part.start, reverse=strand == "-")
    if strand not in {"+", "-"} or not ordered:
        raise ValueError("CDS terminal coordinates require a known strand and CDS")

    def positions(parts, reverse=False):
        descending = (strand == "-") != reverse
        return chain.from_iterable(
            range(p.end, p.start - 1, -1) if descending
            else range(p.start, p.end + 1)
            for p in parts
        )

    start = tuple(islice(positions(ordered), ordered[0].phase, ordered[0].phase + 3))
    stop = tuple(reversed(tuple(islice(positions(reversed(ordered), True), 3))))
    return start, stop


def transcript_cds_structure(candidate, transcript):
    """Compare whole RNA introns with CDS, never clip introns at CDS boundaries."""
    blocks = candidate.intervals
    introns = transcript.introns
    start, stop = terminal_positions(candidate.cds, candidate.strand)
    same = candidate.strand == transcript.strand and candidate.seqid == transcript.seqid
    exonic = overlap_bp(blocks, transcript.exons) if same else 0
    conflicts = tuple(intron for intron in introns
                      if same and overlap_bp(blocks, (intron,)) > 0)
    missing = tuple(sorted(set(candidate.introns) - set(introns)))
    return {
        "cds_exonic_bp": exonic,
        "cds_bp": interval_bp(blocks),
        "cds_bases_inside_rna_introns": overlap_bp(blocks, conflicts),
        "start_codon_exonic_bases": sum(
            same and any(a <= p <= b for a, b in transcript.exons) for p in start
        ),
        "stop_codon_exonic_bases": sum(
            same and any(a <= p <= b for a, b in transcript.exons) for p in stop
        ),
        "candidate_introns_missing_from_rna": missing,
        "conflicting_introns": conflicts,
        "full_cds_compatible": bool(
            same and exonic == interval_bp(blocks) and not conflicts and not missing
            and len(start) == len(stop) == 3
        ),
    }


def _strong(evidence):
    return evidence is not None and (
        evidence.total_hq_reads >= 3 and evidence.sample_count_hq_ge1 >= 2
        and evidence.condition_count_ge2rep_hq_ge1 >= 1
    )


def rna_cds_metrics(candidate, transcripts, junctions):
    """Record all intronic CDS overlaps and the candidate's own RNA evidence.

    Whole-CDS enclosure is observable even without an exonic anchor. A complete
    compatible RNA counters any contradiction. A partial RNA carrying the whole
    coding intron chain can counter INTERNAL exon skipping only when every own
    junction has replicated strong support. It cannot validate coding termini.
    Status describes observations; unopposed conflicts drive the v5 penalty.
    Neither an assembled RNA nor separate junction reads prove a full molecule.
    """
    supporting, conflicting, chain_support = set(), set(), set()
    supporting_sources = set()
    recurrent = {}
    for transcript in transcripts:
        if (transcript.seqid, transcript.strand) != (candidate.seqid, candidate.strand):
            continue
        metrics = transcript_cds_structure(candidate, transcript)
        if metrics["full_cds_compatible"]:
            supporting.add(transcript.transcript_id)
            supporting_sources.update(transcript.sample_sources)
        if (metrics["cds_exonic_bp"] and not metrics["conflicting_introns"]
                and not metrics["candidate_introns_missing_from_rna"]):
            chain_support.add(transcript.transcript_id)
        for start, end in metrics["conflicting_introns"]:
            evidence = junctions.get((candidate.seqid, start, end, candidate.strand))
            if not _strong(evidence):
                continue
            conflicting.add(transcript.transcript_id)
            recurrent[start, end] = {
                "start": start, "end": end,
                "hq_reads": evidence.total_hq_reads,
                "samples": evidence.sample_count_hq_ge1,
                "replicated_conditions": evidence.condition_count_ge2rep_hq_ge1,
            }
    start, stop = terminal_positions(candidate.cds, candidate.strand)
    conflicts = sorted(recurrent)
    self_strong_chain = bool(candidate.introns) and bool(chain_support) and all(
        _strong(junctions.get((candidate.seqid, a, b, candidate.strand)))
        for a, b in candidate.introns
    )
    countered = [
        (a, b) for a, b in conflicts if supporting or (
            self_strong_chain and not any(a <= p <= b for p in (*start, *stop))
            and not overlap_bp((candidate.intervals[0], candidate.intervals[-1]), ((a, b),))
        )
    ]
    unopposed = [key for key in conflicts if key not in countered]
    status = (
        "MIXED" if supporting and conflicts else
        "SUPPORTED" if supporting else
        "CONFLICT" if conflicts else "UNRESOLVED"
    )
    encode = lambda value: json.dumps(value, separators=(",", ":"))
    return {
        **full_cds_sample_metrics(supporting_sources),
        "rna_cds_boundary_version": 2,
        "rna_cds_status": status,
        "rna_cds_support_transcripts_json": encode(sorted(supporting)),
        "rna_cds_conflict_transcripts_json": encode(sorted(conflicting)),
        "rna_cds_recurrent_conflicts_json": encode([recurrent[key] for key in conflicts]),
        "rna_cds_conflict_bp": overlap_bp(candidate.intervals, conflicts),
        "rna_cds_start_conflict": int(any(a <= p <= b for p in start for a, b in conflicts)),
        "rna_cds_stop_conflict": int(any(a <= p <= b for p in stop for a, b in conflicts)),
        "rna_cds_chain_support_transcripts_json": encode(sorted(chain_support)),
        "rna_cds_self_strong_chain": int(self_strong_chain),
        "rna_cds_countered_conflicts_json": encode([recurrent[key] for key in countered]),
        "rna_cds_unopposed_conflicts_json": encode([recurrent[key] for key in unopposed]),
        "rna_cds_unopposed_conflict_fraction": (
            overlap_bp(candidate.intervals, unopposed) / interval_bp(candidate.intervals)
        ),
    }
