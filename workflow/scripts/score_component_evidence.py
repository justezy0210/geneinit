from collections import Counter
import json
from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.evidence.components import (  # noqa: E402
    EvidenceContext,
    iter_candidate_blocks,
    load_junction_aggregate,
    write_evidence_table,
)
from geneinit.evidence.repeats import REPEAT_SCOPE  # noqa: E402
from geneinit.evidence.junctions import write_metric_table  # noqa: E402
from geneinit.models.candidate import read_fasta  # noqa: E402
from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402


sample_gtfs = list(getattr(snakemake.input, "sample_gtfs", ()))
junctions = load_junction_aggregate(Path(snakemake.input.junctions))
genome = read_fasta(Path(snakemake.input.genome))
context = EvidenceContext.from_paths(
    {**{name: getattr(snakemake.input, name) for name in (
        "merged_gtf", "nonhc_trace", "direct_chains", "nonhc_chains",
        "direct_hits", "nonhc_hits", "repeats",
    )}, "sample_gtfs": sample_gtfs},
    genome=genome, junctions=junctions,
    expected_protein_copies=int(snakemake.params.expected_protein_copies),
)

counts = Counter()
orf_counts = Counter()
multiplicity_counts = Counter()
rna_boundary_counts = Counter()


def counted_rows():
    for row in context.score_blocks(
        iter_candidate_blocks(Path(snakemake.input.candidates)),
    ):
        counts["candidate"] += 1
        orf_counts[str(row["orf_status"])] += 1
        multiplicity_counts[str(row["protein_multiplicity_status"])] += 1
        rna_boundary_counts[str(row["rna_cds_status"])] += 1
        counts["predictor_concordant"] += int(row["predictor_tier"]) == 3
        counts["rna_supported"] += int(row["rna_tier"]) >= 2
        counts["protein_supported"] += int(row["protein_tier"]) >= 2
        counts["high_repeat"] += float(row["repeat_cds_fraction"]) >= 0.5
        counts["severe_repeat"] += float(row["repeat_cds_fraction"]) >= 0.8
        yield row


write_evidence_table(Path(snakemake.output.table), counted_rows())

write_metric_table(
    Path(snakemake.output.summary),
    (
        ("assembly_id", snakemake.params.assembly_id),
        ("candidate_count", counts["candidate"]),
        ("merged_rna_transcript_count", context.counts["merged_rna_transcript_count"]),
        ("boundary_rna_distinct_chain_count", context.counts["boundary_rna_distinct_chain_count"]),
        ("boundary_rna_sample_gtf_count", len(sample_gtfs)),
        ("protein_chain_count", context.counts["protein_chain_count"]),
        ("repeat_feature_count", context.counts["repeat_feature_count"]),
        ("classified_te_feature_count", context.counts["classified_te_feature_count"]),
        ("repeat_scope", REPEAT_SCOPE),
        ("junction_aggregate_count", context.counts["junction_aggregate_count"]),
        ("complete_orf_count", orf_counts.get("complete", 0)),
        ("predictor_concordant_count", counts["predictor_concordant"]),
        ("rna_supported_candidate_count", counts["rna_supported"]),
        ("protein_supported_candidate_count", counts["protein_supported"]),
        ("unique_protein_hit_count", multiplicity_counts.get("UNIQUE_PROTEIN_HIT", 0)),
        ("ambiguous_retained_protein_hit_count", multiplicity_counts.get("AMBIGUOUS_PROTEIN_MULTI_HIT", 0)),
        ("high_repeat_candidate_count", counts["high_repeat"]),
        ("severe_repeat_candidate_count", counts["severe_repeat"]),
        *((f"rna_cds_{status.lower()}_count", rna_boundary_counts[status])
          for status in ("SUPPORTED", "MIXED", "CONFLICT", "UNRESOLVED")),
    ),
)

write_completion(
    Path(snakemake.output.complete),
    (
        ("schema_version", 4),
        ("rna_cds_boundary_version", 2),
        ("rna_cds_sample_provenance_version", 1),
        ("sample_gtf_sha256_json", json.dumps(
            {str(p): sha256_file(Path(p)) for p in sample_gtfs}, sort_keys=True)),
        ("repeat_scope", REPEAT_SCOPE),
        ("status", "PASS"),
        ("assembly_id", snakemake.params.assembly_id),
        ("candidate_count", counts["candidate"]),
        ("expected_protein_copies", snakemake.params.expected_protein_copies),
        ("component_evidence_sha256", sha256_file(Path(snakemake.output.table))),
        ("component_summary_sha256", sha256_file(Path(snakemake.output.summary))),
    ),
)
