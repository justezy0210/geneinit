from collections import Counter
from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.evidence.components import (  # noqa: E402
    iter_candidate_blocks,
    load_diamond_targets,
    load_junction_aggregate,
    load_protein_chains,
    load_repeats,
    load_stringtie_transcripts,
    load_trace,
    retained_protein_targets,
    score_candidate_blocks,
    write_evidence_table,
)
from geneinit.evidence.repeats import REPEAT_SCOPE, is_classified_te  # noqa: E402
from geneinit.evidence.junctions import write_metric_table  # noqa: E402
from geneinit.models.candidate import read_fasta  # noqa: E402
from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402


transcripts = load_stringtie_transcripts(Path(snakemake.input.merged_gtf))
trace = load_trace(Path(snakemake.input.nonhc_trace))
chains = load_protein_chains(
    Path(snakemake.input.direct_chains), "direct"
) + load_protein_chains(
    Path(snakemake.input.nonhc_chains), "nonhc", trace=trace
)
seed_targets = load_diamond_targets(
    Path(snakemake.input.direct_hits), Path(snakemake.input.nonhc_hits)
)
retained_targets = retained_protein_targets(chains, seed_targets)
repeats = load_repeats(Path(snakemake.input.repeats))
junctions = load_junction_aggregate(Path(snakemake.input.junctions))
genome = read_fasta(Path(snakemake.input.genome))

counts = Counter()
orf_counts = Counter()
multiplicity_counts = Counter()


def counted_rows():
    for row in score_candidate_blocks(
        iter_candidate_blocks(Path(snakemake.input.candidates)),
        transcripts,
        chains,
        repeats,
        junctions,
        genome,
        retained_targets,
        seed_targets,
        int(snakemake.params.expected_protein_copies),
    ):
        counts["candidate"] += 1
        orf_counts[str(row["orf_status"])] += 1
        multiplicity_counts[str(row["protein_multiplicity_status"])] += 1
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
        ("merged_rna_transcript_count", len(transcripts)),
        ("protein_chain_count", len(chains)),
        ("repeat_feature_count", len(repeats)),
        ("classified_te_feature_count", sum(is_classified_te(r.classification) for r in repeats)),
        ("repeat_scope", REPEAT_SCOPE),
        ("junction_aggregate_count", len(junctions)),
        ("complete_orf_count", orf_counts.get("complete", 0)),
        ("predictor_concordant_count", counts["predictor_concordant"]),
        ("rna_supported_candidate_count", counts["rna_supported"]),
        ("protein_supported_candidate_count", counts["protein_supported"]),
        ("unique_protein_hit_count", multiplicity_counts.get("UNIQUE_PROTEIN_HIT", 0)),
        ("ambiguous_retained_protein_hit_count", multiplicity_counts.get("AMBIGUOUS_PROTEIN_MULTI_HIT", 0)),
        ("high_repeat_candidate_count", counts["high_repeat"]),
        ("severe_repeat_candidate_count", counts["severe_repeat"]),
    ),
)

write_completion(
    Path(snakemake.output.complete),
    (
        ("schema_version", 2),
        ("repeat_scope", REPEAT_SCOPE),
        ("status", "PASS"),
        ("assembly_id", snakemake.params.assembly_id),
        ("candidate_count", counts["candidate"]),
        ("expected_protein_copies", snakemake.params.expected_protein_copies),
        ("component_evidence_sha256", sha256_file(Path(snakemake.output.table))),
        ("component_summary_sha256", sha256_file(Path(snakemake.output.summary))),
    ),
)
