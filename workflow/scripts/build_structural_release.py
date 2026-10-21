from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.release.core import build_structural_release  # noqa: E402


build_structural_release(
    Path(snakemake.input.gff3),
    Path(snakemake.input.transcript_decisions),
    Path(snakemake.input.edge_audit),
    Path(snakemake.input.bridge_reviews),
    Path(snakemake.input.genome),
    Path(snakemake.input.policy),
    [str(value) for value in snakemake.params.seqids],
    str(snakemake.params.assembly_id),
    str(snakemake.params.gene_prefix),
    Path(snakemake.output.gff3),
    Path(snakemake.output.transcripts_fna),
    Path(snakemake.output.cds_fna),
    Path(snakemake.output.proteins_faa),
    Path(snakemake.output.longest_proteins_faa),
    Path(snakemake.output.gene_id_map),
    Path(snakemake.output.transcript_id_map),
    Path(snakemake.output.candidate_lineage),
    Path(snakemake.output.structural_decisions),
    Path(snakemake.output.review_flags),
    Path(snakemake.output.qc_tsv),
    Path(snakemake.output.qc_json),
    Path(snakemake.output.release_status),
    Path(snakemake.output.manifest),
    Path(snakemake.output.checksums),
    Path(snakemake.output.policy),
    Path(snakemake.output.complete),
    Path(snakemake.input.te_evidence[0]) if snakemake.input.te_evidence else None,
    Path(snakemake.input.candidates),
    Path(snakemake.input.selection_audit),
    Path(snakemake.input.generated_te),
)
