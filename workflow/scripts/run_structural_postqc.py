from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.postqc.core import run_postqc  # noqa: E402


run_postqc(
    Path(snakemake.input.gff3),
    Path(snakemake.input.grouping_mapping),
    Path(snakemake.input.candidates),
    Path(snakemake.input.component_evidence),
    Path(snakemake.input.selection_audit),
    Path(snakemake.input.junctions),
    Path(snakemake.input.genome),
    Path(snakemake.input.policy),
    [str(value) for value in snakemake.params.seqids],
    Path(snakemake.output.gff3),
    Path(snakemake.output.transcripts),
    Path(snakemake.output.genes),
    Path(snakemake.output.edges),
    Path(snakemake.output.promotions),
    Path(snakemake.output.bridge_reviews),
    Path(snakemake.output.retirements),
    Path(snakemake.output.summary),
    Path(snakemake.output.complete),
    Path(snakemake.output.iterations),
    Path(snakemake.input.te_evidence[0]) if snakemake.input.te_evidence else None,
    str(snakemake.params.assembly_id),
    Path(snakemake.output.te_classification),
)
