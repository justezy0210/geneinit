from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.utr.grouping import group_coding_models  # noqa: E402


group_coding_models(
    Path(snakemake.input.gff3),
    Path(snakemake.input.mapping),
    Path(snakemake.input.candidates),
    Path(snakemake.input.source_models),
    [str(value) for value in snakemake.params.seqids],
    Path(snakemake.output.gff3),
    Path(snakemake.output.mapping),
    Path(snakemake.output.pairs),
    Path(snakemake.output.summary),
    Path(snakemake.output.complete),
)
