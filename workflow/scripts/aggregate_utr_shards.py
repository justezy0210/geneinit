from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.utr.core import aggregate_utr_shards  # noqa: E402


aggregate_utr_shards(
    [Path(value) for value in snakemake.input.gff3s],
    [Path(value) for value in snakemake.input.mappings],
    Path(snakemake.input.preparation),
    Path(snakemake.output.gff3),
    Path(snakemake.output.mapping),
    Path(snakemake.output.summary),
    Path(snakemake.output.complete),
)
