from pathlib import Path
import json
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.utr.grouping import group_coding_models  # noqa: E402


policy = json.loads(Path(snakemake.input.policy).read_text())

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
    defer_grouping=policy["selector_core"].get("engine") == "unified_v1",
)
