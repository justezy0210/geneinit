"""Snakemake adapter for deterministic cross-sample junction aggregation."""

from pathlib import Path
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.evidence.junctions import (  # noqa: E402
    AGGREGATE_COLUMNS,
    aggregate_rows,
    write_bgzf_indexed,
    write_metric_table,
)
from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402


def fasta_reference_order(path: Path) -> list[str]:
    result: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                result.append(raw[1:].split(maxsplit=1)[0])
    if not result or len(result) != len(set(result)):
        raise ValueError("genome FASTA has no records or duplicate IDs")
    return result


log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
try:
    paths = [Path(value) for value in snakemake.input.tables]
    indices = [Path(value) for value in snakemake.input.indices]
    samples = list(snakemake.params.samples)
    conditions = list(snakemake.params.conditions)
    if not (len(paths) == len(indices) == len(samples) == len(conditions)):
        raise ValueError("junction table/index/sample/condition counts differ")
    for path in [*paths, *indices]:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing or empty junction input: {path}")
    references = fasta_reference_order(Path(snakemake.input.genome))
    rank = {seqid: index for index, seqid in enumerate(references)}
    output = Path(snakemake.output.table)
    write_bgzf_indexed(output, AGGREGATE_COLUMNS, aggregate_rows(paths, conditions, rank))
    write_metric_table(
        Path(snakemake.output.summary),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("sample_count", len(samples)),
            ("condition_count", len(set(conditions))),
            ("sample_order", ",".join(samples)),
            ("condition_order", ",".join(conditions)),
            ("min_mapq", snakemake.params.min_mapq),
            ("min_anchor", snakemake.params.min_anchor),
        ),
    )
    write_completion(
        Path(snakemake.output.complete),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("assembly_id", snakemake.params.assembly_id),
            ("aggregate_sha256", sha256_file(output)),
            ("aggregate_index_sha256", sha256_file(Path(snakemake.output.index))),
            ("summary_sha256", sha256_file(Path(snakemake.output.summary))),
        ),
    )
except Exception:
    log.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log.write_text("junction aggregation PASS\n", encoding="utf-8")
