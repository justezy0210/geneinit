"""Snakemake adapter for one-sample splice-junction counting."""

from pathlib import Path
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.evidence.junctions import (  # noqa: E402
    JUNCTION_COLUMNS,
    count_bam_junctions,
    junction_rows,
    write_bgzf_indexed,
    write_metric_table,
)


log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
try:
    bam = Path(snakemake.input.bam)
    csi = Path(snakemake.input.csi)
    if not bam.is_file() or bam.stat().st_size == 0:
        raise ValueError(f"missing or empty BAM: {bam}")
    if not csi.is_file() or csi.stat().st_size == 0:
        raise ValueError(f"missing or empty BAM CSI: {csi}")
    counts, stats, references = count_bam_junctions(
        bam,
        protocol=str(snakemake.params.strandedness),
        min_mapq=int(snakemake.params.min_mapq),
        min_anchor=int(snakemake.params.min_anchor),
        io_threads=max(1, int(snakemake.threads)),
    )
    write_bgzf_indexed(
        Path(snakemake.output.table), JUNCTION_COLUMNS, junction_rows(counts, references)
    )
    write_metric_table(
        Path(snakemake.output.summary),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("sample", snakemake.params.sample),
            ("condition", snakemake.params.condition),
            ("layout", snakemake.params.layout),
            ("strandedness", snakemake.params.strandedness),
            ("min_mapq", snakemake.params.min_mapq),
            ("min_anchor", snakemake.params.min_anchor),
            *sorted(stats.items()),
        ),
    )
except Exception:
    log.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log.write_text("junction extraction PASS\n", encoding="utf-8")
