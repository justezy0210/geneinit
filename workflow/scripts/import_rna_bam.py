"""Stage content-locked precomputed RNA alignments into this execution."""

from pathlib import Path
import sys
import traceback

sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.alignments import stage_alignment
from geneinit.contracts.samples import load_samples

log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
try:
    records = load_samples(Path(snakemake.input.samples), input_mode="bam")
    sample = next(row for row in records if row.sample == snakemake.wildcards.sample)
    if sample.layout != snakemake.wildcards.layout:
        raise ValueError("import layout differs from the locked sample metadata")
    stage_alignment(sample, Path(snakemake.input.manifest), Path(snakemake.output.bam),
                    Path(snakemake.output.csi), Path(snakemake.output.audit))
except Exception:
    log.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log.write_text(f"GeneInIt RNA BAM import PASS: {sample.sample}\n", encoding="utf-8")
