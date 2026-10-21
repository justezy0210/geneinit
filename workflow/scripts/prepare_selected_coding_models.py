from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402
from geneinit.utr.core import prepare_selected_coding_models  # noqa: E402


gff3 = Path(snakemake.output.gff3)
mapping = Path(snakemake.output.mapping)
summary = Path(snakemake.output.summary)
rows = prepare_selected_coding_models(
    Path(snakemake.input.candidates),
    Path(snakemake.input.selection_audit),
    [str(value) for value in snakemake.params.seqids],
    gff3,
    mapping,
)
summary.parent.mkdir(parents=True, exist_ok=True)
summary.write_text(
    "metric\tvalue\n"
    "status\tPASS_MIGRATION_PENDING\n"
    f"selected_coding_model_count\t{len(rows)}\n"
    f"coding_gff3_sha256\t{sha256_file(gff3)}\n"
    f"coding_mapping_sha256\t{sha256_file(mapping)}\n",
    encoding="utf-8",
)
write_completion(
    Path(snakemake.output.complete),
    (
        ("schema_version", 1),
        ("status", "PASS_MIGRATION_PENDING"),
        ("selected_coding_model_count", len(rows)),
        ("coding_gff3_sha256", sha256_file(gff3)),
        ("coding_mapping_sha256", sha256_file(mapping)),
        ("summary_sha256", sha256_file(summary)),
    ),
)
