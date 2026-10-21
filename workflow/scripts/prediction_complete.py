"""Create the predictor-stage aggregate completion contract."""

from pathlib import Path
import sys


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.predictor.common import require_nonempty, sha256_file, write_completion  # noqa: E402


for value in snakemake.input:
    require_nonempty(Path(value), "prediction-stage input")
write_completion(
    Path(snakemake.output.complete),
    (
        ("schema_version", 1),
        ("status", "PASS"),
        ("assembly_id", snakemake.params.assembly_id),
        ("genemark_complete_sha256", sha256_file(Path(snakemake.input.genemark))),
        ("hints_complete_sha256", sha256_file(Path(snakemake.input.hints))),
        ("augustus_complete_sha256", sha256_file(Path(snakemake.input.augustus))),
        *(("masking_rescue_complete_sha256", sha256_file(Path(path)))
          for path in getattr(snakemake.input, "masking_rescue", [])),
    ),
)
log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text("prediction stage PASS\n", encoding="utf-8")
