"""Validate the external licensed predictor image before expensive work."""

from pathlib import Path
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.deployment import (  # noqa: E402
    inspect_container,
    load_predictor_deployment,
    validate_predictor_versions,
)
from geneinit.predictor.common import write_completion  # noqa: E402


log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
try:
    deployment = load_predictor_deployment(dict(snakemake.config))
    versions = inspect_container(deployment)
    parsed_versions = validate_predictor_versions(versions)
    versions_path = Path(snakemake.output.versions)
    versions_path.parent.mkdir(parents=True, exist_ok=True)
    versions_path.write_text(versions, encoding="utf-8")
    write_completion(
        Path(snakemake.output.complete),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("deployment_mode", deployment.mode),
            ("predictor_sif_sha256", deployment.sif_sha256),
            ("genemark_license_acknowledged", "true"),
            ("augustus_version", parsed_versions["AUGUSTUS"]),
            ("genemark_etp_version", parsed_versions["GeneMark-ETP"]),
        ),
    )
except Exception:
    log.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log.write_text("predictor container preflight PASS\n", encoding="utf-8")
