"""Create a run-local writable AUGUSTUS configuration from the locked image."""

from pathlib import Path
import os
import stat
import subprocess
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.deployment import (  # noqa: E402
    apptainer_command,
    load_predictor_deployment,
)
from geneinit.predictor.common import PredictorError, write_completion  # noqa: E402


log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
try:
    deployment = load_predictor_deployment(dict(snakemake.config))
    target = Path(snakemake.output.config).resolve()
    target.mkdir(parents=True, exist_ok=True)
    command = apptainer_command(
        deployment,
        ("/bin/cp", "-a", "/opt/Augustus/config/.", f"{target}/"),
        (target,),
    )
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise PredictorError(
            f"copying AUGUSTUS config failed: {(completed.stderr or completed.stdout).strip()}"
        )
    required = target / "species/generic/generic_parameters.cfg"
    if not required.is_file() or required.stat().st_size == 0:
        raise PredictorError(f"incomplete AUGUSTUS config copy: {required}")
    for root, directories, files in os.walk(target):
        for name in directories:
            path = Path(root) / name
            path.chmod(path.stat().st_mode | stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
        for name in files:
            path = Path(root) / name
            path.chmod(path.stat().st_mode | stat.S_IRUSR | stat.S_IWUSR)
    write_completion(
        Path(snakemake.output.complete),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("container_sha256", deployment.sif_sha256),
            ("augustus_config", target),
        ),
    )
except Exception:
    log.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log.write_text("writable AUGUSTUS config PASS\n", encoding="utf-8")
