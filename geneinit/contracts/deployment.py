"""Fail-closed deployment checks for the direct predictor container."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


REQUIRED_CONTAINER_PATHS = (
    "/usr/bin/perl",
    "/opt/ETP/bin/gmetp.pl",
    "/opt/ETP/tools/diamond",
    "/opt/ETP/bin/format_back.pl",
    "/opt/ETP/bin/GeneMarkSTFiltering/proteins_from_gtf.pl",
    "/opt/Augustus/bin/augustus",
    "/opt/Augustus/bin/etraining",
    "/opt/Augustus/scripts/new_species.pl",
    "/opt/Augustus/scripts/gff2gbSmallDNA.pl",
    "/opt/Augustus/scripts/optimize_augustus.pl",
    "/opt/Augustus/scripts/join_mult_hints.pl",
    "/opt/Augustus/scripts/join_aug_pred.pl",
    "/opt/Augustus/config/species/generic/generic_parameters.cfg",
)
REQUIRED_PREDICTOR_VERSIONS = {
    "AUGUSTUS": "3.5.0",
    "GeneMark-ETP": "1.02",
}


class DeploymentError(ValueError):
    """Raised when a licensed predictor deployment is missing or inconsistent."""


@dataclass(frozen=True)
class PredictorDeployment:
    mode: str
    sif: Path
    sif_sha256: str
    apptainer: Path | None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_predictor_deployment(config: Mapping[str, Any]) -> PredictorDeployment:
    raw = config.get("deployment")
    if not isinstance(raw, Mapping):
        raise DeploymentError(
            "prediction requires a deployment mapping; see config/config.example.yaml"
        )
    if raw.get("genemark_license_acknowledged") is not True:
        raise DeploymentError(
            "prediction requires deployment.genemark_license_acknowledged: true"
        )

    mode = str(raw.get("mode", "apptainer"))
    if mode not in {"apptainer", "in_container"}:
        raise DeploymentError(
            "deployment.mode must be either 'apptainer' or 'in_container'"
        )

    sif_raw = raw.get("predictor_sif")
    if not isinstance(sif_raw, str) or not sif_raw.strip():
        raise DeploymentError("deployment.predictor_sif must be a non-empty path")
    sif = Path(sif_raw).expanduser().resolve()
    if not sif.is_file() or not os.access(sif, os.R_OK) or sif.stat().st_size == 0:
        raise DeploymentError(f"predictor SIF is absent, unreadable, or empty: {sif}")

    expected = str(raw.get("predictor_sif_sha256", "")).lower()
    if len(expected) != 64 or any(base not in "0123456789abcdef" for base in expected):
        raise DeploymentError("deployment.predictor_sif_sha256 must contain 64 hex digits")
    observed = sha256_file(sif)
    if observed != expected:
        raise DeploymentError(
            f"predictor SIF SHA-256 mismatch: expected={expected} observed={observed}"
        )

    executable = None
    if mode == "apptainer":
        executable_raw = str(raw.get("apptainer_executable", "apptainer"))
        executable = shutil.which(executable_raw)
        if executable is None:
            raise DeploymentError(f"Apptainer executable was not found: {executable_raw}")
    return PredictorDeployment(
        mode=mode,
        sif=sif,
        sif_sha256=observed,
        apptainer=Path(executable).resolve() if executable is not None else None,
    )


def bind_arguments(paths: Sequence[Path]) -> list[str]:
    """Return deterministic same-path binds without broadening file inputs."""

    roots = sorted(
        {
            str(path.resolve() if path.exists() else path.resolve().parent)
            for path in paths
        }
    )
    selected: list[str] = []
    for root_path in roots:
        if any(
            root_path == root or root_path.startswith(root.rstrip("/") + "/")
            for root in selected
        ):
            continue
        selected = [
            root
            for root in selected
            if not (root == root_path or root.startswith(root_path.rstrip("/") + "/"))
        ]
        selected.append(root_path)
    arguments: list[str] = []
    for parent in sorted(selected):
        arguments.extend(("--bind", f"{parent}:{parent}"))
    return arguments


def apptainer_command(
    deployment: PredictorDeployment,
    command: Sequence[str],
    bind_paths: Sequence[Path],
    pwd: Path | None = None,
) -> list[str]:
    # A child may already create a process/thread pool using --cores/--threads.
    # Set this inside the container too, because --cleanenv drops host settings.
    worker_environment = [
        f"{name}=1" for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
        )
    ]
    command = ["/usr/bin/env", *worker_environment, *map(str, command)]
    if deployment.mode == "in_container":
        result = [str(value) for value in command]
        if pwd is not None:
            result = ["/usr/bin/env", f"--chdir={pwd.resolve()}", *result]
        return result

    if deployment.apptainer is None:
        raise DeploymentError("Apptainer deployment has no executable")
    result = [str(deployment.apptainer), "exec", "--cleanenv"]
    result.extend(bind_arguments(bind_paths))
    if pwd is not None:
        result.extend(("--pwd", str(pwd.resolve())))
    result.append(str(deployment.sif))
    result.extend(str(value) for value in command)
    return result


def inspect_container(deployment: PredictorDeployment) -> str:
    checks = " && ".join(f"test -e {path}" for path in REQUIRED_CONTAINER_PATHS)
    perl_checks = (
        "/usr/bin/perl -MHash::Merge -MMCE::Mutex -MMath::Utils "
        "-MParallel::ForkManager -MYAML -MYAML::XS -e 1 && "
        "while IFS= read -r -d '' script; do "
        "/usr/bin/perl -c \"$script\" >/dev/null 2>&1; "
        "done < <(find /opt/ETP/bin -type f -name '*.pl' -print0)"
    )
    version_commands = (
        "/opt/Augustus/bin/augustus --version; "
        "grep -m1 'my .version' /opt/ETP/bin/gmetp.pl"
    )
    command = apptainer_command(
        deployment,
        ("/bin/bash", "-lc", f"set -e; {checks}; {perl_checks}; {version_commands}"),
        (),
    )
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise DeploymentError(
            f"predictor container inspection failed ({completed.returncode}): {detail}"
        )
    output = (completed.stdout + completed.stderr).strip()
    if not output:
        raise DeploymentError("predictor container returned no version information")
    return output + "\n"


def validate_predictor_versions(version_text: str) -> dict[str, str]:
    """Require the exact predictor stack used by the frozen command contracts."""

    patterns = {
        "AUGUSTUS": r"AUGUSTUS[^\n]*?([0-9]+(?:\.[0-9]+){2})",
        "GeneMark-ETP": r"(?:my\s+)?\$?version\s*=\s*[\"']([0-9]+(?:\.[0-9]+)+)",
    }
    observed: dict[str, str] = {}
    for tool, pattern in patterns.items():
        match = re.search(pattern, version_text, flags=re.IGNORECASE)
        if match is None:
            raise DeploymentError(f"could not parse {tool} version from container")
        observed[tool] = match.group(1)
        expected = REQUIRED_PREDICTOR_VERSIONS[tool]
        if observed[tool] != expected:
            raise DeploymentError(
                f"unsupported {tool} version: expected={expected} observed={observed[tool]}"
            )
    return observed
