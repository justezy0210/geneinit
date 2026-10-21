from pathlib import Path
import os
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTAINERS = REPO_ROOT / "containers"


def _read(relative: str) -> str:
    return (CONTAINERS / relative).read_text(encoding="utf-8")


def test_reusable_container_files_are_present() -> None:
    required = {
        "README.md",
        "predictor/Dockerfile",
        "predictor/build_sif.sh",
        "allinone/Dockerfile",
        "allinone/build_sif.sh",
        "allinone/environment-main.yaml",
        "allinone/environment-alignment.yaml",
        "allinone/environment-stringtie.yaml",
        "allinone/environment-utr.yaml",
        "allinone/geneinit-run",
        "allinone/geneinit-software-versions",
    }
    assert all((CONTAINERS / path).is_file() for path in required)


def test_reusable_recipes_have_no_project_or_user_paths() -> None:
    recipe_files = [
        path
        for path in CONTAINERS.rglob("*")
        if path.is_file() and path.name != "manifest.tsv"
    ]
    forbidden = ("/10Gdata/", "/home/ezy", "jinju-station")
    for path in recipe_files:
        text = path.read_text(encoding="utf-8")
        assert not any(token in text for token in forbidden), path


def test_images_never_embed_genemark_license_key() -> None:
    for relative in ("predictor/Dockerfile", "allinone/Dockerfile"):
        assert 'org.geneinit.license-key.embedded="false"' in _read(relative)


def test_every_etp_perl_entry_point_is_pinned_and_checked() -> None:
    for relative in ("predictor/Dockerfile", "allinone/Dockerfile"):
        text = _read(relative)
        assert "find /opt/ETP/bin -type f -name '*.pl'" in text
        assert "#!/usr/bin/perl" in text
        assert "^#!/usr/bin/env perl$" in text
        assert "/usr/bin/perl -c" in text


def test_alignment_and_stringtie_dependencies_stay_separate() -> None:
    alignment = _read("allinone/environment-alignment.yaml")
    stringtie = _read("allinone/environment-stringtie.yaml")
    main = _read("allinone/environment-main.yaml")
    assert "samtools=1.24" in alignment
    assert "stringtie" not in alignment
    assert "stringtie=3.0.3" in stringtie
    assert "samtools" not in stringtie
    assert "samtools" not in main
    assert "stringtie" not in main


def test_builds_use_temporary_home_for_license_smoke() -> None:
    for relative in ("predictor/build_sif.sh", "allinone/build_sif.sh"):
        text = _read(relative)
        assert 'install -m 600 "$genemark_key" "$smoke_home/.gm_key"' in text
        assert '--home "$smoke_home"' in text


@pytest.mark.parametrize("tool", ["hisat2", "samtools"])
@pytest.mark.parametrize("exit_code", [0, 7])
def test_version_summary_drains_output_and_preserves_failure(tmp_path, tool, exit_code) -> None:
    command = next(
        line for line in _read("allinone/geneinit-software-versions").splitlines()
        if line.startswith(f"{tool} --version")
    )
    executable = tmp_path / tool
    executable.write_text(
        f"#!{sys.executable}\n"
        "import signal, sys\n"
        "signal.signal(signal.SIGPIPE, signal.SIG_DFL)\n"
        "sys.stdout.write('tool version 1.0\\n' + 'details\\n' * 131072)\n"
        "sys.stdout.flush()\n"
        f"sys.exit({exit_code})\n"
    )
    executable.chmod(0o755)
    result = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", command],
        env={**os.environ, "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"]},
        capture_output=True, text=True, timeout=10,
    )
    assert result.returncode == exit_code, result.stderr
    assert result.stdout == "tool version 1.0\n"
