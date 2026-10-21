import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from geneinit import __version__
from geneinit.cli import parse_arguments


@pytest.mark.parametrize("arguments", [[], ["--threads", "0"], ["--threads", "-2"],
                                      ["--threads", "all"], ["--threads", "2", "--cores", "3"]])
def test_invalid_or_ambiguous_budget_fails(arguments):
    with pytest.raises(SystemExit) as error:
        parse_arguments(arguments)
    assert error.value.code == 2


@pytest.mark.parametrize("option", ["--threads", "--cores", "-c"])
@pytest.mark.parametrize("budget", [1, 2, 34])
def test_launcher_passes_budget_and_affinity_to_child(tmp_path, option, budget):
    executable = tmp_path / "snakemake"
    executable.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "print(json.dumps({'args':sys.argv[1:], 'affinity':sorted(os.sched_getaffinity(0)),"
        "'blas':os.environ['OPENBLAS_NUM_THREADS']}))\n"
    )
    executable.chmod(0o755)
    environment = dict(os.environ, PATH=str(tmp_path) + os.pathsep + os.environ["PATH"],
                       OPENBLAS_NUM_THREADS="99")
    result = subprocess.run(
        [sys.executable, "-m", "geneinit.cli", option, str(budget),
         "--configfile", "path with spaces.yaml", "--dry-run"],
        env=environment, capture_output=True, text=True, check=True,
    )
    child = json.loads(result.stdout)
    effective = min(budget, len(os.sched_getaffinity(0)))
    assert len(child["affinity"]) == effective
    assert set(child["affinity"]) <= os.sched_getaffinity(0)
    assert child["args"][child["args"].index("--cores") + 1] == str(effective)
    assert child["args"][-3:] == ["--configfile", "path with spaces.yaml", "--dry-run"]
    assert child["blas"] == "1"


def test_alias_equals_and_help():
    assert parse_arguments(["--threads=34", "--cores=34", "-n"]) == (34, ["-n"])
    assert parse_arguments(["-c2", "all"]) == (2, ["all"])
    assert parse_arguments(["--version"]) == (None, ["--version"])
    with pytest.raises(SystemExit) as error:
        parse_arguments(["--help"])
    assert error.value.code == 0


@pytest.mark.parametrize("budget", [[], ["--threads", "2"]])
def test_version_reports_geneinit_without_launching_snakemake(tmp_path, budget):
    executable = tmp_path / "snakemake"
    executable.write_text(f"#!{sys.executable}\nraise SystemExit(93)\n")
    executable.chmod(0o755)
    environment = dict(os.environ, PATH=str(tmp_path))
    result = subprocess.run(
        [sys.executable, "-m", "geneinit.cli", *budget, "--version"],
        env=environment, capture_output=True, text=True, check=True,
    )
    assert result.stdout.strip() == f"GeneInIt {__version__}"
    assert result.stderr == ""
