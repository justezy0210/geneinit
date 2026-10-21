"""Regression checks for Python files executed through Snakemake ``script:``."""

from pathlib import Path


def test_workflow_scripts_compile_after_snakemake_preamble() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    scripts = sorted((repo_root / "workflow/scripts").glob("*.py"))
    assert scripts

    preamble = "snakemake = object()\n"
    for script in scripts:
        compile(preamble + script.read_text(encoding="utf-8"), str(script), "exec")
