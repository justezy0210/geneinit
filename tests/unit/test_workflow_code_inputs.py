"""Imported helpers must invalidate their consumers, without global reruns."""

import ast
from importlib.util import resolve_name
import os
from pathlib import Path
import subprocess
import sys

import pytest

from geneinit.contracts.workflow import SOURCE_GROUPS, source_inputs


SOURCE = Path(__file__).resolve().parents[2]
SCRIPT_GROUPS = {
    "preflight.py": "preflight",
    "import_rna_bam.py": "rna_import",
    "predictor_container_preflight.py": "predictor_container",
    "run_genemark_etp.py": "genemark",
    "build_etp_hints.py": "prediction_hints",
    "prepare_augustus_config.py": "augustus_config",
    "run_augustus_etp.py": "augustus",
    "prediction_complete.py": "prediction_complete",
    "run_masking_rescue.py": "masking_rescue",
    "extract_junctions.py": "junction_extract",
    "aggregate_junctions.py": "junction_aggregate",
    "build_candidate_universe.py": "models",
    "score_component_evidence.py": "evidence",
    "select_coding_models.py": "selection",
    "prepare_selected_coding_models.py": "coding_prepare",
    "group_coding_models.py": "coding_group",
    "run_structural_postqc.py": "postqc",
    "prepare_postqc_utr_shards.py": "utr_prepare",
    "refine_utr_shard.py": "utr_refine",
    "aggregate_utr_shards.py": "utr_aggregate",
    "build_structural_release.py": "release",
}


def _imported_local_files(path):
    """Test-only check of explicit groups; never used to build the workflow."""
    relative = path.relative_to(SOURCE)
    if relative.parts[0] == "geneinit":
        package = ".".join(relative.with_suffix("").parts[:-1])
    else:
        package = ""
    for node in ast.walk(ast.parse(path.read_text())):
        names = []
        if isinstance(node, ast.ImportFrom):
            name = node.module or ""
            if node.level:
                name = resolve_name("." * node.level + name, package)
            names.append(name)
        elif isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
        for name in names:
            if name != "geneinit" and not name.startswith("geneinit."):
                continue
            module = SOURCE.joinpath(*name.split("."))
            imported = module.with_suffix(".py")
            if not imported.is_file():
                imported = module / "__init__.py"
            assert imported.is_file(), name
            yield imported
            for directory in imported.parents:
                if directory == SOURCE:
                    break
                yield directory / "__init__.py"


@pytest.mark.parametrize("script,group", SCRIPT_GROUPS.items())
def test_explicit_group_covers_direct_and_transitive_local_imports(script, group):
    pending = [SOURCE / "workflow/scripts" / script]
    visited = set(pending)
    imported = set()
    while pending:
        for path in _imported_local_files(pending.pop()):
            imported.add(path)
            if path not in visited:
                visited.add(path)
                pending.append(path)
    declared = set(map(Path, source_inputs(SOURCE, group)))
    assert imported <= declared, sorted(str(path.relative_to(SOURCE)) for path in imported - declared)


def test_source_inputs_are_existing_sorted_unique_and_narrow():
    for group in SOURCE_GROUPS:
        paths = source_inputs(SOURCE, group, group)
        assert paths == sorted(set(paths))
        assert all(Path(path).is_absolute() and Path(path).is_file() for path in paths)
    selection = source_inputs(SOURCE, "selection")
    assert str(SOURCE / "geneinit/selection/unified.py") in selection
    assert str(SOURCE / "geneinit/predictor/training.py") not in selection
    assert str(SOURCE / "geneinit/contracts/preflight.py") not in selection
    assert str(SOURCE / "geneinit/selection/unified.py") not in source_inputs(SOURCE, "evidence")
    assert str(SOURCE / "geneinit/release/core.py") not in source_inputs(SOURCE, "postqc")
    with pytest.raises(KeyError):
        source_inputs(SOURCE, "misspelled_stage")


def test_source_inputs_fail_for_missing_code(tmp_path):
    with pytest.raises(FileNotFoundError, match="workflow source input is missing"):
        source_inputs(tmp_path, "selection")


def test_evidence_import_does_not_load_preflight_and_public_contract_api_is_preserved():
    code = """
import sys
import geneinit.evidence.components
assert 'geneinit.contracts.preflight' not in sys.modules
assert 'geneinit.contracts.deployment' not in sys.modules
import geneinit.contracts as contracts
from geneinit.contracts import SampleRecord, SampleContractError, load_samples
from geneinit.contracts.samples import SampleRecord as original_record
assert SampleRecord is original_record
assert 'geneinit.contracts.preflight' not in sys.modules
from geneinit.contracts import PreflightError, run_preflight
from geneinit.contracts.preflight import run_preflight as original_preflight
assert run_preflight is original_preflight
assert set(contracts.__all__) <= set(dir(contracts))
try:
    contracts.not_a_public_contract
except AttributeError as error:
    assert "has no attribute 'not_a_public_contract'" in str(error)
else:
    raise AssertionError('missing contract attribute did not fail')
"""
    run = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True,
        env=dict(os.environ, PYTHONPATH=str(SOURCE), PYTHONDONTWRITEBYTECODE="1"),
    )
    assert run.returncode == 0, run.stdout + run.stderr
