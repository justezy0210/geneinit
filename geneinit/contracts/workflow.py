"""Explicit local code inputs for Snakemake's incremental reruns.

Scripts are tracked by Snakemake itself. Imported project modules must also be
inputs: editing a selector must invalidate selection, without remapping RNA.
Keep these small dependency groups explicit when adding local imports.
"""

from pathlib import Path


_SAMPLES = ("contracts/samples.py",)
_ALIGNMENTS = _SAMPLES + ("contracts/alignments.py",)
_PREFLIGHT = _ALIGNMENTS + ("contracts/preflight.py", "contracts/deployment.py")
_COMMON = ("predictor/common.py",)
_MODEL = ("models/candidate.py",)
_STREAMING = ("contracts/streaming.py",)
_METRICS = ("evidence/metrics.py",)
_RNA = _METRICS + ("evidence/rna_cds.py",)
_REPEATS = ("evidence/repeats.py",)
_TE = ("evidence/te_classification.py",)
_JUNCTIONS = ("evidence/junctions.py",)
_COMPONENTS = _MODEL + _STREAMING + _RNA + _REPEATS + _JUNCTIONS + (
    "evidence/components.py",
)
_PREDICTION = _COMMON + (
    "predictor/commands.py", "predictor/parallel.py", "predictor/resume.py",
)
_MASKING = _MODEL + _METRICS + _REPEATS + _PREDICTION + ("predictor/masking.py",)
_SELECTION = _MODEL + _STREAMING + _RNA + _REPEATS + _TE + (
    "selection/core.py", "selection/balanced.py", "selection/unified.py",
)
_UTR = _MODEL + _STREAMING + ("utr/core.py",)
_POSTQC = _COMPONENTS + _SELECTION + _UTR + ("postqc/core.py", "postqc/unified.py")

SOURCE_GROUPS = {
    # preflight.py imports deployment even in its input-only execution mode.
    "input_preflight": _PREFLIGHT,
    "preflight": _PREFLIGHT + _TE,
    "rna_import": _ALIGNMENTS,
    "predictor_container": ("contracts/deployment.py",) + _COMMON,
    "genemark": ("contracts/deployment.py", "predictor/commands.py") + _COMMON,
    "prediction_hints": ("contracts/deployment.py",) + _COMMON,
    "augustus_config": ("contracts/deployment.py",) + _COMMON,
    "augustus": ("contracts/deployment.py", "predictor/training.py") + _MODEL + _PREDICTION,
    "prediction_complete": _COMMON,
    "masking_rescue": ("contracts/deployment.py",) + _COMPONENTS + _MASKING,
    "junction_extract": _JUNCTIONS,
    "junction_aggregate": _JUNCTIONS + _COMMON,
    "models": _MODEL + _MASKING,
    "evidence": _COMPONENTS + _COMMON,
    "selection": _SELECTION + _JUNCTIONS + _COMMON,
    "coding_prepare": _UTR + _COMMON,
    "coding_group": _UTR + ("utr/grouping.py",),
    "postqc": _POSTQC,
    "utr_prepare": _UTR + _COMMON,
    "utr_refine": _UTR + ("utr/ranked.py",),
    "utr_aggregate": _UTR,
    "release": _POSTQC + ("release/core.py", "release/te.py"),
}


def source_inputs(repo_root, *groups):
    """Return sorted, existing modules and package initializers for the groups."""
    root = Path(repo_root).resolve()
    paths = {root / "geneinit/contracts/workflow.py"}
    for group in groups:
        for relative in SOURCE_GROUPS[group]:
            path = root / "geneinit" / relative
            paths.add(path)
            for directory in path.parents:
                if directory == root:
                    break
                paths.add(directory / "__init__.py")
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"workflow source input is missing: {path}")
    return sorted(str(path) for path in paths)
