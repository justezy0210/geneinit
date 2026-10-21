"""Thin Snakemake adapter for the shared preflight implementation."""

from pathlib import Path
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.preflight import run_input_preflight, run_preflight  # noqa: E402


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
try:
    output_dir = Path(snakemake.output.complete).parent
    if snakemake.params.input_only:
        summary = run_input_preflight(dict(snakemake.config), output_dir)
    else:
        summary = run_preflight(dict(snakemake.config), output_dir,
                                input_validation=Path(snakemake.input.inputs))
except Exception:
    log_path.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log_path.write_text(
        "GeneInIt preflight PASS\n"
        f"validation_sha256\t{summary.get('semantic_config_sha256', summary.get('input_sha256'))}\n",
        encoding="utf-8",
    )
