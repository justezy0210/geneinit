"""Thin Snakemake adapter for the shared preflight implementation."""

from pathlib import Path
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.preflight import run_preflight  # noqa: E402


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
try:
    output_dir = Path(snakemake.output.complete).parent
    summary = run_preflight(dict(snakemake.config), output_dir)
except Exception:
    log_path.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log_path.write_text(
        "GeneInIt preflight PASS\n"
        f"semantic_config_sha256\t{summary['semantic_config_sha256']}\n",
        encoding="utf-8",
    )
