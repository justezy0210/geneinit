"""Launch the local workflow with one explicit CPU budget."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import sys

from geneinit import __version__


def positive_threads(value: str) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("threads must be a positive integer") from error
    if number < 1:
        raise argparse.ArgumentTypeError("threads must be a positive integer")
    return number


def parse_arguments(argv: list[str]) -> tuple[int | None, list[str]]:
    parser = argparse.ArgumentParser(
        prog="geneinit-run", allow_abbrev=False,
        description="Run GeneInIt with a total CPU budget shared by all jobs.",
        epilog="Remaining arguments are passed to Snakemake, e.g. --configfile config.yaml.",
    )
    parser.add_argument(
        "--threads", "--cores", "-c", type=positive_threads, action="append",
        metavar="N", help="maximum CPUs for the entire run (required; --cores/-c are aliases)",
    )
    parsed, remaining = parser.parse_known_args(argv)
    values = parsed.threads or []
    if not values and remaining == ["--version"]:
        return None, remaining
    if not values:
        parser.error("specify --threads N (or --cores N)")
    if len(set(values)) != 1:
        parser.error("conflicting thread/CPU budgets")
    return values[0], remaining


def limit_cpus(requested: int) -> tuple[int, list[int] | None]:
    """Preserve scheduler restrictions and pass the CPU mask to all descendants."""
    if hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"):
        cpus = sorted(os.sched_getaffinity(0))[:requested]
        if not cpus:
            raise RuntimeError("no CPUs are available in the current affinity mask")
        os.sched_setaffinity(0, cpus)
        return len(cpus), cpus
    return min(requested, os.cpu_count() or 1), None


def main() -> None:
    requested, remaining = parse_arguments(sys.argv[1:])
    if remaining == ["--version"]:
        print(f"GeneInIt {__version__}")
        return
    root = Path(__file__).resolve().parents[1]
    executable = shutil.which("snakemake")
    if executable is None:
        sibling = Path(sys.executable).with_name("snakemake")
        if sibling.is_file() and os.access(sibling, os.X_OK):
            executable = str(sibling)
    if executable is None:
        raise SystemExit("geneinit-run: snakemake is not installed in this environment")
    command = [executable, "--snakefile", str(root / "Snakefile")]
    profile = "allinone" if root == Path("/opt/GeneInIt") else "local"
    command += ["--profile", str(root / "profiles" / profile)]
    if requested is not None:
        effective, cpus = limit_cpus(requested)
        # Tool-specific worker arguments control parallelism. Prevent each worker
        # from also creating a full BLAS/OpenMP pool.
        for variable in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
        ):
            os.environ[variable] = "1"
        print(
            f"GeneInIt CPU budget: requested={requested}, effective={effective}, "
            f"affinity={','.join(map(str, cpus)) if cpus is not None else 'unavailable'}",
            file=sys.stderr, flush=True,
        )
        command += ["--cores", str(effective)]
    os.execv(executable, [*command, *remaining])


if __name__ == "__main__":
    main()
