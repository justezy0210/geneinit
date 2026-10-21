"""Shared deterministic predictor helpers."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Sequence


class PredictorError(RuntimeError):
    """Raised when a predictor stage violates its completion contract."""


def require_nonempty(path: Path, label: str | None = None) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise PredictorError(f"missing or empty {label or 'file'}: {path}")
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_copy(source: Path, target: Path) -> None:
    require_nonempty(source, "source file")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp.{os.getpid()}")
    shutil.copy2(source, temporary)
    if sha256_file(source) != sha256_file(temporary):
        temporary.unlink(missing_ok=True)
        raise PredictorError(f"staged copy checksum mismatch: {source} -> {target}")
    temporary.replace(target)


def run_logged(command: Sequence[str], log: Path, cwd: Path | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("command\t" + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(
            command,
            cwd=cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        raise PredictorError(
            f"command failed with exit code {completed.returncode}; inspect {log}"
        )


def write_completion(path: Path, rows: Iterable[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("field", "value"))
        writer.writerows(rows)
    temporary.replace(path)


def validate_gff_coordinates(path: Path, lengths: dict[str, int]) -> int:
    require_nonempty(path, "GFF/GTF")
    count = 0
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise PredictorError(f"{path}:{number} does not have 9 columns")
            if fields[0] not in lengths:
                raise PredictorError(f"{path}:{number} seqid absent from genome: {fields[0]}")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise PredictorError(f"{path}:{number} has non-integer coordinates") from error
            if start < 1 or end < start or end > lengths[fields[0]]:
                raise PredictorError(
                    f"{path}:{number} coordinate out of bounds: {fields[0]}:{start}-{end}"
                )
            if fields[6] not in {"+", "-", "."}:
                raise PredictorError(f"{path}:{number} has invalid strand: {fields[6]}")
            count += 1
    if count == 0:
        raise PredictorError(f"no feature records found: {path}")
    return count


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split(maxsplit=1)[0]
                if current in lengths:
                    raise PredictorError(f"duplicate genome FASTA ID: {current}")
                lengths[current] = 0
            elif current is None:
                raise PredictorError("genome sequence precedes FASTA header")
            else:
                lengths[current] += len(line)
    if not lengths:
        raise PredictorError(f"no genome records: {path}")
    return lengths
