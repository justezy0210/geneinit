"""Predict complete reference sequences concurrently, without new gene boundaries."""

from __future__ import annotations

from collections import OrderedDict, deque
from dataclasses import dataclass
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import tempfile
import time
from typing import Callable, Sequence, TextIO

from geneinit.predictor.common import (
    PredictorError, require_nonempty, run_logged, validate_gff_coordinates,
)


@dataclass
class PredictionShard:
    seqid: str
    length: int
    fasta: Path
    hints: Path
    output: Path
    log: Path


def split_prediction_inputs(genome: Path, hints: Path, directory: Path) -> list[PredictionShard]:
    """Stream FASTA and hints; numeric filenames never derive from reference IDs."""
    directory.mkdir(parents=True, exist_ok=True)
    shards: list[PredictionShard] = []
    by_id: dict[str, PredictionShard] = {}
    handle = None
    try:
        with genome.open(encoding="utf-8") as source:
            for raw in source:
                if not raw.strip():
                    continue
                if raw.startswith(">"):
                    if handle is not None:
                        handle.close()
                    parts = raw[1:].split()
                    if not parts or parts[0] in by_id:
                        raise PredictorError("empty or duplicate genome FASTA ID")
                    prefix = directory / f"{len(shards) + 1:06d}"
                    shard = PredictionShard(parts[0], 0, prefix.with_suffix(".fa"),
                                            prefix.with_suffix(".hints.gff"),
                                            prefix.with_suffix(".gff3"),
                                            prefix.with_suffix(".log"))
                    shards.append(shard)
                    by_id[shard.seqid] = shard
                    shard.hints.write_text("", encoding="utf-8")
                    handle = shard.fasta.open("w", encoding="utf-8")
                elif handle is None:
                    raise PredictorError("genome sequence precedes FASTA header")
                else:
                    shards[-1].length += len(raw.strip())
                handle.write(raw if raw.endswith("\n") else raw + "\n")
    finally:
        if handle is not None:
            handle.close()
    if not shards or any(shard.length == 0 for shard in shards):
        raise PredictorError("empty genome or empty FASTA sequence")

    # Bound open file descriptors even for assemblies with many scaffolds.
    handles: OrderedDict[str, TextIO] = OrderedDict()
    try:
        with hints.open(encoding="utf-8") as source:
            for number, raw in enumerate(source, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[0] not in by_id:
                    raise PredictorError(f"invalid hints record or unknown sequence at line {number}")
                shard = by_id[fields[0]]
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise PredictorError(f"non-integer hints coordinates at line {number}") from error
                if start < 1 or end < start or end > shard.length:
                    raise PredictorError(f"hints coordinates outside reference at line {number}")
                if shard.seqid in handles:
                    output = handles.pop(shard.seqid)
                else:
                    if len(handles) == 32:
                        handles.popitem(last=False)[1].close()
                    output = shard.hints.open("a", encoding="utf-8")
                handles[shard.seqid] = output
                output.write(raw if raw.endswith("\n") else raw + "\n")
    finally:
        for output in handles.values():
            output.close()
    return shards


def run_prediction_jobs(
    shards: Sequence[PredictionShard], threads: int,
    command_for: Callable[[PredictionShard], Sequence[str]],
    progress_log: Path,
) -> int:
    """Bound active AUGUSTUS processes; cancel their groups if any worker fails."""
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise PredictorError("prediction threads must be a positive integer")
    pending = deque(sorted(shards, key=lambda shard: -shard.length))
    running: dict[subprocess.Popen, tuple[PredictionShard, TextIO]] = {}
    maximum_active = 0
    finished = 0
    previous_sigterm = signal.getsignal(signal.SIGTERM)

    def interrupted(signum, _frame):
        raise PredictorError(f"AUGUSTUS prediction interrupted by signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    try:
        with progress_log.open("w", encoding="utf-8") as progress:
            progress.write(f"workers\t{min(threads, len(shards))}\nsequences\t{len(shards)}\n")
            progress.flush()
            while pending or running:
                while pending and len(running) < threads:
                    shard = pending.popleft()
                    shard.output.unlink(missing_ok=True)
                    command = list(command_for(shard))
                    handle = shard.log.open("w", encoding="utf-8")
                    handle.write("command\t" + shlex.join(command) + "\n")
                    handle.flush()
                    try:
                        process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                                   start_new_session=True)
                    except BaseException:
                        handle.close()
                        raise
                    running[process] = (shard, handle)
                    maximum_active = max(maximum_active, len(running))
                for process, (shard, handle) in list(running.items()):
                    code = process.poll()
                    if code is None:
                        continue
                    handle.close()
                    if code:
                        raise PredictorError(f"AUGUSTUS failed ({code}) on {shard.seqid}; inspect {shard.log}")
                    require_nonempty(shard.output, "AUGUSTUS shard output")
                    del running[process]
                    finished += 1
                    progress.write(f"completed\t{finished}/{len(shards)}\t{shard.seqid}\n")
                    progress.flush()
                if running:
                    time.sleep(0.05)
    finally:
        # Includes exceptions during launch and SIGINT/SIGTERM, not only nonzero exits.
        for process in running:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for process, (_shard, handle) in running.items():
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            handle.close()
        signal.signal(signal.SIGTERM, previous_sigterm)
    return maximum_active


def predict_parallel(
    *, genome: Path, hints: Path, work: Path, output: Path, threads: int,
    command_for: Callable[[PredictionShard], Sequence[str]],
    wrap_command: Callable[[Sequence[str], Path], Sequence[str]],
    log: Path, completed_sequences: list[dict] | None = None,
) -> dict:
    directory = Path(tempfile.mkdtemp(prefix="prediction_shards.", dir=work))
    shards = split_prediction_inputs(genome, hints, directory)
    completed_sequences = completed_sequences or []
    if [(s.seqid, s.length) for s in shards[:len(completed_sequences)]] != [
        (s["seqid"], s["length"]) for s in completed_sequences
    ]:
        raise PredictorError("resumed sequences must be a complete reference prefix")
    for shard, saved in zip(shards, completed_sequences):
        shard.output = Path(saved["output"])
    maximum = run_prediction_jobs(shards[len(completed_sequences):], threads, command_for, log)
    concatenated = directory / "concatenated.gff3"
    with concatenated.open("wb") as destination:
        for shard in shards:
            with shard.output.open("rb") as source:
                # Preserve AUGUSTUS headers: the official merger uses them to
                # identify runs and renumbers all gene/transcript references.
                shutil.copyfileobj(source, destination)
    temporary = directory / "joined.gff3"
    join_command = (
        "/bin/bash", "-c", 'exec /usr/bin/perl "$1" < "$2" > "$3"',
        "geneinit-join-augustus", "/opt/Augustus/scripts/join_aug_pred.pl",
        str(concatenated), str(temporary),
    )
    run_logged(wrap_command(join_command, directory), directory / "join.log")
    require_nonempty(temporary, "joined AUGUSTUS output")
    validate_gff_coordinates(temporary, {shard.seqid: shard.length for shard in shards})
    if completed_sequences:
        from geneinit.predictor.resume import require_preserved_features
        require_preserved_features(temporary, completed_sequences)
    temporary.replace(output)
    report = {
        "method": "whole_reference_sequences", "threads": threads,
        "sequence_count": len(shards), "maximum_active_workers": maximum,
        "reused_sequence_count": len(completed_sequences),
        "predicted_sequence_count": len(shards) - len(completed_sequences),
        "shards_directory": str(directory),
        "sequences": [{"seqid": s.seqid, "length": s.length,
                       "output": str(s.output), "log": str(s.log)} for s in shards],
    }
    (work / "prediction_parallel.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
