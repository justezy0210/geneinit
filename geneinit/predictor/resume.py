"""Validated checkpoints for resuming a stopped serial AUGUSTUS prediction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import shutil

from geneinit.predictor.common import PredictorError, sha256_file, validate_gff_coordinates
from geneinit.predictor.commands import augustus_prediction_command


SEQUENCE_HEADER = re.compile(
    r"^# ----- prediction on sequence number (\d+) \(length = (\d+), name = (\S+)\) -----"
)
METADATA_FILES = ("selection.json", "selection.tsv", "heldout_evaluation.json")


def prediction_signature() -> list[str]:
    return list(augustus_prediction_command(
        species="SPECIES", genome=Path("GENOME"), hints=Path("HINTS"),
        extrinsic_config=Path("EXTRINSIC"), augustus_config=Path("CONFIG"), output=Path("OUTPUT"),
    ))


def file_identity(path: Path) -> dict:
    path = path.resolve()
    before = path.stat()
    digest = sha256_file(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino) != (
        after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino
    ):
        raise PredictorError(f"prediction resume input changed while hashing: {path}")
    return {"path": str(path), "size": after.st_size, "mtime_ns": after.st_mtime_ns,
            "sha256": digest}


def model_hashes(directory: Path) -> dict:
    result = {str(p.relative_to(directory)): sha256_file(p)
              for p in sorted(directory.rglob("*")) if p.is_file()}
    if not result:
        raise PredictorError(f"empty trained species: {directory}")
    return result


def feature_hashes(path: Path, seqids) -> dict:
    digests = {seqid: hashlib.sha256() for seqid in seqids}
    with path.open("rb") as source:
        for line in source:
            if line.strip() and not line.startswith(b"#"):
                seqid = line.split(b"\t", 1)[0].decode()
                if seqid in digests:
                    digests[seqid].update(line)
    return {seqid: digest.hexdigest() for seqid, digest in digests.items()}


def extract_completed_prefix(raw: Path, lengths: dict[str, int], directory: Path) -> list[dict]:
    """Only the next reference header proves the previous reference finished.

    The last reference is always discarded, even if its last gene is complete.
    Source files must belong to a stopped process before calling this function.
    """
    directory.mkdir(parents=True, exist_ok=False)
    expected = list(lengths.items())
    preamble = []
    completed = []
    current = None
    handle = None
    active_gene = None
    genes = 0

    def finish():
        handle.close()
        if active_gene is not None:
            raise PredictorError("incomplete gene before next reference header")
        path = Path(current["output"])
        if genes:
            validate_gff_coordinates(path, {current["seqid"]: current["length"]})
        current.update(genes=genes, sha256=sha256_file(path),
                       feature_sha256=feature_hashes(path, [current["seqid"]])[current["seqid"]])
        completed.append(current)

    try:
        with raw.open(encoding="utf-8") as source:
            for line in source:
                match = SEQUENCE_HEADER.match(line)
                if match:
                    number, length, seqid = int(match[1]), int(match[2]), match[3]
                    if number != len(completed) + (2 if current else 1):
                        raise PredictorError("nonconsecutive serial prediction headers")
                    if number > len(expected) or expected[number-1] != (seqid, length):
                        raise PredictorError("serial prediction reference order/length mismatch")
                    if handle is not None:
                        finish()
                    elif not any(s.startswith("# This output was generated with AUGUSTUS") for s in preamble):
                        raise PredictorError("missing AUGUSTUS run header")
                    output = directory / f"{number:06d}.gff3"
                    current = {"seqid": seqid, "length": length, "output": str(output.resolve())}
                    handle = output.open("w", encoding="utf-8")
                    handle.writelines(preamble)
                    active_gene, genes = None, 0
                if handle is None:
                    preamble.append(line)
                    continue
                if line.startswith("# start gene "):
                    if active_gene is not None:
                        raise PredictorError("nested AUGUSTUS gene")
                    active_gene = line.split()[3]
                elif line.startswith("# end gene "):
                    if active_gene != line.split()[3]:
                        raise PredictorError("mismatched AUGUSTUS gene boundary")
                    active_gene = None
                    genes += 1
                elif line.strip() and not line.startswith("#"):
                    if active_gene is None:
                        raise PredictorError("AUGUSTUS feature outside complete gene block")
                handle.write(line)
    finally:
        if handle is not None and not handle.closed:
            handle.close()
    if current is not None:
        Path(current["output"]).unlink()
    if not completed:
        raise PredictorError("no completed reference sequences in serial prediction")
    return completed


def load_prediction_resume(
    manifest: Path, *, genome: Path, hints: Path, training_gtf: Path,
    extrinsic: Path, species: str, species_dir: Path, work: Path, settings: dict,
) -> dict:
    """Fail before altering the species or running tools if provenance differs."""
    data = json.loads(manifest.read_text())
    if data.get("schema_version") != 2:
        raise PredictorError(
            "prediction resume requires schema_version=2 with input SHA-256 hashes; "
            "recreate the checkpoint from verified original inputs or run prediction anew"
        )
    if data.get("species") != species:
        raise PredictorError("prediction resume schema/species mismatch")
    if data["settings"] != settings:
        raise PredictorError("prediction resume training settings mismatch")
    if data["prediction_signature"] != prediction_signature():
        raise PredictorError("prediction resume AUGUSTUS options changed")
    for label, path in (("genome", genome), ("hints", hints), ("training_gtf", training_gtf)):
        identity = data.get("inputs", {}).get(label, {})
        if not isinstance(identity, dict) or not re.fullmatch(r"[0-9a-f]{64}", str(identity.get("sha256", ""))):
            raise PredictorError(
                f"prediction resume input lacks a valid SHA-256: {label}; "
                "recreate the checkpoint from verified original inputs or run prediction anew"
            )
        if identity != file_identity(path):
            raise PredictorError(f"prediction resume input changed: {label}")
    if data["extrinsic_sha256"] != sha256_file(extrinsic):
        raise PredictorError("prediction resume extrinsic config changed")
    backup = Path(data["trained_species"])
    if model_hashes(backup) != data["model_sha256"]:
        raise PredictorError("prediction resume archived model changed")
    if species_dir.exists() and model_hashes(species_dir) != data["model_sha256"]:
        raise PredictorError("prediction resume current model changed")
    training = Path(data["training_metadata"])
    for name in METADATA_FILES:
        if sha256_file(training / name) != data["metadata_sha256"][name]:
            raise PredictorError(f"prediction resume training metadata changed: {name}")
    if not data["completed_sequences"]:
        raise PredictorError("prediction resume has no completed sequences")
    for sequence in data["completed_sequences"]:
        if sha256_file(Path(sequence["output"])) != sequence["sha256"]:
            raise PredictorError(f"prediction resume fragment changed: {sequence['seqid']}")
    if not species_dir.exists():
        shutil.copytree(backup, species_dir)
    destination = work / "training"
    destination.mkdir(parents=True, exist_ok=True)
    for name in METADATA_FILES:
        shutil.copyfile(training / name, destination / name)
    (work / "prediction_resume.json").write_text(json.dumps(data, indent=2) + "\n")
    return data


def require_preserved_features(output: Path, completed: list[dict]) -> None:
    expected = {item["seqid"]: item["feature_sha256"] for item in completed}
    if feature_hashes(output, expected) != expected:
        raise PredictorError("merged prediction changed preserved features, IDs, or scores")
