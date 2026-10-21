import json
import os
from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest

from geneinit.predictor.common import PredictorError, sha256_file
from geneinit.predictor.resume import (
    METADATA_FILES, extract_completed_prefix, file_identity, load_prediction_resume,
    model_hashes, prediction_signature, require_preserved_features,
)


def header(number, seqid):
    return f"# ----- prediction on sequence number {number} (length = 10, name = {seqid}) -----\n"


def gene(seqid, number):
    return (f"# start gene g{number}\n{seqid}\tAUGUSTUS\tgene\t1\t9\t0.8\t+\t.\tID=g{number}\n"
            f"# end gene g{number}\n")


def raw_fixture(tmp_path):
    raw = tmp_path / "serial.gff3"
    raw.write_text("##gff-version 3\n# This output was generated with AUGUSTUS (version 3.5.0).\n"
                   + header(1, "chr1") + gene("chr1", 1) + header(2, "empty")
                   + header(3, "unfinished") + "# start gene g2\nunfinished\tAUGUSTUS")
    return raw


def test_extract_only_proven_complete_sequences_including_empty(tmp_path):
    raw = raw_fixture(tmp_path)
    before = raw.read_bytes()
    result = extract_completed_prefix(raw, {"chr1": 10, "empty": 10, "unfinished": 10}, tmp_path / "out")
    assert [item["seqid"] for item in result] == ["chr1", "empty"]
    assert [item["genes"] for item in result] == [1, 0]
    assert raw.read_bytes() == before
    assert not (tmp_path / "out/000003.gff3").exists()
    require_preserved_features(raw, result)
    raw.write_text(raw.read_text().replace("0.8", "0.7"))
    with pytest.raises(PredictorError, match="preserved features"):
        require_preserved_features(raw, result)


@pytest.mark.parametrize("old,new", [("name = chr1", "name = wrong"),
    ("length = 10", "length = 11"), ("# end gene g1", "# end gene g9"),
    ("# end gene g1\n", ""), ("sequence number 2", "sequence number 4")])
def test_extraction_rejects_inconsistent_serial_output(tmp_path, old, new):
    raw = raw_fixture(tmp_path)
    raw.write_text(raw.read_text().replace(old, new))
    with pytest.raises(PredictorError):
        extract_completed_prefix(raw, {"chr1": 10, "empty": 10, "unfinished": 10}, tmp_path / "out")


def checkpoint(tmp_path):
    raw = raw_fixture(tmp_path)
    completed = extract_completed_prefix(raw, {"chr1": 10, "empty": 10, "unfinished": 10}, tmp_path / "out")
    inputs = {}
    for name in ("genome", "hints", "training_gtf", "extrinsic"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name)
    backup = tmp_path / "backup"
    backup.mkdir()
    (backup / "test_parameters.cfg").write_text("trained")
    training = tmp_path / "metadata"
    training.mkdir()
    for name in METADATA_FILES:
        (training / name).write_text("{}")
    (training / "selection.json").write_text(json.dumps({"splits": {"onlytrain": {"loci": 5000}, "optimize": {"loci": 500}}}))
    (training / "heldout_evaluation.json").write_text('{"initial": {}, "optimized": {}}')
    settings = {"rounds": 5, "kfold": 34, "training_flank": 1000, "training_selection": {}, "genetic_code": 1}
    data = {"schema_version": 2, "species": "test", "settings": settings,
        "inputs": {name: file_identity(inputs[name]) for name in ("genome", "hints", "training_gtf")},
        "extrinsic_sha256": sha256_file(inputs["extrinsic"]), "trained_species": str(backup),
        "model_sha256": model_hashes(backup), "training_metadata": str(training),
        "metadata_sha256": {name: sha256_file(training/name) for name in METADATA_FILES},
        "completed_sequences": completed, "prediction_signature": prediction_signature(),
        "kfold": 34, "optimization_cpus": 34}
    manifest = tmp_path / "resume.json"
    manifest.write_text(json.dumps(data))
    work = tmp_path / "work"
    work.mkdir()
    args = dict(**inputs, species="test", species_dir=tmp_path / "config/species/test", work=work, settings=settings)
    return manifest, args


@pytest.mark.parametrize("changed", ["genome", "hints", "training_gtf"])
def test_resume_rejects_changed_content_with_preserved_size_and_mtime(tmp_path, changed):
    manifest, args = checkpoint(tmp_path)
    path = args[changed]
    before = path.stat()
    path.write_text(path.read_text().upper())
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert path.stat().st_size == before.st_size
    assert path.stat().st_mtime_ns == before.st_mtime_ns
    with pytest.raises(PredictorError, match=f"input changed: {changed}"):
        load_prediction_resume(manifest, **args)
    assert not args["species_dir"].exists()


@pytest.mark.parametrize("legacy", ["schema", "missing_hash"])
def test_resume_refuses_unhashed_legacy_checkpoints_with_recovery_guidance(tmp_path, legacy):
    manifest, args = checkpoint(tmp_path)
    data = json.loads(manifest.read_text())
    if legacy == "schema":
        data["schema_version"] = 1
    else:
        del data["inputs"]["genome"]["sha256"]
    manifest.write_text(json.dumps(data))
    with pytest.raises(PredictorError, match="recreate the checkpoint from verified original inputs"):
        load_prediction_resume(manifest, **args)
    assert not args["species_dir"].exists()


@pytest.mark.parametrize("changed", ["genome", "extrinsic", "model", "metadata", "fragment", "settings"])
def test_resume_refuses_changed_provenance_before_creating_model(tmp_path, changed):
    manifest, args = checkpoint(tmp_path)
    data = json.loads(manifest.read_text())
    if changed in ("genome", "extrinsic"):
        args[changed].write_text("changed")
    elif changed == "model":
        (Path(data["trained_species"]) / "test_parameters.cfg").write_text("changed")
    elif changed == "metadata":
        (Path(data["training_metadata"]) / "selection.tsv").write_text("changed")
    elif changed == "fragment":
        Path(data["completed_sequences"][0]["output"]).write_text("changed")
    else:
        args["settings"]["kfold"] = 8
    with pytest.raises(PredictorError):
        load_prediction_resume(manifest, **args)
    assert not args["species_dir"].exists()


def test_resume_copies_metadata_without_writable_archive_alias(tmp_path):
    manifest, args = checkpoint(tmp_path)
    data = load_prediction_resume(manifest, **args)
    assert model_hashes(args["species_dir"]) == data["model_sha256"]
    (args["work"] / "training/selection.tsv").write_text("modified work")
    assert (Path(data["training_metadata"]) / "selection.tsv").read_text() == "{}"
    (args["species_dir"] / "test_parameters.cfg").write_text("changed")
    with pytest.raises(PredictorError, match="current model changed"):
        load_prediction_resume(manifest, **args)


def test_workflow_resume_reaches_prediction_without_any_training_command(tmp_path, monkeypatch):
    manifest, args = checkpoint(tmp_path)
    import geneinit.contracts.deployment as deployment
    import geneinit.predictor.common as common
    import geneinit.predictor.parallel as parallel
    monkeypatch.setattr(deployment, "load_predictor_deployment", lambda config: object())
    def no_training(*args, **kwargs):
        pytest.fail("resume invoked a training command")
    monkeypatch.setattr(common, "run_logged", no_training)
    class ReachedPrediction(Exception):
        pass
    def prediction(**kwargs):
        assert len(kwargs["completed_sequences"]) == 2
        assert kwargs["threads"] == 34
        raise ReachedPrediction()
    monkeypatch.setattr(parallel, "predict_parallel", prediction)
    fake = SimpleNamespace(config={}, threads=34, log=[str(tmp_path/"rule.log")],
        input=SimpleNamespace(genome=str(args["genome"]), training=str(args["training_gtf"]),
            hints=str(args["hints"]), extrinsic_config=str(args["extrinsic"]),
            augustus_config=str(args["species_dir"].parents[1]), prediction_resume=[str(manifest)]),
        output=SimpleNamespace(work=str(args["work"])),
        params=SimpleNamespace(species="test", rounds=5, kfold=34, training_flank=1000, training_selection={}))
    script = Path(__file__).resolve().parents[2]/"workflow/scripts/run_augustus_etp.py"
    with pytest.raises(ReachedPrediction):
        runpy.run_path(str(script), init_globals={"snakemake": fake})
