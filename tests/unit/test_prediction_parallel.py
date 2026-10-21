import json
from pathlib import Path
import sys
import time

import pytest

from geneinit.predictor.common import PredictorError
from geneinit.predictor.parallel import PredictionShard, run_prediction_jobs, split_prediction_inputs


def test_partition_preserves_masking_coordinates_and_hint_order(tmp_path):
    genome = tmp_path / "input.fa"
    genome.write_text(">chr2 description\nAcgTN\n>../chr1\nAaTtCc\n>empty_hints\nGGG\n")
    hints = tmp_path / "input.gff"
    first = "chr2\tE\tintron\t2\t4\t.\t+\t.\tgrp=a;src=E;\n"
    second = "../chr1\tP\tCDSpart\t1\t3\t.\t-\t0\tgrp=b;src=P;\n"
    third = "chr2\tRM\tnonexonpart\t1\t5\t.\t.\t.\tsrc=RM;\n"
    hints.write_text("# ignored comment\n" + first + second + third)
    shards = split_prediction_inputs(genome, hints, tmp_path / "shards")
    assert [s.seqid for s in shards] == ["chr2", "../chr1", "empty_hints"]
    assert [s.length for s in shards] == [5, 6, 3]
    assert shards[0].fasta.read_text() == ">chr2 description\nAcgTN\n"
    assert shards[0].hints.read_text() == first + third
    assert shards[1].hints.read_text() == second
    assert shards[2].hints.read_text() == ""
    assert all(s.fasta.parent == tmp_path / "shards" for s in shards)


@pytest.mark.parametrize("genome_text,hint_text", [
    (">x\nAAA\n>x\nCCC\n", ""), (">x\n", ""), ("AAA\n", ""),
    (">x\nAAA\n", "missing\tE\tintron\t1\t2\t.\t+\t.\tsrc=E;\n"),
    (">x\nAAA\n", "x\tE\tintron\t1\t4\t.\t+\t.\tsrc=E;\n"),
])
def test_partition_rejects_invalid_input(tmp_path, genome_text, hint_text):
    genome, hints = tmp_path / "genome.fa", tmp_path / "hints.gff"
    genome.write_text(genome_text)
    hints.write_text(hint_text)
    with pytest.raises(PredictorError):
        split_prediction_inputs(genome, hints, tmp_path / "shards")


def make_shards(directory, number):
    return [PredictionShard(f"chr{i}", number - i, directory / f"{i}.fa",
                            directory / f"{i}.hints", directory / f"{i}.gff",
                            directory / f"{i}.log") for i in range(number)]


@pytest.mark.parametrize("budget", [1, 2, 34])
def test_actual_process_overlap_respects_budget(tmp_path, budget):
    shards = make_shards(tmp_path, max(4, budget + 1))
    events = tmp_path / "events.jsonl"
    worker = tmp_path / "worker.py"
    worker.write_text(
        "import fcntl,json,sys,time\nfrom pathlib import Path\n"
        "def event(change):\n"
        " with open(sys.argv[1], 'a') as f:\n"
        "  fcntl.flock(f, fcntl.LOCK_EX)\n"
        "  f.write(json.dumps({'t':time.monotonic(),'change':change})+'\\n')\n"
        "event(1)\ntime.sleep(0.2)\n"
        "Path(sys.argv[2]).write_text('# finished\\n')\nevent(-1)\n"
    )
    maximum = run_prediction_jobs(shards, budget,
        lambda s: [sys.executable, str(worker), str(events), str(s.output)],
        tmp_path / "progress.log")
    active = peak = 0
    for event in sorted(map(json.loads, events.read_text().splitlines()), key=lambda e: e["t"]):
        active += event["change"]
        peak = max(peak, active)
    assert active == 0
    assert peak <= budget
    assert maximum == budget
    if budget > 1:
        assert peak > 1
    assert all(s.output.exists() for s in shards)


def test_worker_failure_stops_running_children_and_pending_jobs(tmp_path):
    shards = make_shards(tmp_path, 4)
    worker = tmp_path / "worker.py"
    worker.write_text(
        "import subprocess,sys,time\nfrom pathlib import Path\n"
        "if sys.argv[1]=='chr0':\n"
        " p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'])\n"
        " Path(sys.argv[2]).write_text(str(p.pid))\n time.sleep(30)\n"
        "else:\n time.sleep(0.3)\n sys.exit(7)\n"
    )
    child_pid = tmp_path / "child.pid"
    started = time.monotonic()
    with pytest.raises(PredictorError, match=r"failed \(7\)"):
        run_prediction_jobs(shards, 2,
            lambda s: [sys.executable, str(worker), s.seqid, str(child_pid)],
            tmp_path / "progress.log")
    assert time.monotonic() - started < 5
    assert not shards[2].log.exists()
    stat = Path("/proc") / child_pid.read_text() / "stat"
    assert not stat.exists() or stat.read_text().split(") ")[1].startswith("Z")
