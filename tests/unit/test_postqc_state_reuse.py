"""Preserve post-QC state transitions while avoiding redundant whole-set work."""

import csv
from copy import deepcopy

import pytest

import geneinit.postqc.core as postqc
from geneinit.selection.unified import UnifiedSelector
from test_postqc_core import node
from test_unified_postqc import make_inputs
from test_unified_selection import engine, snapshot, supported


def test_convergence_reuses_previous_digest_with_identical_state_receipts(monkeypatch):
    nodes = {"t1": node("t1", "C1", 100, 180)}
    original_hash = postqc._node_state_sha256
    expected = [original_hash(nodes)]
    future = deepcopy(nodes)
    future["t1"].retired = True
    expected.append(original_hash(future))
    future["t1"].active_fusion = True
    expected.append(original_hash(future))
    hashed = []

    def counted_hash(current):
        digest = original_hash(current)
        hashed.append(digest)
        return digest

    calls = 0

    def apply_iteration():
        nonlocal calls
        calls += 1
        if calls == 1:
            nodes["t1"].retired = True
            return [], [], [], [{"transcript_id": "t1"}]
        if calls == 2:
            nodes["t1"].active_fusion = True
        return [], [], [], []

    monkeypatch.setattr(postqc, "_node_state_sha256", counted_hash)
    _, _, retired, rows, converged = postqc.iterate_structural_qc(
        nodes, 4, apply_iteration,
    )
    assert converged == calls == 3
    assert hashed == [*expected, expected[-1]]
    assert [row["state_before_sha256"] for row in rows] == expected
    assert [row["state_after_sha256"] for row in rows] == expected[1:] + expected[-1:]
    assert [row["changed"] for row in rows] == [1, 1, 0]
    assert retired == [{"transcript_id": "t1"}]


def test_iteration_exception_propagates_without_an_after_digest(monkeypatch):
    nodes = {"t1": node("t1", "C1", 100, 180)}
    hashes = []
    monkeypatch.setattr(postqc, "_node_state_sha256", lambda current: hashes.append(current) or "state")

    def fail():
        raise RuntimeError("repair failed")

    with pytest.raises(RuntimeError, match="repair failed"):
        postqc.iterate_structural_qc(nodes, 3, fail)
    assert len(hashes) == 1


def test_invalid_iteration_limit_does_not_hash_or_apply(monkeypatch):
    def unexpected(*args):
        pytest.fail("invalid configuration must fail before doing work")

    monkeypatch.setattr(postqc, "_node_state_sha256", unexpected)
    with pytest.raises(postqc.PostQcError, match="must be positive"):
        postqc.iterate_structural_qc({}, 0, unexpected)


def test_selected_keys_tracks_updates_and_does_not_expose_internal_state():
    selector = engine([supported("a"), supported("far", 2000)])
    keys = selector.selected_keys
    assert keys == {record.key for record in selector.selected}
    keys.clear()
    assert selector.selected_keys == {"a", "far"}
    selector.update([], removed=["a"])
    assert selector.selected_keys == {"far"}
    selector.update([supported("new")])
    assert selector.selected_keys == {"new", "far"}
    assert snapshot(selector) == snapshot(engine([supported("new"), supported("far", 2000)]))


def test_postqc_membership_never_sorts_selected_records(tmp_path, monkeypatch):
    args, _, _ = make_inputs(tmp_path, "unified_v10_single_exon50_corroborated_multisample.json")

    def unnecessary_sort(selector):
        pytest.fail("post-QC only needs selected candidate membership")

    monkeypatch.setattr(UnifiedSelector, "selected", property(unnecessary_sort))
    postqc.run_postqc(*args)
    with (tmp_path / "postqc_summary.tsv").open() as handle:
        summary = {row["metric"]: row["value"]
                   for row in csv.DictReader(handle, delimiter="\t")}
    assert summary["output_transcript_count"] == "1"
    assert summary["retired_transcript_count"] == "2"
    assert summary["converged_iteration"] == "2"
