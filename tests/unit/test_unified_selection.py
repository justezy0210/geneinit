from copy import deepcopy
import itertools
import json
from pathlib import Path

import pytest

from geneinit.selection.core import SelectionError, audit_row
from geneinit.selection.unified import UnifiedSelector
from test_selection_core import make_record as _make_record


POLICY = json.loads(
    (Path(__file__).resolve().parents[2] / "policies/unified_v3.json").read_text()
)["selector_core"]


def make_record(*args, **kwargs):
    record = _make_record(*args, **kwargs)
    record.row["protein_phase_comparable_parts"] = "0"
    record.row["protein_chain_id"] = "chain" if kwargs.get("protein_tier", 0) else "NA"
    return record


def engine(records=()):
    result = UnifiedSelector("ASM", POLICY, sample_count=6, replicates_per_condition=2)
    result.update(list(records))
    return result


def snapshot(selector):
    return {key: audit_row(record, POLICY["policy_id"])
            for key, record in selector.records.items()}


def supported(key, offset=0):
    return make_record(key, ((100+offset, 150+offset, 0),
                             (201+offset, 300+offset, 0)))


def test_hard_qc_still_excludes_doddi_like_short_single_exon():
    record = make_record("72aa", ((100, 318, 0),), protein_aa=72,
                         protein_tier=4, protein_unique=True,
                         protein_precision=1, terminals=2)
    selector = engine([record])
    assert not selector.selected
    assert selector.records["72aa"].reasons == ["HARD_QC_FAIL"]


def test_supplemental_evidence_is_admitted_without_a_second_selector():
    record = make_record("full", ((100, 150, 0), (201, 250, 0), (301, 400, 0)),
                         hq=[1, 1], samples=[1, 1], replicate=[0, 0], rna_tier=2)
    selector = engine([record])
    assert [r.key for r in selector.selected] == ["full"]
    assert selector.records["full"].routes == ["SUPPORTED_FULL_RNA"]
    assert not selector.records["full"].base_supplement


def test_supplemental_models_do_not_grandfather_conflicts():
    records = [
        make_record(key, ((start, 150, 0), (201, 250, 0), (301, 400, 0)),
                    hq=[1, 1], samples=[1, 1], replicate=[0, 0], rna_tier=2)
        for key, start in (("a", 100), ("b", 90))
    ]
    selector = engine(records)
    assert len(selector.selected) == 1
    assert sum(r.decision == "REJECT_INCOMPATIBLE_GRAPH"
               for r in selector.records.values()) == 1


def test_te_and_general_models_never_block_each_other():
    general, te = supported("gene"), supported("te")
    general.row["candidate_pool"] = "GENERAL"
    te.row["candidate_pool"] = "TE"
    selector = engine([te, general])
    assert [r.key for r in selector.selected] == ["gene"]
    assert selector.records["te"].decision == "TE_SELECTED"


def test_no_change_reuses_evidence_and_unrelated_loci():
    selector = engine([supported("a"), supported("b", 1000), supported("c", 2000)])
    before = snapshot(selector)
    stats = selector.update([supported("a")])
    assert stats.evaluated_candidates == 0
    assert stats.reused_candidates == 1
    assert stats.untouched_candidates == 2
    assert snapshot(selector) == before


def test_new_competitor_invalidates_neighbors_but_reuses_their_evidence():
    a, far = supported("a"), supported("far", 2000)
    selector = engine([a, far])
    new = make_record("new", ((90, 150, 0), (201, 320, 0)), protein_aa=240)
    stats = selector.update([new])
    assert stats.evaluated_candidates == 1
    assert stats.reused_candidates == 1
    assert stats.untouched_candidates == 1
    assert snapshot(selector) == snapshot(engine([a, far, new]))


def test_removal_can_restore_a_previously_rejected_neighbor():
    first = supported("a")
    other = make_record("b", ((90, 150, 0), (201, 320, 0)), protein_aa=240)
    selector = engine([first, other, supported("far", 2000)])
    winner = selector.selected[0].key
    remaining = [r for r in [first, other, supported("far", 2000)] if r.key != winner]
    selector.update([], removed=[winner])
    assert snapshot(selector) == snapshot(engine(remaining))


def test_moving_a_candidate_rechecks_both_old_and_new_loci():
    records = [supported("a"), supported("b"), supported("c", 1000),
               supported("d", 5000)]
    selector = engine(records)
    moved = supported("a", 1000)
    selector.update([moved])
    assert selector.stats.untouched_candidates == 1
    assert snapshot(selector) == snapshot(engine([moved] + records[1:]))


def test_changed_evidence_invalidates_only_the_changed_candidates_metrics():
    a, b = supported("a"), supported("b")
    selector = engine([a, b])
    updated = deepcopy(a)
    updated.row["orf_status"] = "missing_start"
    selector.update([updated])
    assert selector.stats.evaluated_candidates == 1
    assert selector.stats.reused_candidates == 1
    assert snapshot(selector) == snapshot(engine([updated, b]))


def test_order_does_not_change_fusion_or_selection_decisions():
    records = [
        supported("left"),
        supported("right", 1000),
        make_record("center", ((90, 150, 0), (201, 300, 0),
                              (1100, 1150, 0), (1201, 1320, 0))),
    ]
    expected = snapshot(engine(records))
    for order in itertools.permutations(records):
        assert snapshot(engine(order)) == expected


def test_repeated_updates_do_not_accumulate_old_decisions():
    a, b = supported("a"), supported("b")
    selector = engine([a, b])
    for _ in range(3):
        selector.update([a, b])
    assert snapshot(selector) == snapshot(engine([a, b]))


def test_rejects_ambiguous_update_requests_and_implicit_policy():
    selector = engine([supported("a")])
    with pytest.raises(SelectionError, match="duplicate"):
        selector.update([supported("b"), supported("b")])
    with pytest.raises(SelectionError, match="removed and updated"):
        selector.update([supported("a")], removed=["a"])
    with pytest.raises(SelectionError, match="unknown"):
        selector.update([], removed=["absent"])
    policy = deepcopy(POLICY)
    del policy["engine"]
    with pytest.raises(SelectionError, match="explicit"):
        UnifiedSelector("ASM", policy, sample_count=6, replicates_per_condition=2)


def test_restored_checkpoint_matches_a_fresh_selector_without_recomputing(monkeypatch):
    import geneinit.selection.unified as unified
    records = [supported("a"), supported("b"), supported("far", 2000)]
    original = engine(records)
    restored = engine()
    audits = snapshot(original)
    expected = engine(records + [supported("new", 2000)])
    real_evaluate = unified.evaluate_record
    calls = []

    def counted(record, *args, **kwargs):
        calls.append(record.key)
        return real_evaluate(record, *args, **kwargs)

    monkeypatch.setattr(unified, "evaluate_record", counted)
    restored.restore(records, audits)
    assert calls == []
    assert snapshot(restored) == snapshot(original)
    restored.update([supported("new", 2000)])
    assert calls == ["new"]
    assert snapshot(restored) == snapshot(expected)
