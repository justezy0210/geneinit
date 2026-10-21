from copy import deepcopy
import itertools
import json
from pathlib import Path

import pytest

from geneinit.evidence.components import (
    CandidateRecord, JunctionAggregate, RnaTranscript, format_value,
)
from geneinit.evidence.rna_cds import rna_cds_metrics
from geneinit.selection.core import SelectionError, audit_row
from geneinit.selection.unified import UnifiedSelector
from test_rna_cds_boundary import frozen_jinju_records, ranked_model
from test_rna_competition import eight_records
from test_unified_selection import make_record


ROOT = Path(__file__).resolve().parents[2]
MODES = ("minimal", "coherent", "protein_first")
POLICIES = {
    mode: json.loads((ROOT / f"policies/unified_v6_single_exon50_{mode}.json").read_text())[
        "selector_core"
    ] for mode in MODES
}


def select(records=(), mode="coherent", policy=None):
    selector = UnifiedSelector("ASM", policy or POLICIES[mode],
                               sample_count=65, replicates_per_condition=1)
    selector.update(records)
    return selector


def snapshot(selector):
    return {key: audit_row(record, selector._policy["policy_id"])
            for key, record in selector.records.items()}


def boundary(record, *, full=False, chain=False, conflict=False):
    candidate = CandidateRecord(
        record.key, record.row["structure_id"], record.seqid, record.strand,
        record.cds, ("augustus",), sum(part.length for part in record.cds), "complete",
    )
    junctions = {
        (record.seqid, start, end, record.strand): JunctionAggregate(hq, samples, reps, 0)
        for (start, end), hq, samples, reps in zip(
            record.introns, json.loads(record.row["rna_hq_read_counts"]),
            json.loads(record.row["rna_sample_hq_ge1_counts"]),
            json.loads(record.row["rna_condition_ge2rep_hq_ge1_counts"]),
        )
    }
    transcripts = []
    if full or chain:
        exons = list(record.intervals)
        if not full:
            exons[0] = exons[0][0] + 3, exons[0][1]
            exons[-1] = exons[-1][0], exons[-1][1] - 3
        transcripts.append(RnaTranscript("merged_witness", record.seqid, record.strand, tuple(exons)))
    if conflict:
        start, end = record.end - 1, record.end + 30
        junctions[record.seqid, start, end, record.strand] = JunctionAggregate(8, 2, 1, 0)
        transcripts.append(RnaTranscript(
            "opposing_rna", record.seqid, record.strand,
            ((record.start, start - 1), (end + 1, end + 50)),
        ))
    record.row.update({key: format_value(value) for key, value in
                       rna_cds_metrics(candidate, transcripts, junctions).items()})
    return record


def spliced(key, *, protein=False, full=False, chain=False):
    return boundary(make_record(
        key, ((100, 198, 0), (300, 398, 0), (500, 598, 0)),
        protein_tier=4 if protein else 0, protein_unique=protein,
        protein_precision=1 if protein else 0,
        protein_intron_precision=1 if protein else 0,
        protein_intron_recall=1 if protein else 0,
        terminals=2 if protein else 0, protein_aa=98,
    ), full=full, chain=chain)


@pytest.mark.parametrize("mode", MODES)
def test_single_exon_complete_model_is_not_displaced_by_rna_only_fragment(mode):
    complete = ranked_model("complete", 430, full_rna=True)
    fragment = spliced("fragment", full=True)
    selector = select([fragment, complete], mode)
    assert {record.key for record in selector.selected} == {"complete"}
    observed = selector.records["complete"]
    assert not observed.strong_rna  # Admission remains intron-based.
    route = json.loads(snapshot(selector)["complete"]["selection_rank_json"])[0]
    assert route == (-8 if mode == "protein_first" else -9)
    if mode != "protein_first":
        assert "SINGLE_EXON_RNA_PROTEIN_DUAL_RANK" in observed.reasons


@pytest.mark.parametrize("mode", MODES)
def test_partial_rna_admission_remains_permissive_but_dual_bonus_is_explicit(mode):
    record = spliced("partial", protein=True)
    record.row["rna_hq_read_counts"] = "[8,0]"
    record.row["rna_sample_hq_ge1_counts"] = "[2,0]"
    record.row["rna_condition_ge2rep_hq_ge1_counts"] = "[1,0]"
    selector = select([record], mode)
    result = selector.records["partial"]
    assert result.selected and result.strong_rna and result.strong_unique_full_protein
    assert "RNA50_ADMISSION_WITHOUT_WHOLE_RNA_WITNESS" in result.reasons
    assert json.loads(snapshot(selector)["partial"]["selection_rank_json"])[0] == (
        -8 if mode == "coherent" else -9
    )


@pytest.mark.parametrize("mode", MODES)
def test_coding_chain_can_support_dual_rank_without_claiming_complete_ends(mode):
    selector = select([spliced("chain", protein=True, chain=True)], mode)
    result = selector.records["chain"]
    assert result.row["rna_cds_support_transcripts_json"] == "[]"
    assert result.row["rna_cds_self_strong_chain"] == "1"
    assert json.loads(snapshot(selector)["chain"]["selection_rank_json"])[0] == -9


@pytest.mark.parametrize("mode", MODES)
def test_better_supported_shorter_model_can_win(mode):
    shorter = spliced("shorter", protein=True, full=True)
    longer = ranked_model("longer", 430)
    assert [record.key for record in select([longer, shorter], mode).selected] == ["shorter"]


def fusion_models(*, full=False, chain=False, conflict=False):
    left = boundary(make_record(
        "left", ((100, 198, 0), (300, 398, 0)),
        protein_tier=4, protein_unique=True, protein_precision=1,
        protein_intron_precision=1, protein_intron_recall=1, terminals=2,
    ), full=True)
    right = boundary(make_record(
        "right", ((1000, 1098, 0), (1200, 1298, 0)),
        protein_tier=4, protein_unique=True, protein_precision=1,
        protein_intron_precision=1, protein_intron_recall=1, terminals=2,
    ), full=True)
    center = boundary(make_record(
        "center", ((90, 198, 0), (300, 500, 0), (900, 1098, 0), (1200, 1320, 0)),
        protein_tier=2, protein_unique=True, protein_precision=1,
        protein_intron_precision=1, protein_intron_recall=1, terminals=1,
    ), full=full, chain=chain, conflict=conflict)
    return [left, right, center]


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("support,precedence", [
    ({}, False), ({"full": True}, True), ({"chain": True}, True),
    ({"chain": True, "conflict": True}, False),
])
def test_fusion_precedence_needs_whole_chain_support_without_unopposed_conflict(
    mode, support, precedence,
):
    records = fusion_models(**support)
    selector = select(records, mode)
    center = selector.records["center"]
    assert center.active_fusion and center.eligible
    assert center.fusion_center_precedence is precedence
    assert (center.selected is precedence)
    if not precedence:
        assert "FUSION_PRECEDENCE_WITHHELD_WHOLE_RNA_OR_CONFLICT" in center.reasons
    expected = snapshot(selector)
    for order in itertools.permutations(records):
        assert snapshot(select(order, mode)) == expected


@pytest.mark.parametrize("mode", MODES)
def test_focal_and_eight_loci_replay_restore_and_update_agree(mode):
    records, _ = eight_records()
    records += frozen_jinju_records()
    original = select(records, mode)
    assert original.records["GCAN-f25132743f098767e8652d6a"].selected
    assert not original.records["GCAN-d77d0aa4cf68a62558eb51bc"].selected
    assert not original.records["GCAN-26e832edfbbd3413438ba65a"].selected
    expected = snapshot(original)
    assert snapshot(select(list(reversed(records)), mode)) == expected
    restored = select(mode=mode)
    restored.restore(records, expected)
    assert snapshot(restored) == expected
    assert all(record.whole_model_ranking == mode for record in restored.records.values())
    changed = deepcopy(records[0])
    changed.row["orf_status"] = "missing_start"
    restored.update([changed])
    assert snapshot(restored) == snapshot(select([changed, *records[1:]], mode))


@pytest.mark.parametrize("value", [None, "", "typo", False, 1])
def test_unknown_experiment_fails_closed(value):
    policy = deepcopy(POLICIES["minimal"])
    policy["whole_model_ranking"] = value
    with pytest.raises(SelectionError, match="unknown whole_model"):
        select(policy=policy)


@pytest.mark.parametrize("key", ["rna_cds_boundary_ranking", "rna_cds_competition_ranking"])
def test_whole_model_policy_requires_competition_mode(key):
    policy = deepcopy(POLICIES["minimal"])
    del policy[key]
    with pytest.raises(SelectionError, match="requires RNA CDS competition"):
        select(policy=policy)


@pytest.mark.parametrize("damage", ["missing", "stale", "status", "chain"])
def test_missing_or_inconsistent_witnesses_fail_closed(damage):
    record = spliced("gene", protein=True, full=True)
    if damage == "missing":
        del record.row["rna_cds_support_transcripts_json"]
    elif damage == "stale":
        record.row["rna_cds_boundary_version"] = "1"
    elif damage == "status":
        record.row["rna_cds_status"] = "UNRESOLVED"
    else:
        record.row["rna_cds_chain_support_transcripts_json"] = "[]"
    with pytest.raises(SelectionError, match="recompute components"):
        select([record])


def test_new_experiments_do_not_change_admission_or_legacy_policy():
    old = json.loads((ROOT / "policies/unified_v5_single_exon50_competition.json").read_text())
    for mode in MODES:
        new = json.loads((ROOT / f"policies/unified_v6_single_exon50_{mode}.json").read_text())
        for key in ("policy_id", "policy_version", "note"):
            new[key] = old[key]
        new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
        assert new["selector_core"].pop("whole_model_ranking") == mode
        assert new == old
    records, _ = eight_records()
    legacy = select(records, policy=old["selector_core"])
    assert all(record.whole_model_ranking == "" for record in legacy.records.values())
