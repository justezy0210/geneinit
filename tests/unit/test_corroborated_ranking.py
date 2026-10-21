from copy import deepcopy
import itertools
import json

import pytest

from private_fixtures import load_private_json

from geneinit.evidence.components import parse_cds_signature
from geneinit.selection.core import Record, SelectionError
from geneinit.selection.unified import UnifiedSelector
from test_rna_cds_boundary import frozen_jinju_records, ranked_model
from test_rna_competition import eight_records
from test_unified_selection import make_record
from test_whole_model_ranking import ROOT, boundary, fusion_models, snapshot, spliced


POLICY_DOCUMENT = json.loads(
    (ROOT / "policies/unified_v7_single_exon50_corroborated.json").read_text()
)
POLICY = POLICY_DOCUMENT["selector_core"]
V5 = json.loads(
    (ROOT / "policies/unified_v5_single_exon50_competition.json").read_text()
)["selector_core"]


def select(records=(), policy=None):
    selector = UnifiedSelector("ASM", policy or POLICY,
                               sample_count=65, replicates_per_condition=1)
    selector.update(records)
    return selector


def route_rank(selector, key):
    return -json.loads(snapshot(selector)[key]["selection_rank_json"])[0]


def multihit_model():
    record = make_record(
        "multihit_725", ((100, 999, 0), (1100, 1999, 0), (2100, 2477, 0)),
        protein_aa=725, protein_tier=4, protein_precision=0.933,
        protein_intron_precision=1, protein_intron_recall=1, terminals=2,
        rna_tier=2, hq=[1, 1], samples=[1, 1], replicate=[0, 0],
    )
    record.row.update(
        protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
        protein_phase_fraction="1", protein_chain_complete="1",
    )
    return boundary(record, full=True)


def test_rna_and_multihit_protein_agreement_can_outweigh_unique_fragment():
    full = multihit_model()
    short = ranked_model("unique_137", 137)
    assert [r.key for r in select([full, short], V5).selected] == ["unique_137"]
    result = select([short, full])
    assert [r.key for r in result.selected] == ["multihit_725"]
    assert route_rank(result, full.key) == 9
    observed = result.records[full.key]
    assert observed.ranking_protein_fit
    assert not observed.strong_unique_protein
    assert "RNA_PROTEIN_CDS_CONCORDANCE" in observed.reasons
    assert "MULTI_HIT_FIT_USED_FOR_RANK_ONLY" in observed.reasons
    assert observed.routes == select([full], V5).records[full.key].routes


def test_single_exon_430_with_full_rna_beats_rna_only_70():
    complete = ranked_model("complete_430", 430, full_rna=True)
    fragment = boundary(make_record(
        "fragment_70", ((100, 198, 0), (300, 413, 0)), protein_aa=70,
    ), full=True)
    assert [r.key for r in select([complete, fragment], V5).selected] == ["fragment_70"]
    result = select([fragment, complete])
    assert [r.key for r in result.selected] == ["complete_430"]
    assert route_rank(result, complete.key) == 9
    assert not result.records[complete.key].strong_rna
    assert "SINGLE_EXON_RNA_PROTEIN_DUAL_RANK" in result.records[complete.key].reasons


@pytest.mark.parametrize("terminals", [(0, 0), (0, 1), (1, 0)])
def test_core_and_whole_rna_do_not_require_exact_ends_or_complete_chain(terminals):
    complete = spliced("core", protein=True, full=True)
    complete.row.update(
        protein_start_match=str(terminals[0]), protein_stop_match=str(terminals[1]),
        protein_chain_complete="0", protein_tier="3",
        rna_tier="2", rna_hq_read_counts="[1,1]",
        rna_sample_hq_ge1_counts="[2,2]", rna_condition_ge2rep_hq_ge1_counts="[0,0]",
    )
    boundary(complete, full=True)
    fragment = spliced("partial")
    fragment.row.update(
        rna_hq_read_counts="[8,0]", rna_sample_hq_ge1_counts="[2,0]",
        rna_condition_ge2rep_hq_ge1_counts="[1,0]",
    )
    boundary(fragment, conflict=True)
    assert [r.key for r in select([complete, fragment], V5).selected] == ["partial"]
    result = select([fragment, complete])
    assert [r.key for r in result.selected] == ["core"]
    assert route_rank(result, complete.key) == 9
    assert not result.records[complete.key].strong_unique_full_protein
    assert "RNA_PROTEIN_CDS_CONCORDANCE" in result.records[complete.key].reasons


def test_self_strong_chain_can_corroborate_multihit_without_full_cds_rna():
    record = spliced("chain", protein=True, chain=True)
    record.row.update(
        protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
        protein_start_match="0", protein_stop_match="0",
    )
    result = select([record])
    assert result.records[record.key].selected
    assert route_rank(result, record.key) == 9
    assert record.row["rna_cds_support_transcripts_json"] == "[]"
    assert "DUAL_RANK_STRONG_CODING_CHAIN_WITNESS" in result.records[record.key].reasons


@pytest.mark.parametrize("metric", [
    "protein_model_precision", "protein_chain_recall", "protein_cds_f1",
    "protein_phase_fraction", "protein_intron_precision",
])
@pytest.mark.parametrize("value", ["NA", "0.79"])
def test_partial_or_unphased_protein_does_not_gain_rna_protein_bonus(metric, value):
    record = spliced("partial", protein=True, full=True)
    record.row[metric] = value
    result = select([record])
    observed = result.records[record.key]
    assert observed.selected and observed.strong_rna
    assert not observed.ranking_protein_fit
    assert route_rank(result, record.key) == 8
    assert "RNA_PROTEIN_CDS_CONCORDANCE" not in observed.reasons


def test_tier_two_and_absent_whole_rna_do_not_gain_new_bonus():
    for protein_tier, full in [(2, True), (3, False)]:
        record = spliced("gene", protein=True, full=full)
        record.row.update(protein_tier=str(protein_tier), protein_stop_match="0")
        result = select([record])
        assert route_rank(result, record.key) == 8
        assert "RNA_PROTEIN_CDS_CONCORDANCE" not in result.records[record.key].reasons


def test_ranking_fit_cannot_admit_an_unqualified_single_exon_multihit():
    record = ranked_model("unqualified", 725, full_rna=True)
    record.row.update(
        protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
        protein_start_match="0", protein_stop_match="0",
    )
    result = select([record])
    assert result.records[record.key].ranking_protein_fit
    assert not result.records[record.key].eligible
    assert result.records[record.key].routes == []
    assert not result.selected


def test_legacy_dual_and_rna_first_tie_break_remain():
    legacy_dual = spliced("legacy", protein=True)
    result = select([legacy_dual])
    assert route_rank(result, legacy_dual.key) == 9
    assert "RNA50_ADMISSION_WITHOUT_WHOLE_RNA_WITNESS" in result.records[legacy_dual.key].reasons
    rna = spliced("rna")
    protein = ranked_model("protein", 430)
    for policy in (V5, POLICY):
        result = select([protein, rna], policy)
        assert route_rank(result, "rna") == route_rank(result, "protein") == 8
        assert [r.key for r in result.selected] == ["rna"]


def test_shorter_with_strong_rna_still_wins_when_both_have_corroboration():
    shorter = boundary(make_record(
        "shorter_704", ((100, 999, 0), (1100, 2314, 0)),
        protein_aa=704, protein_tier=4, protein_unique=True, protein_precision=1,
        protein_intron_precision=1, protein_intron_recall=1, terminals=2,
    ), full=True)
    longer = boundary(make_record(
        "longer_733", ((100, 999, 0), (1100, 2401, 0)),
        protein_aa=733, protein_tier=4, protein_unique=True, protein_precision=1,
        protein_intron_precision=1, protein_intron_recall=1, terminals=2,
        hq=[1], samples=[1], replicate=[0], rna_tier=2,
    ), full=True)
    result = select([longer, shorter])
    assert route_rank(result, shorter.key) == route_rank(result, longer.key) == 9
    assert [r.key for r in result.selected] == [shorter.key]


@pytest.mark.parametrize("support,precedence", [
    ({}, False), ({"full": True}, True), ({"chain": True}, True),
    ({"chain": True, "conflict": True}, False),
])
def test_fusion_requires_whole_rna_without_unopposed_conflict(support, precedence):
    records = fusion_models(**support)
    result = select(records)
    center = result.records["center"]
    assert center.active_fusion and center.eligible
    assert center.fusion_center_precedence is precedence
    assert center.selected is precedence


def test_focal_434_and_c06_legacy_choice_survive():
    records, pairs = eight_records()
    result = select([*records, *frozen_jinju_records()])
    assert result.records["GCAN-f25132743f098767e8652d6a"].selected
    c06 = next(pair for pair in pairs if pair["case"] == "C06")
    assert result.records[c06["lost_candidate"]].selected
    assert not result.records[c06["selected_blocker"]].selected


@pytest.mark.parametrize("pair_index", [0, 1])
def test_real_multihit_725_and_one_terminal_750_are_not_displaced(pair_index):
    fixture = load_private_json("jinju_corroborated_regressions.json")
    preferred, fragment = fixture["pairs"][pair_index]
    records = [
        Record(dict(row), parse_cds_signature(row["cds_signature"]))
        for row in fixture["models"] if row["candidate_id"] in {preferred, fragment}
    ]
    v6_policy = json.loads(
        (ROOT / "policies/unified_v6_single_exon50_minimal.json").read_text()
    )["selector_core"]
    assert [r.key for r in select(records, v6_policy).selected] == [fragment]
    result = select(records)
    assert [r.key for r in result.selected] == [preferred]
    assert route_rank(result, preferred) == 9
    assert "RNA_PROTEIN_CDS_CONCORDANCE" in result.records[preferred].reasons
    assert snapshot(select(list(reversed(records)))) == snapshot(result)
    restored = select()
    restored.restore(records, snapshot(result))
    assert snapshot(restored) == snapshot(result)
    assert restored.records[preferred].ranking_protein_fit


def test_checkpoint_recomputes_fit_and_updates_match_fresh_selection():
    records = [multihit_model(), ranked_model("fragment", 137)]
    original = select(records)
    expected = snapshot(original)
    for order in itertools.permutations(records):
        assert snapshot(select(order)) == expected
    restored = select()
    restored.restore(records, expected)
    assert snapshot(restored) == expected
    assert restored.records["multihit_725"].ranking_protein_fit
    assert restored._evidence["multihit_725"][1].ranking_protein_fit
    changed = deepcopy(records[0])
    changed.row["protein_phase_fraction"] = "NA"
    restored.update([changed])
    fresh = select([changed, records[1]])
    assert snapshot(restored) == snapshot(fresh)
    assert not restored.records[changed.key].ranking_protein_fit
    assert [r.key for r in restored.selected] == ["fragment"]
    restored.update([records[0]])
    assert snapshot(restored) == expected


@pytest.mark.parametrize("damage", ["missing", "stale"])
def test_corrupted_rna_evidence_cannot_supply_concordance(damage):
    record = multihit_model()
    if damage == "missing":
        del record.row["rna_cds_support_transcripts_json"]
    else:
        record.row["rna_cds_boundary_version"] = "1"
    with pytest.raises(SelectionError, match="recompute components"):
        select([record])


def test_policy_changes_only_ranking_and_experiment_metadata():
    old = json.loads(
        (ROOT / "policies/unified_v5_single_exon50_competition.json").read_text()
    )
    new = deepcopy(POLICY_DOCUMENT)
    for key in ("policy_id", "policy_version", "note"):
        new[key] = old[key]
    new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
    assert new["selector_core"].pop("whole_model_ranking") == "corroborated"
    assert new == old
