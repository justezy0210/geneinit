from copy import deepcopy
import itertools
import json

import pytest

from private_fixtures import load_private_json, private_pairs

from geneinit.evidence.components import parse_cds_signature
from geneinit.selection.core import Record
from geneinit.selection.unified import UnifiedSelector
from test_corroborated_ranking import POLICY as V7, route_rank
from test_directional_ranking import POLICY as V8, real_pair as directional_pair
from test_rna_cds_boundary import frozen_jinju_records, ranked_model
from test_rna_competition import eight_records
from test_whole_model_ranking import ROOT, fusion_models, snapshot, spliced


DOCUMENT = json.loads(
    (ROOT / "policies/unified_v9_single_exon50_corroborated_guarded.json").read_text()
)
POLICY = DOCUMENT["selector_core"]
FIXTURE_NAME = "jinju_guarded_regressions.json"


def select(records=(), policy=None):
    selector = UnifiedSelector("ASM", policy or POLICY,
                               sample_count=65, replicates_per_condition=1)
    selector.update(records)
    return selector


def records_for(*keys):
    return [
        Record(dict(row), parse_cds_signature(row["cds_signature"]))
        for row in load_private_json(FIXTURE_NAME)["models"] if row["candidate_id"] in keys
    ]


@pytest.mark.parametrize("pair", private_pairs(FIXTURE_NAME))
def test_real_pair_gate_restores_previous_model_without_claiming_complete_locus(pair):
    records = records_for(pair["old"], pair["new"])
    previous = select(records, V8)
    result = select(records)
    assert [record.key for record in previous.selected] == [pair["new"]]
    assert [record.key for record in result.selected] == [pair["old"]]
    candidate = result.records[pair["new"]]
    assert not candidate.introns and not candidate.ranking_protein_fit
    assert route_rank(previous, candidate.key) == 9
    assert route_rank(result, candidate.key) < 9
    assert "SINGLE_EXON_CONCORDANCE_TERMINAL_GUARD" in candidate.reasons
    for key in result.records:
        assert result.records[key].eligible == previous.records[key].eligible
        assert result.records[key].routes == previous.records[key].routes
    # The available audit subset omits an unchanged overlapping 946-aa model
    # at the 1054-aa locus; this pair test is not full-locus restoration proof.


def test_actual_r43_documents_loss_of_well_supported_one_terminal_single_exon():
    pair = load_private_json(FIXTURE_NAME)["tradeoff"]
    records = records_for(pair["old"], pair["new"])
    previous = select(records, V8)
    result = select(records)
    assert [record.key for record in previous.selected] == [pair["old"]]
    assert [record.key for record in result.selected] == [pair["new"]]
    old = result.records[pair["old"]]
    assert old.row["translated_protein_aa"] == "364"
    assert not old.introns and old.row["rna_cds_status"] == "SUPPORTED"
    assert old.eligible
    assert route_rank(previous, old.key) == 9
    assert route_rank(result, old.key) < 9


@pytest.mark.parametrize("terminals", [(0, 0), (1, 0), (0, 1), (1, 1)])
@pytest.mark.parametrize("unique", [False, True])
def test_single_exon_bonus_requires_both_ends_without_adding_uniqueness(terminals, unique):
    record = ranked_model("single", 50, full_rna=True)
    record.row.update(
        protein_start_match=str(terminals[0]), protein_stop_match=str(terminals[1]),
        protein_multiplicity_status=(
            "UNIQUE_PROTEIN_HIT" if unique else "AMBIGUOUS_PROTEIN_MULTI_HIT"
        ),
    )
    result = select([record])
    observed = result.records[record.key]
    assert observed.hard_qc  # Existing SE50 floor.
    assert observed.ranking_protein_fit is (sum(terminals) == 2)
    assert (route_rank(result, record.key) == 9) is (sum(terminals) == 2)
    if sum(terminals):
        assert observed.eligible
    if sum(terminals) == 2:
        assert "SINGLE_EXON_RNA_PROTEIN_DUAL_RANK" in observed.reasons
        if not unique:
            assert "MULTI_HIT_FIT_USED_FOR_RANK_ONLY" in observed.reasons


@pytest.mark.parametrize("terminals", [(0, 0), (1, 0), (0, 1)])
def test_spliced_core_concordance_keeps_one_or_zero_terminal_models(terminals):
    record = spliced("spliced", protein=True, full=True)
    record.row.update(
        protein_start_match=str(terminals[0]), protein_stop_match=str(terminals[1]),
        protein_chain_complete="0",
    )
    result = select([record])
    assert route_rank(result, record.key) == 9
    assert result.records[record.key].ranking_protein_fit
    assert "SINGLE_EXON_CONCORDANCE_TERMINAL_GUARD" not in result.records[record.key].reasons


def test_terminal_guard_uses_existing_policy_requirement():
    record = ranked_model("single", 66, full_rna=True)
    record.row["protein_stop_match"] = "0"
    policy = deepcopy(POLICY)
    policy["protein"]["required_terminal_matches"] = 1
    result = select([record], policy)
    assert result.records[record.key].ranking_protein_fit
    assert route_rank(result, record.key) == 9


def test_directional_hypothesis_remains_only_in_v8_and_v7_examples_survive():
    records, (supported, fragment) = directional_pair()
    result = select(records)
    assert [record.key for record in select(records, V8).selected] == [supported]
    assert [record.key for record in select(records, V7).selected] == [fragment]
    assert [record.key for record in result.selected] == [fragment]
    assert not result.records[supported].ranking_directional_protein_fit
    assert "RNA_PROTEIN_DIRECTIONAL_CONCORDANCE" not in result.records[supported].reasons
    fixture = load_private_json("jinju_corroborated_regressions.json")
    records = [
        Record(dict(row), parse_cds_signature(row["cds_signature"]))
        for row in fixture["models"]
    ]
    result = select(records)
    for preferred, fragment in fixture["pairs"]:
        assert result.records[preferred].selected
        assert not result.records[fragment].selected
    complete = ranked_model("complete_430", 430, full_rna=True)
    result = select([complete, spliced("rna_fragment", full=True)])
    assert [record.key for record in result.selected] == [complete.key]


def test_legacy_dual_focal_434_and_c06_remain():
    legacy = spliced("legacy", protein=True)
    result = select([legacy])
    assert route_rank(result, legacy.key) == 9
    records, pairs = eight_records()
    result = select([*records, *frozen_jinju_records()])
    assert result.records["GCAN-f25132743f098767e8652d6a"].selected
    c06 = next(pair for pair in pairs if pair["case"] == "C06")
    assert result.records[c06["lost_candidate"]].selected
    assert not result.records[c06["selected_blocker"]].selected


@pytest.mark.parametrize("support,precedence", [
    ({}, False), ({"full": True}, True), ({"chain": True}, True),
    ({"chain": True, "conflict": True}, False),
])
def test_guarded_mode_preserves_fusion_gate(support, precedence):
    result = select(fusion_models(**support))
    center = result.records["center"]
    assert center.active_fusion and center.eligible
    assert center.fusion_center_precedence is precedence
    assert center.selected is precedence


def test_guarded_checkpoint_update_and_permutation_match_fresh_selection():
    pair = load_private_json(FIXTURE_NAME)["pairs"][0]
    records = records_for(pair["old"], pair["new"])
    result = select(records)
    expected = snapshot(result)
    for order in itertools.permutations(records):
        assert snapshot(select(order)) == expected
    restored = select()
    restored.restore(records, expected)
    assert snapshot(restored) == expected
    assert not restored.records[pair["new"]].ranking_protein_fit
    changed = deepcopy(next(record for record in records if record.key == pair["new"]))
    changed.row.update(protein_start_match="1", protein_stop_match="1")
    restored.update([changed])
    assert restored.records[pair["new"]].ranking_protein_fit
    fresh = select([
        changed, *(record for record in records if record.key != pair["new"]),
    ])
    assert snapshot(restored) == snapshot(fresh)
    assert restored.records[pair["new"]].selected
    restored.update([next(record for record in records if record.key == pair["new"])])
    assert snapshot(restored) == expected


def test_policy_changes_only_metadata_and_mode():
    old = json.loads(
        (ROOT / "policies/unified_v7_single_exon50_corroborated.json").read_text()
    )
    new = deepcopy(DOCUMENT)
    for key in ("policy_id", "policy_version", "note"):
        new[key] = old[key]
    new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
    new["selector_core"]["whole_model_ranking"] = old["selector_core"]["whole_model_ranking"]
    assert new == old
