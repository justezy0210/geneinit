from copy import deepcopy
import itertools
import json

import pytest

from private_fixtures import load_private_json

from geneinit.evidence.components import parse_cds_signature
from geneinit.selection.core import Record, SelectionError
from geneinit.selection.unified import UnifiedSelector
from test_corroborated_ranking import POLICY as V7, route_rank
from test_rna_cds_boundary import frozen_jinju_records, ranked_model
from test_rna_competition import eight_records
from test_whole_model_ranking import ROOT, boundary, fusion_models, snapshot, spliced


DOCUMENT = json.loads(
    (ROOT / "policies/unified_v8_single_exon50_corroborated_directional.json").read_text()
)
POLICY = DOCUMENT["selector_core"]


def select(records=(), policy=None):
    selector = UnifiedSelector("ASM", policy or POLICY,
                               sample_count=65, replicates_per_condition=1)
    selector.update(records)
    return selector


def directional_model(key="directional", *, full=True, chain=False):
    record = spliced(key, protein=True)
    record.row.update(
        protein_tier="2", protein_model_precision="0.87", protein_phase_fraction="1",
        protein_chain_recall="0.53", protein_cds_f1="0.66",
        protein_intron_precision="1",
    )
    return boundary(record, full=full, chain=chain)


def real_pair():
    fixture = load_private_json("jinju_directional_regression.json")
    return [
        Record(dict(row), parse_cds_signature(row["cds_signature"]))
        for row in fixture["models"]
    ], fixture["pair"]


def test_real_341_candidate_can_outweigh_66_after_directional_corroboration():
    records, (supported, fragment) = real_pair()
    assert [r.key for r in select(records, V7).selected] == [fragment]
    result = select(records)
    assert [r.key for r in result.selected] == [supported]
    assert route_rank(result, supported) == route_rank(result, fragment) == 9
    observed = result.records[supported]
    assert observed.strong_rna and observed.ranking_directional_protein_fit
    assert not observed.ranking_protein_fit
    assert "RNA_PROTEIN_DIRECTIONAL_CONCORDANCE" in observed.reasons
    assert "CANDIDATE_PRECISION_RANK_ONLY_REFERENCE_RECALL_NOT_REQUIRED" in observed.reasons
    previous = select(records, V7)
    for key in result.records:
        assert result.records[key].routes == previous.records[key].routes
        assert result.records[key].eligible == previous.records[key].eligible


def test_directional_bonus_requires_full_cds_rna_and_all_own_strong_junctions():
    full = directional_model()
    assert route_rank(select([full]), full.key) == 9
    chain_only = directional_model(full=False, chain=True)
    assert chain_only.row["rna_cds_self_strong_chain"] == "1"
    assert route_rank(select([chain_only]), chain_only.key) == 8
    partial = directional_model()
    partial.row.update(
        rna_hq_read_counts="[8,0]", rna_sample_hq_ge1_counts="[2,0]",
        rna_condition_ge2rep_hq_ge1_counts="[1,0]",
    )
    boundary(partial, full=True)
    result = select([partial])
    assert result.records[partial.key].strong_rna  # RNA50 admission still passes.
    assert partial.row["rna_cds_self_strong_chain"] == "0"
    assert route_rank(result, partial.key) == 8
    assert "RNA_PROTEIN_DIRECTIONAL_CONCORDANCE" not in result.records[partial.key].reasons


def test_directional_bonus_still_requires_strong_rna_tier():
    record = directional_model()
    record.row["rna_tier"] = "2"
    result = select([record])
    assert result.records[record.key].ranking_directional_protein_fit
    assert record.row["rna_cds_self_strong_chain"] == "1"
    assert not result.records[record.key].strong_rna
    assert route_rank(result, record.key) < 9


@pytest.mark.parametrize("metric", [
    "protein_model_precision", "protein_phase_fraction", "protein_intron_precision",
])
@pytest.mark.parametrize("value", ["NA", "0.79", "nan"])
def test_directional_candidate_fit_fails_on_missing_or_weak_required_metrics(metric, value):
    record = directional_model()
    record.row[metric] = value
    result = select([record])
    assert result.records[record.key].selected
    assert not result.records[record.key].ranking_directional_protein_fit
    assert route_rank(result, record.key) == 8


def test_weak_protein_tier_cannot_gain_directional_bonus():
    record = directional_model()
    record.row["protein_tier"] = "1"
    result = select([record])
    assert not result.records[record.key].ranking_directional_protein_fit
    assert route_rank(result, record.key) == 8


def test_directional_fit_does_not_require_unique_hit_or_exact_ends():
    record = directional_model()
    record.row.update(
        protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
        protein_start_match="0", protein_stop_match="0", protein_chain_complete="0",
    )
    result = select([record])
    assert route_rank(result, record.key) == 9
    assert not result.records[record.key].strong_unique_protein
    assert "MULTI_HIT_FIT_USED_FOR_RANK_ONLY" in result.records[record.key].reasons


def test_directional_extension_does_not_promote_single_exon_partial_reference():
    record = ranked_model("single", 66, full_rna=True)
    record.row.update(
        protein_tier="2", protein_model_precision="0.87", protein_phase_fraction="1",
        protein_chain_recall="0.70", protein_cds_f1="0.75",
    )
    result = select([record])
    assert result.records[record.key].eligible  # Existing supplemental admission.
    assert record.row["rna_cds_self_strong_chain"] == "0"
    assert route_rank(result, record.key) < 9
    assert "RNA_PROTEIN_DIRECTIONAL_CONCORDANCE" not in result.records[record.key].reasons


def test_directional_fit_never_bypasses_hard_qc():
    record = directional_model()
    record.row["orf_status"] = "missing_start"
    result = select([record])
    assert result.records[record.key].ranking_directional_protein_fit
    assert not result.records[record.key].eligible
    assert result.records[record.key].routes == []
    assert not result.selected


def test_shorter_model_with_v7_protein_rna_agreement_can_still_win():
    records, (longer_id, shorter_id) = real_pair()
    longer = next(record for record in records if record.key == longer_id)
    for column in (
        "rna_hq_read_counts", "rna_sample_hq_ge1_counts",
        "rna_condition_ge2rep_hq_ge1_counts",
    ):
        counts = json.loads(longer.row[column])
        counts[0] = 0
        longer.row[column] = json.dumps(counts)
    boundary(longer, full=True)
    result = select(records)
    assert [r.key for r in result.selected] == [shorter_id]
    assert route_rank(result, longer_id) == 8
    assert route_rank(result, shorter_id) == 9


@pytest.mark.parametrize("support,precedence", [
    ({}, False), ({"full": True}, True), ({"chain": True}, True),
    ({"chain": True, "conflict": True}, False),
])
def test_directional_mode_preserves_fusion_gate(support, precedence):
    result = select(fusion_models(**support))
    center = result.records["center"]
    assert center.active_fusion and center.eligible
    assert center.fusion_center_precedence is precedence
    assert center.selected is precedence


def test_v7_multihit_and_one_terminal_examples_and_focal_choices_survive():
    fixture = load_private_json("jinju_corroborated_regressions.json")
    records = [
        Record(dict(row), parse_cds_signature(row["cds_signature"]))
        for row in fixture["models"]
    ]
    result = select(records)
    for preferred, fragment in fixture["pairs"]:
        assert result.records[preferred].selected
        assert not result.records[fragment].selected
    records, pairs = eight_records()
    result = select([*records, *frozen_jinju_records()])
    assert result.records["GCAN-f25132743f098767e8652d6a"].selected
    c06 = next(pair for pair in pairs if pair["case"] == "C06")
    assert result.records[c06["lost_candidate"]].selected
    assert not result.records[c06["selected_blocker"]].selected


def test_directional_checkpoint_update_and_permutation_recompute_flags():
    records, (supported, fragment) = real_pair()
    result = select(records)
    expected = snapshot(result)
    for order in itertools.permutations(records):
        assert snapshot(select(order)) == expected
    restored = select()
    restored.restore(records, expected)
    assert snapshot(restored) == expected
    assert restored.records[supported].ranking_directional_protein_fit
    assert restored._evidence[supported][1].ranking_directional_protein_fit
    changed = deepcopy(next(record for record in records if record.key == supported))
    changed.row["protein_phase_fraction"] = "NA"
    restored.update([changed])
    assert not restored.records[supported].ranking_directional_protein_fit
    assert snapshot(restored) == snapshot(select([
        changed, *(record for record in records if record.key != supported),
    ]))
    assert [record.key for record in restored.selected] == [fragment]


def test_directional_mode_still_requires_fresh_rna_evidence():
    record = directional_model()
    record.row["rna_cds_boundary_version"] = "1"
    with pytest.raises(SelectionError, match="recompute components"):
        select([record])


def test_v8_policy_changes_only_experiment_metadata_and_ranking_mode():
    old = json.loads(
        (ROOT / "policies/unified_v7_single_exon50_corroborated.json").read_text()
    )
    new = deepcopy(DOCUMENT)
    for key in ("policy_id", "policy_version", "note"):
        new[key] = old[key]
    new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
    new["selector_core"]["whole_model_ranking"] = old["selector_core"]["whole_model_ranking"]
    assert new == old
