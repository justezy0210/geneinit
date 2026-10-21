from copy import deepcopy
import json
from pathlib import Path

import pytest

from geneinit.selection.unified import UnifiedSelector
from test_unified_selection import make_record


ROOT = Path(__file__).resolve().parents[2]
POLICY = json.loads((ROOT / "policies/unified_v3_single_exon50.json").read_text())


def select(records):
    selector = UnifiedSelector(
        "ASM", POLICY["selector_core"], sample_count=6, replicates_per_condition=2,
    )
    selector.update(records)
    return selector


def supported(length, key="small"):
    return make_record(
        key, ((100, 100 + 3 * (length + 1) - 1, 0),), protein_aa=length,
        protein_tier=4, protein_unique=True, protein_precision=1, terminals=2,
    )


def test_only_single_exon_floor_changes_besides_policy_identity_and_note():
    old = json.loads((ROOT / "policies/unified_v3.json").read_text())
    new = deepcopy(POLICY)
    for key in ("policy_id", "note"):
        new[key] = old[key]
    new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
    assert new["selector_core"]["hard_qc"]["min_protein_aa"] == 50
    assert new["selector_core"]["hard_qc"]["single_exon_min_protein_aa"] == 50
    new["selector_core"]["hard_qc"]["single_exon_min_protein_aa"] = 100
    assert new == old


@pytest.mark.parametrize("length,accepted", [(49, False), (50, True), (72, True),
                                           (99, True), (100, True)])
def test_supported_short_models_obey_the_50aa_boundary(length, accepted):
    selector = select([supported(length)])
    assert bool(selector.selected) is accepted


def test_short_models_still_need_independent_evidence():
    record = make_record("unsupported72", ((100, 318, 0),), protein_aa=72)
    selector = select([record])
    assert selector.records[record.key].hard_qc
    assert not selector.selected
    assert selector.records[record.key].reasons == ["NO_INDEPENDENT_EVIDENCE_ROUTE"]


def test_short_models_still_need_a_complete_orf():
    record = supported(72)
    record.row["orf_status"] = "missing_start"
    selector = select([record])
    assert not selector.selected
    assert selector.records[record.key].reasons == ["HARD_QC_FAIL"]


def test_short_models_still_compete_with_longer_supported_models():
    short, long = supported(72), supported(120, "long")
    selector = select([short, long])
    assert [record.key for record in selector.selected] == ["long"]
    assert selector.records[short.key].decision == "REJECT_INCOMPATIBLE_GRAPH"
