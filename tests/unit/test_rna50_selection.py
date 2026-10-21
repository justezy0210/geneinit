"""The default RNA fraction admits complete models without changing their CDS."""
import copy
import json
from pathlib import Path
import runpy

import pytest

from geneinit.selection.core import audit_row, build_loci_and_edges, evaluate_and_select


ROOT = Path(__file__).resolve().parents[2]
POLICY = json.loads((ROOT / "policies/production_v2.json").read_text())["selector_core"]
helpers = runpy.run_path(str(Path(__file__).with_name("test_selection_core.py")))


def record(supported, total=10, *, hq=3, samples=2, replicate=1):
    result = helpers["make_record"](
        "RNA50", [(100 + i * 200, 198 + i * 200, 0) for i in range(total + 1)],
        hq=[hq] * supported + [0] * (total - supported),
        samples=[samples] * supported + [0] * (total - supported),
        replicate=[replicate] * supported + [0] * (total - supported),
    )
    result.row.update(candidate_pool="GENERAL", repeat_cds_fraction="1")
    return result


def select(result, policy=POLICY):
    loci, edges = build_loci_and_edges("TEST", [result])
    return evaluate_and_select(loci, edges, policy, sample_count=6, replicates_per_condition=2)


@pytest.mark.parametrize("supported,total,expected", [(4, 10, False), (5, 10, True), (7, 9, True)])
def test_fraction_boundary_and_multiple_unsupported_connections(supported, total, expected):
    result = record(supported, total)
    before = (result.cds, result.row["cds_signature"])
    assert bool(select(result)) is expected
    assert (result.cds, result.row["cds_signature"]) == before
    if expected:
        assert "RECURRENT_RNA_50" in result.routes
        assert result.repeat_override_reason == "OVERLAP_REPORT_ONLY"


@pytest.mark.parametrize("hq,samples,replicate,expected", [(1, 1, 0, False), (3, 1, 0, False), (2, 2, 0, True)])
def test_existing_recurrent_connection_quality(hq, samples, replicate, expected):
    result = record(5, hq=hq, samples=samples, replicate=replicate)
    result.row["rna_tier"] = "2"
    assert bool(select(result)) is expected


@pytest.mark.parametrize("field,value", [("orf_status", "missing_stop"), ("common_qc_pass", "0"), ("rna_tier", "1")])
def test_structural_qc_and_rna_tier_remain_required(field, value):
    result = record(5)
    result.row[field] = value
    assert select(result) == []


def test_single_exon_requires_existing_non_rna_route():
    result = record(0, 0)
    assert select(result) == []
    assert audit_row(result, POLICY["policy_id"])["rna_recurrent_intron_fraction"] == "NA"


def test_confirmed_te_stays_separate_at_half_rna_support():
    result = record(5)
    result.row["candidate_pool"] = "TE"
    assert select(result) == []
    assert result.row["te_selected"] == "1"
    assert result.decision == "TE_SELECTED"


def test_audit_records_actual_fraction_and_missing_coordinates():
    result = record(7, 9)
    assert select(result)
    audit = audit_row(result, POLICY["policy_id"])
    assert audit["rna_recurrent_intron_count"] == 7
    assert float(audit["rna_recurrent_intron_fraction"]) == pytest.approx(7 / 9)
    assert json.loads(audit["rna_nonrecurrent_introns_json"]) == [list(x) for x in result.introns[7:]]


def test_explicit_legacy_cap_and_route_remain_effective():
    legacy = copy.deepcopy(POLICY)
    legacy["rna"].update(min_recurrent_fraction=0.8, max_nonrecurrent_introns=1)
    assert select(record(8, 10), legacy) == []
    result = record(9, 10)
    assert select(result, legacy)
    assert "RECURRENT_RNA_80" in result.routes
