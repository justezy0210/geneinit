from copy import deepcopy
from dataclasses import replace
import itertools
import json

import pytest

from private_fixtures import load_private_json, private_pairs

from geneinit.evidence.components import (
    CandidateRecord, RnaTranscript, format_value, load_rna_cds_transcripts,
    parse_cds_signature,
)
from geneinit.evidence.rna_cds import (
    RNA_CDS_COLUMNS, RNA_CDS_SAMPLE_COLUMNS, full_cds_sample_metrics, rna_cds_metrics,
)
from geneinit.models.candidate import CdsPart
from geneinit.postqc.core import load_raw_gff3, run_postqc
from geneinit.postqc.unified import RepairEvidence
from geneinit.selection.core import Record, SelectionError
from geneinit.selection.unified import UnifiedSelector
from test_corroborated_ranking import POLICY as V7, route_rank
from test_directional_ranking import real_pair as directional_pair
from test_guarded_ranking import POLICY as V9
from test_rna_cds_boundary import frozen_jinju_records, ranked_model
from test_rna_competition import eight_records
from test_unified_postqc import (
    make_inputs, test_workflow_selection_and_uniform_repair_reach_a_fixed_point as run_fixture,
)
from test_whole_model_ranking import ROOT, fusion_models, snapshot, spliced


POLICY_NAME = "unified_v10_single_exon50_corroborated_multisample.json"
DOCUMENT = json.loads((ROOT / "policies" / POLICY_NAME).read_text())
POLICY = DOCUMENT["selector_core"]
FIXTURE_NAME = "jinju_multisample_regressions.json"


def with_sources(record, sources=()):
    record = deepcopy(record)
    record.row.update({
        key: format_value(value)
        for key, value in full_cds_sample_metrics(sources).items()
    })
    return record


def select(records=(), policy=None):
    selector = UnifiedSelector("ASM", policy or POLICY,
                               sample_count=65, replicates_per_condition=1)
    selector.update(records)
    return selector


def write_gtf(path, transcripts):
    path.write_text("".join(
        f'chr1\ttest\texon\t{a}\t{b}\t.\t+\t.\ttranscript_id "{name}";\n'
        for name, exons in transcripts for a, b in exons
    ))
    return path


def candidate():
    return CandidateRecord(
        "candidate", "structure", "chr1", "+", (CdsPart(100, 399, 0),),
        (), 300, "complete",
    )


@pytest.mark.parametrize("pair", private_pairs(FIXTURE_NAME))
def test_real_pairs_keep_multisample_single_exons_and_reject_sparse_new_concordance(pair):
    records = [
        Record(dict(row), parse_cds_signature(row["cds_signature"]))
        for row in load_private_json(FIXTURE_NAME)["models"] if row["candidate_id"] in {pair["old"], pair["new"]}
    ]
    v7, v9, result = select(records, V7), select(records, V9), select(records)
    assert [record.key for record in result.selected] == [pair["old"]]
    if pair["case"].startswith("R"):
        assert [record.key for record in v7.selected] == [pair["old"]]
        assert [record.key for record in v9.selected] == [pair["new"]]
        supported = result.records[pair["old"]]
        assert not supported.introns
        assert supported.ranking_multisample_rna
        assert "MULTIPLE_SAMPLE_GTFS_SUPPORT_FULL_CDS" in supported.reasons
    else:
        assert [record.key for record in v7.selected] == [pair["new"]]
        assert [record.key for record in v9.selected] == [pair["old"]]
        sparse = result.records[pair["new"]]
        assert not sparse.introns and not sparse.ranking_multisample_rna
        assert not sparse.ranking_protein_fit
    for key in result.records:
        assert (result.records[key].eligible, result.records[key].routes) == (
            v7.records[key].eligible, v7.records[key].routes,
        )
    # These are two-model checks, not a claim about every alternative at the locus.


@pytest.mark.parametrize("include_merged", [False, True])
def test_identical_chains_retain_two_file_sources_and_original_representative(tmp_path, include_merged):
    exons = ((90, 410),)
    merged = [RnaTranscript("merged", "chr1", "+", exons)] if include_merged else []
    a = write_gtf(tmp_path / "a.gtf", [("a", exons)])
    b = write_gtf(tmp_path / "b.gtf", [("b", exons)])
    alias = tmp_path / "alias.gtf"
    alias.symlink_to(a)
    first = load_rna_cds_transcripts(merged, [b, a, alias, a])
    second = load_rna_cds_transcripts(merged, [a, b])
    assert first == second and len(first) == 1
    assert first[0].transcript_id == ("merged" if include_merged else f"{a.resolve()}::a")
    assert first[0].sample_sources == (str(a.resolve()), str(b.resolve()))
    metrics = rna_cds_metrics(candidate(), first, {})
    assert metrics["rna_cds_support_sample_count"] == 2
    # Adding provenance must not change any existing compatibility observation.
    old_shape = [replace(tx, sample_sources=()) for tx in first]
    old_metrics = rna_cds_metrics(candidate(), old_shape, {})
    assert {key: metrics[key] for key in RNA_CDS_COLUMNS} == {
        key: old_metrics[key] for key in RNA_CDS_COLUMNS
    }


def test_one_file_many_full_isoforms_counts_once_and_partial_file_does_not_count(tmp_path):
    sample = write_gtf(tmp_path / "sample.gtf", [
        ("full1", ((90, 410),)), ("full2", ((80, 420),)),
        ("full_duplicate", ((90, 410),)),
    ])
    partial = write_gtf(tmp_path / "partial.gtf", [("partial", ((200, 410),))])
    transcripts = load_rna_cds_transcripts([], [sample, partial])
    metrics = rna_cds_metrics(candidate(), transcripts, {})
    assert len(json.loads(metrics["rna_cds_support_transcripts_json"])) == 2
    assert metrics["rna_cds_support_sample_count"] == 1
    assert json.loads(metrics["rna_cds_support_sample_sources_json"]) == [str(sample.resolve())]


def test_merged_only_support_is_not_a_sample_file():
    merged = [RnaTranscript("merged", "chr1", "+", ((90, 410),))]
    metrics = rna_cds_metrics(candidate(), load_rna_cds_transcripts(merged), {})
    assert metrics["rna_cds_status"] == "SUPPORTED"
    assert metrics["rna_cds_support_sample_count"] == 0
    assert metrics["rna_cds_support_sample_sources_json"] == "[]"


@pytest.mark.parametrize("sample_files", [0, 1, 2, 3])
@pytest.mark.parametrize("terminals", [(0, 0), (1, 0), (0, 1), (1, 1)])
def test_single_exon_concordance_uses_terminals_or_distinct_sample_files(sample_files, terminals):
    record = ranked_model("single", 100, full_rna=True)
    record.row.update(protein_start_match=str(terminals[0]), protein_stop_match=str(terminals[1]))
    record = with_sources(record, [f"/samples/sample{i}.gtf" for i in range(sample_files)])
    result = select([record])
    observed = result.records[record.key]
    expected = sum(terminals) == 2 or sample_files >= 2
    assert observed.ranking_protein_fit is expected
    assert (route_rank(result, record.key) == 9) is expected
    old = select([record], V7).records[record.key]
    assert (observed.eligible, observed.routes) == (old.eligible, old.routes)
    assert not observed.ranking_directional_protein_fit
    if sample_files >= 2 and observed.eligible:
        assert "MULTIPLE_SAMPLE_GTFS_SUPPORT_FULL_CDS" in observed.reasons


def test_policy_minimum_is_explicit_and_not_hardcoded():
    record = with_sources(ranked_model("single", 100, full_rna=True),
                          ["/samples/a.gtf", "/samples/b.gtf"])
    record.row["protein_stop_match"] = "0"
    policy = deepcopy(POLICY)
    policy["single_exon_min_full_cds_sample_sources"] = 3
    assert not select([record], policy).records[record.key].ranking_protein_fit
    for value in (None, True, 1, "2"):
        policy["single_exon_min_full_cds_sample_sources"] = value
        with pytest.raises(SelectionError, match="explicit sample-GTF"):
            select([], policy)


@pytest.mark.parametrize("missing", RNA_CDS_SAMPLE_COLUMNS)
def test_missing_provenance_fails_closed_even_for_two_terminal_or_spliced_candidate(missing):
    for record in (ranked_model("single", 100, full_rna=True), spliced("spliced")):
        record = with_sources(record)
        del record.row[missing]
        with pytest.raises(SelectionError, match="recompute"):
            select([record])
        assert select([record], V9).records[record.key].eligible


@pytest.mark.parametrize("changes", [
    {"rna_cds_sample_provenance_version": "0"},
    {"rna_cds_support_sample_count": "3"},
    {"rna_cds_support_sample_sources_json": '["/a.gtf","/a.gtf"]'},
    {"rna_cds_support_sample_sources_json": '["relative.gtf"]', "rna_cds_support_sample_count": "1"},
    {"rna_cds_support_sample_sources_json": '["/a/../b.gtf"]', "rna_cds_support_sample_count": "1"},
    {"rna_cds_support_sample_sources_json": "{}"},
])
def test_stale_inconsistent_or_noncanonical_provenance_is_an_error(changes):
    record = with_sources(ranked_model("single", 100, full_rna=True), ["/a.gtf", "/b.gtf"])
    record.row.update(changes)
    with pytest.raises(SelectionError, match="recompute"):
        select([record])


def test_samples_cannot_claim_full_support_when_no_compatible_witness_exists():
    record = with_sources(ranked_model("single", 100), ["/a.gtf", "/b.gtf"])
    with pytest.raises(SelectionError, match="recompute"):
        select([record])


def test_multisample_support_does_not_rescue_partial_or_unphased_protein_or_bad_admission():
    for metric in ("protein_model_precision", "protein_chain_recall",
                   "protein_cds_f1", "protein_phase_fraction"):
        record = with_sources(ranked_model("single", 100, full_rna=True), ["/a.gtf", "/b.gtf"])
        record.row.update(protein_stop_match="0", **{metric: "NA"})
        assert not select([record]).records[record.key].ranking_protein_fit
    record = with_sources(ranked_model("no_admission", 100, full_rna=True), ["/a.gtf", "/b.gtf"])
    record.row.update(protein_multiplicity_status="AMBIGUOUS_PROTEIN_MULTI_HIT",
                      protein_start_match="0", protein_stop_match="0")
    result = select([record])
    assert result.records[record.key].ranking_protein_fit
    assert not result.records[record.key].eligible and not result.selected


def test_spliced_concordance_legacy_focal_and_directional_exclusion_remain():
    record = with_sources(spliced("spliced", protein=True, full=True))
    record.row.update(protein_start_match="0", protein_stop_match="0", protein_chain_complete="0")
    assert route_rank(select([record]), record.key) == 9
    records, pairs = eight_records()
    result = select([with_sources(record) for record in [*records, *frozen_jinju_records()]])
    assert result.records["GCAN-f25132743f098767e8652d6a"].selected
    c06 = next(pair for pair in pairs if pair["case"] == "C06")
    assert result.records[c06["lost_candidate"]].selected
    records, (supported, fragment) = directional_pair()
    result = select([with_sources(record) for record in records])
    assert [record.key for record in result.selected] == [fragment]
    assert not result.records[supported].ranking_directional_protein_fit


@pytest.mark.parametrize("support,precedence", [
    ({}, False), ({"full": True}, True), ({"chain": True}, True),
    ({"chain": True, "conflict": True}, False),
])
def test_multisample_mode_preserves_existing_fusion_gate(support, precedence):
    result = select([with_sources(record) for record in fusion_models(**support)])
    assert result.records["center"].fusion_center_precedence is precedence
    assert result.records["center"].selected is precedence


def test_provenance_update_restore_and_input_permutation_are_deterministic():
    single = with_sources(ranked_model("single", 100, full_rna=True), ["/a.gtf"])
    single.row["protein_stop_match"] = "0"
    rival = with_sources(spliced("rival"))
    result = select([single, rival])
    expected = snapshot(result)
    for records in itertools.permutations([single, rival]):
        assert snapshot(select(records)) == expected
    restored = select()
    restored.restore([single, rival], expected)
    assert snapshot(restored) == expected
    strengthened = with_sources(single, ["/a.gtf", "/b.gtf"])
    restored.update([strengthened])
    assert restored.records[single.key].ranking_multisample_rna
    assert snapshot(restored) == snapshot(select([strengthened, rival]))
    restored.update([single])
    assert snapshot(restored) == expected
    missing = deepcopy(single)
    del missing.row["rna_cds_sample_provenance_version"]
    with pytest.raises(SelectionError, match="recompute"):
        select().restore([missing, rival], expected)


def test_v10_actual_selection_postqc_reentry_and_structural_release(tmp_path):
    run_fixture(tmp_path, POLICY_NAME)


def test_runtime_synthetic_candidate_scores_the_same_full_cds_sample_provenance(tmp_path):
    args, genome, junctions = make_inputs(tmp_path, POLICY_NAME)
    run_postqc(*args)
    _, models = load_raw_gff3(tmp_path / "postqc.gff3")
    model = next(iter(models.values()))
    exons = tuple((f.start, f.end) for f in sorted(model.cds, key=lambda f: f.start))
    a = write_gtf(tmp_path / "a.gtf", [("full", exons)])
    b = write_gtf(tmp_path / "b.gtf", [("full", exons)])
    merged = write_gtf(tmp_path / "merged.gtf", [("merged_full", exons)])
    paths = dict(args[-1]["paths"], sample_gtfs=[str(a), str(b)], merged_gtf=str(merged))
    record = RepairEvidence.from_paths(
        paths, genome=genome, junctions=junctions, expected_protein_copies=2,
    ).score(model)
    assert record.key.startswith("GRNA-")
    assert record.row["rna_cds_support_sample_count"] == "2"
    assert json.loads(record.row["rna_cds_support_sample_sources_json"]) == [str(a), str(b)]
    assert record.row["rna_cds_support_transcripts_json"] == '["merged_full"]'
    assert select([record]).records[record.key].ranking_multisample_rna
