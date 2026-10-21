import json
from pathlib import Path

import pytest

from private_fixtures import load_private_json

from geneinit.evidence.components import (
    CandidateRecord, JunctionAggregate, RnaTranscript, best_rna_metrics,
    format_value, load_rna_cds_transcripts, parse_cds_signature,
)
from geneinit.evidence.metrics import IntervalIndex
from geneinit.evidence.rna_cds import RNA_CDS_COLUMNS, RNA_CDS_SAMPLE_COLUMNS, rna_cds_metrics
from geneinit.models.candidate import CdsPart
from geneinit.selection.core import Record, SelectionError, audit_row
from geneinit.selection.unified import UnifiedSelector
from test_rna_cds_boundary import frozen_jinju_records, ranked_model


ROOT = Path(__file__).resolve().parents[2]
POLICY = json.loads(
    (ROOT / "policies/unified_v5_single_exon50_competition.json").read_text()
)["selector_core"]


def select(records):
    selector = UnifiedSelector("ASM", POLICY, sample_count=65, replicates_per_condition=1)
    selector.update(records)
    return selector


def tx(exons, name="rna", strand="+"):
    return RnaTranscript(name, "chr1", strand, tuple(exons))


def index(transcripts):
    return IntervalIndex(transcripts, seqid=lambda t: t.seqid,
                         start=lambda t: t.start, end=lambda t: t.end, strand=lambda t: t.strand)


def inclusion_fixture(strand="+"):
    candidate = CandidateRecord(
        "include", "include", "chr1", strand,
        tuple(CdsPart(a, b, 0) for a, b in ((100, 198), (300, 398), (500, 598))),
        (), 297, "complete")
    # The ends are incomplete, but the entire coding intron chain is present.
    partial = tx(((150, 198), (300, 398), (500, 550)), "include", strand)
    skipping = tx(((90, 198), (500, 650)), "skip", strand)
    evidence = {("chr1", a, b, strand): JunctionAggregate(10, 2, 1, 0)
                for a, b in ((199, 299), (399, 499), (199, 499))}
    return candidate, partial, skipping, evidence


@pytest.mark.parametrize("strand", ["+", "-"])
def test_wholly_intronic_candidate_cannot_escape_observation(strand):
    c = CandidateRecord("c", "c", "chr1", strand, (CdsPart(200, 499, 0),), (), 300, "complete")
    r = tx(((50, 99), (600, 700)), strand=strand)
    evidence = {("chr1", 100, 599, strand): JunctionAggregate(10, 2, 1, 0)}
    metrics = rna_cds_metrics(c, [r], evidence)
    assert metrics["rna_cds_status"] == "CONFLICT"
    assert metrics["rna_cds_unopposed_conflict_fraction"] == 1
    assert metrics["rna_cds_start_conflict"] == metrics["rna_cds_stop_conflict"] == 1
    nested = tx(((190, 550),), "nested", strand)
    metrics = rna_cds_metrics(c, [r, nested], evidence)
    assert metrics["rna_cds_status"] == "MIXED"
    assert metrics["rna_cds_unopposed_conflict_fraction"] == 0


@pytest.mark.parametrize("strand", ["+", "-"])
def test_replicated_inclusion_chain_counters_internal_exon_skipping(strand):
    c, partial, skipping, evidence = inclusion_fixture(strand)
    metrics = rna_cds_metrics(c, [partial, skipping], evidence)
    assert metrics["rna_cds_status"] == "CONFLICT"  # Retain the observed alternative.
    assert metrics["rna_cds_support_transcripts_json"] == "[]"  # No full CDS claim.
    assert metrics["rna_cds_self_strong_chain"] == 1
    assert metrics["rna_cds_unopposed_conflicts_json"] == "[]"
    assert json.loads(metrics["rna_cds_countered_conflicts_json"])[0]["start"] == 199


@pytest.mark.parametrize("missing", ["partial", "replication", "own_junction"])
def test_separate_or_weak_junction_evidence_does_not_claim_a_supported_chain(missing):
    c, partial, skipping, evidence = inclusion_fixture()
    transcripts = [partial, skipping]
    if missing == "partial":
        transcripts.remove(partial)
    elif missing == "replication":
        evidence["chr1", 199, 299, "+"] = JunctionAggregate(100, 5, 0, 0)
    else:
        del evidence["chr1", 199, 299, "+"]
    metrics = rna_cds_metrics(c, transcripts, evidence)
    assert metrics["rna_cds_self_strong_chain"] == 0
    assert metrics["rna_cds_unopposed_conflict_fraction"] > 0


@pytest.mark.parametrize("start,end", [(100, 149), (120, 149), (551, 598)])
def test_partial_chain_cannot_validate_first_last_exon_or_coding_termini(start, end):
    c, partial, _, evidence = inclusion_fixture()
    conflict = tx(((1, start - 1), (end + 1, 700)), "terminal")
    evidence["chr1", start, end, "+"] = JunctionAggregate(30, 3, 1, 0)
    metrics = rna_cds_metrics(c, [partial, conflict], evidence)
    assert metrics["rna_cds_self_strong_chain"] == 1
    assert metrics["rna_cds_unopposed_conflict_fraction"] > 0


def test_stronger_evidence_class_is_not_discarded_for_an_unresolved_weak_candidate():
    strong = ranked_model("strong", 250, conflict=True)
    weak = ranked_model("weak", 220, strong_class=False)
    for records in ([strong, weak], [weak, strong]):
        result = select(records)
        assert [r.key for r in result.selected] == ["strong"]
        assert result.records["strong"].rna_cds_conflict


def test_equal_class_uses_conflict_and_positive_support_before_length():
    strong = ranked_model("strong", 250, conflict=True)
    alternative = ranked_model("alternative", 220)
    assert [r.key for r in select([strong, alternative]).selected] == ["alternative"]
    short = ranked_model("short", 100, full_rna=True)
    long = ranked_model("long", 250)
    assert [r.key for r in select([short, long]).selected] == ["short"]


def test_v5_rejects_old_boundary_evidence_and_inconsistent_fraction():
    record = ranked_model("gene", 250)
    record.row["rna_cds_boundary_version"] = "1"
    with pytest.raises(SelectionError, match="recompute"):
        select([record])
    record.row["rna_cds_boundary_version"] = "2"
    record.row["rna_cds_unopposed_conflict_fraction"] = "1"
    with pytest.raises(SelectionError, match="inconsistent"):
        select([record])


def test_sample_gtf_adds_boundary_support_without_changing_admission_metrics(tmp_path):
    c, partial, skipping, evidence = inclusion_fixture()
    sample = tmp_path / "sample.gtf"
    sample.write_text("".join(
        f'chr1\ttest\texon\t{a}\t{b}\t.\t+\t.\ttranscript_id "full";\n'
        for a, b in ((90, 198), (300, 398), (500, 650))))
    combined = load_rna_cds_transcripts([partial, skipping], [sample, sample])
    assert len(combined) == 3  # Duplicate files/identical chains add no votes.
    before = best_rna_metrics(c, index([partial, skipping]), evidence)
    after = best_rna_metrics(c, index([partial, skipping]), evidence,
                             boundary_index=index(combined))
    assert before["rna_cds_status"] == "CONFLICT"
    assert after["rna_cds_status"] == "MIXED"
    assert str(sample.resolve()) in after["rna_cds_support_transcripts_json"]
    boundary_columns = RNA_CDS_COLUMNS + RNA_CDS_SAMPLE_COLUMNS
    assert {k: v for k, v in before.items() if k not in boundary_columns} == {
        k: v for k, v in after.items() if k not in boundary_columns}


def eight_records(*, samples=True):
    data = load_private_json("jinju_eight_competition.json")
    transcripts = [RnaTranscript(
        r["sample"] + "::" + r["transcript_id"], r["seqid"], r["strand"],
        tuple(tuple(part) for part in json.loads(r["exons"])),
    ) for r in data["transcripts"] if samples or r["sample"] == "merged"]
    evidence = {tuple(r["key"]): JunctionAggregate(
        r["hq_reads"], r["samples"], r["conditions"], r["conditions3"])
        for r in data["junctions"]}
    records = []
    for source in data["models"]:
        row = dict(source)
        parts = parse_cds_signature(row["cds_signature"])
        candidate = CandidateRecord(
            row["candidate_id"], row["structure_id"], row["seqid"], row["strand"],
            parts, tuple(row["predictors"].split(",")), sum(p.length for p in parts), row["orf_status"])
        row.update({k: format_value(v) for k, v in
                    rna_cds_metrics(candidate, transcripts, evidence).items()})
        records.append(Record(row, parts))
    return records, data["pairs"]


def test_real_case_inclusion_survives_without_erasing_skipping_evidence():
    records, pairs = eight_records()
    pair = next(p for p in pairs if p["case"] == "C07")
    result = select(records)
    old = result.records[pair["lost_candidate"]]
    assert old.row["rna_cds_status"] == "CONFLICT"
    assert old.row["rna_cds_support_transcripts_json"] == "[]"
    assert old.row["rna_cds_self_strong_chain"] == "1"
    assert not old.rna_cds_conflict
    assert old.selected
    assert "RNA_CDS_CONFLICT_COUNTERED_BY_OWN_RNA" in old.reasons


def test_real_sample_isoforms_and_enclosed_models_are_both_visible():
    records, pairs = eight_records()
    by_id = {r.key: r for r in records}
    for pair in pairs:
        new = by_id[pair["selected_blocker"]].row
        if pair["case"] in {"C02", "C03", "C04", "C05"}:
            assert float(new["rna_cds_conflict_bp"]) > 0
        if pair["case"] in {"C03", "C05"}:
            assert new["rna_cds_status"] == "MIXED"
            assert json.loads(new["rna_cds_support_transcripts_json"])
            assert new["rna_cds_unopposed_conflicts_json"] == "[]"
        if pair["case"] in {"C02", "C04"}:
            assert new["rna_cds_status"] == "CONFLICT"
            assert float(new["rna_cds_unopposed_conflict_fraction"]) == 1


def test_focal_434_model_remains_selected_by_general_rules():
    result = select(frozen_jinju_records())
    assert result.records["GCAN-f25132743f098767e8652d6a"].selected
    assert not result.records["GCAN-d77d0aa4cf68a62558eb51bc"].selected
    assert not result.records["GCAN-26e832edfbbd3413438ba65a"].selected


def test_replay_order_checkpoint_and_evidence_updates_agree_with_fresh_selection():
    records, _ = eight_records()
    result = select(records)
    audit = {k: audit_row(r, POLICY["policy_id"]) for k, r in result.records.items()}
    reverse = select(list(reversed(records)))
    assert audit == {k: audit_row(r, POLICY["policy_id"]) for k, r in reverse.records.items()}
    restored = select([])
    restored.restore(records, audit)
    old_only, _ = eight_records(samples=False)
    restored.update(old_only)
    fresh = select(old_only)
    assert {k: audit_row(r, POLICY["policy_id"]) for k, r in restored.records.items()} == {
        k: audit_row(r, POLICY["policy_id"]) for k, r in fresh.records.items()}


def test_v5_changes_only_boundary_competition_policy():
    old = json.loads((ROOT / "policies/unified_v4_single_exon50_boundary.json").read_text())
    new = json.loads((ROOT / "policies/unified_v5_single_exon50_competition.json").read_text())
    for key in ("policy_id", "policy_version", "note"):
        new[key] = old[key]
    new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
    assert new["selector_core"].pop("rna_cds_competition_ranking") is True
    assert new == old
