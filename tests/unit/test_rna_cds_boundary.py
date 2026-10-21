from copy import deepcopy
import json
from pathlib import Path

import pytest

from private_fixtures import load_private_json

from geneinit.evidence.components import (
    CandidateRecord, JunctionAggregate, RnaTranscript, best_rna_metrics,
    format_value, parse_cds_signature,
)
from geneinit.evidence.metrics import IntervalIndex
from geneinit.evidence.rna_cds import (
    RNA_CDS_COLUMNS, rna_cds_metrics, terminal_positions, transcript_cds_structure,
)
from geneinit.models.candidate import CdsPart
from geneinit.selection.core import Record, SelectionError, audit_row
from geneinit.selection.unified import UnifiedSelector
from test_unified_selection import make_record


ROOT = Path(__file__).resolve().parents[2]
POLICY = json.loads(
    (ROOT / "policies/unified_v4_single_exon50_boundary.json").read_text()
)["selector_core"]


def candidate(strand="+"):
    return CandidateRecord("gene", "gene", "chr1", strand,
                           (CdsPart(100, 399, 0),), ("augustus",), 300, "complete")


def transcript(exons, strand="+", name="rna"):
    return RnaTranscript(name, "chr1", strand, tuple(exons))


def junctions(start=399, end=449, strand="+", reads=3, samples=2, conditions=1):
    return {("chr1", start, end, strand): JunctionAggregate(reads, samples, conditions, 0)}


def record_with_rna(rnas, evidence, strand="+"):
    record = make_record("gene", ((100, 399, 0),), protein_aa=99,
                         protein_tier=4, protein_unique=True,
                         protein_precision=1, terminals=2)
    record.row["strand"] = strand
    record.row.update({k: format_value(v) for k, v in
                       rna_cds_metrics(candidate(strand), rnas, evidence).items()})
    return record


def select(records, policy=POLICY):
    selector = UnifiedSelector("ASM", policy, sample_count=6, replicates_per_condition=2)
    selector.update(records)
    return selector


@pytest.mark.parametrize("strand,terminal", [("+", "stop"), ("-", "start")])
def test_one_base_boundary_crossing_is_not_hidden_by_99_percent_overlap(strand, terminal):
    rna = transcript(((50, 398), (450, 550)), strand)
    metrics = transcript_cds_structure(candidate(strand), rna)
    assert metrics["cds_exonic_bp"] == 299
    assert metrics["cds_bases_inside_rna_introns"] == 1
    assert not metrics["full_cds_compatible"]
    record = record_with_rna([rna], junctions(strand=strand), strand)
    assert record.row[f"rna_cds_{terminal}_conflict"] == "1"
    result = select([record]).records["gene"]
    assert result.hard_qc
    assert result.eligible and result.selected
    assert result.decision == "SELECTED"
    assert result.rna_cds_conflict
    assert "SELECTED_WITH_UNRESOLVED_RNA_CDS_CONFLICT" in result.reasons
    # A sole eligible model can be returned, with its contradiction exposed.
    assert result.strong_unique_full_protein


@pytest.mark.parametrize("strand,terminal", [("+", "start"), ("-", "stop")])
def test_intron_crossing_low_boundary_is_detected(strand, terminal):
    rna = transcript(((1, 49), (102, 450)), strand)
    metrics = rna_cds_metrics(candidate(strand), [rna], junctions(50, 101, strand))
    assert metrics["rna_cds_status"] == "CONFLICT"
    assert metrics[f"rna_cds_{terminal}_conflict"] == 1


@pytest.mark.parametrize("strand", ["+", "-"])
def test_utr_splice_does_not_conflict_with_single_coding_exon(strand):
    rna = transcript(((1, 49), (100, 450)), strand)
    record = record_with_rna([rna], junctions(50, 99, strand), strand)
    assert record.row["rna_cds_status"] == "SUPPORTED"
    assert record.row["rna_cds_conflict_bp"] == "0"
    assert select([record]).selected


@pytest.mark.parametrize("rnas", [
    [], [transcript(((150, 250),))],
    [transcript(((1, 49), (102, 450)), ".")],
    [transcript(((1, 49), (102, 450)), "-")],
])
def test_missing_partial_unknown_or_opposite_strand_rna_does_not_reject(rnas):
    record = record_with_rna(rnas, junctions(50, 101))
    assert record.row["rna_cds_status"] == "UNRESOLVED"
    assert select([record]).selected


@pytest.mark.parametrize("reads,samples,conditions", [(2, 2, 1), (10, 1, 0), (10, 3, 0)])
def test_weak_or_unreplicated_junction_does_not_create_a_penalty(reads, samples, conditions):
    rna = transcript(((50, 398), (450, 550)))
    record = record_with_rna([rna], junctions(reads=reads, samples=samples, conditions=conditions))
    assert record.row["rna_cds_status"] == "UNRESOLVED"
    assert select([record]).selected


def test_full_alternative_isoform_preserves_both_support_and_conflict():
    bad = transcript(((50, 398), (450, 550)), name="spliced")
    good = transcript(((50, 550),), name="retained")
    record = record_with_rna([bad, good], junctions())
    assert record.row["rna_cds_status"] == "MIXED"
    assert json.loads(record.row["rna_cds_support_transcripts_json"]) == ["retained"]
    assert json.loads(record.row["rna_cds_conflict_transcripts_json"]) == ["spliced"]
    assert select([record]).selected
    assert rna_cds_metrics(candidate(), [good, bad], junctions()) == rna_cds_metrics(
        candidate(), [bad, good], junctions()
    )


def test_split_codons_follow_transcription_order_on_both_strands():
    parts = (CdsPart(100, 101, 0), CdsPart(200, 206, 0))
    assert terminal_positions(parts, "+") == ((100, 101, 200), (204, 205, 206))
    assert terminal_positions(parts, "-") == ((206, 205, 204), (200, 101, 100))


def test_cds_intron_missing_from_rna_is_not_full_structure_support():
    c = CandidateRecord("c", "c", "chr1", "+",
                        (CdsPart(100, 150, 0), CdsPart(201, 299, 0)),
                        ("augustus",), 150, "complete")
    rna = transcript(((90, 350),))
    metrics = transcript_cds_structure(c, rna)
    assert metrics["cds_exonic_bp"] == 150
    assert metrics["candidate_introns_missing_from_rna"] == ((151, 200),)
    assert not metrics["full_cds_compatible"]


def test_old_component_tables_fail_closed_only_under_the_new_policy():
    record = record_with_rna([], {})
    for key in RNA_CDS_COLUMNS:
        del record.row[key]
    with pytest.raises(SelectionError, match="recompute components"):
        select([record])
    legacy = deepcopy(POLICY)
    legacy.pop("rna_cds_boundary_ranking")
    assert select([record], legacy).selected


def test_conflict_flag_survives_checkpoint_and_disappears_with_updated_evidence():
    bad = record_with_rna([transcript(((50, 398), (450, 550)))], junctions())
    initial = select([bad])
    audits = {k: audit_row(r, POLICY["policy_id"]) for k, r in initial.records.items()}
    restored = select([])
    restored.restore([bad], audits)
    restored.update([bad])
    assert restored.records["gene"].rna_cds_conflict
    assert restored.selected
    good = record_with_rna([transcript(((50, 550),))], {})
    restored.update([good])
    assert restored.selected
    assert not restored.records["gene"].rna_cds_conflict
    assert restored.stats.evaluated_candidates == 1


def test_obsolete_review_only_policy_is_not_silently_reinterpreted():
    obsolete = deepcopy(POLICY)
    obsolete["rna_cds_boundary_review"] = obsolete.pop("rna_cds_boundary_ranking")
    with pytest.raises(SelectionError, match="obsolete"):
        select([], obsolete)


def test_policy_keeps_length_orf_te_and_admission_rules_unchanged():
    old = json.loads((ROOT / "policies/unified_v3_single_exon50.json").read_text())
    new = json.loads((ROOT / "policies/unified_v4_single_exon50_boundary.json").read_text())
    for key in ("policy_id", "policy_version", "note"):
        new[key] = old[key]
    new["selector_core"]["policy_id"] = old["selector_core"]["policy_id"]
    assert new["selector_core"].pop("rna_cds_boundary_ranking") is True
    assert new == old


def frozen_jinju_records():
    fixture = load_private_json("jinju_rna_boundary.json")
    rnas = [RnaTranscript(
        r["transcript_id"], r["seqid"], r["strand"],
        tuple(tuple(p) for p in json.loads(r["exons"])),
    ) for r in fixture["transcripts"]]
    index = IntervalIndex(rnas, seqid=lambda r: r.seqid, start=lambda r: r.start,
                          end=lambda r: r.end, strand=lambda r: r.strand)
    evidence = {
        tuple(r["key"]): JunctionAggregate(r["hq_reads"], r["samples"], r["conditions"], 0)
        for r in fixture["junctions"]
    }
    records = []
    for source in fixture["models"]:
        row = dict(source)
        parts = parse_cds_signature(row["cds_signature"])
        c = CandidateRecord(row["candidate_id"], row["structure_id"], row["seqid"],
                            row["strand"], parts, tuple(row["predictors"].split(",")),
                            sum(p.length for p in parts), row["orf_status"])
        # Keep the original admission metrics; only add the new boundary evidence.
        metrics = best_rna_metrics(c, index, evidence)
        row.update({k: format_value(metrics[k]) for k in RNA_CDS_COLUMNS})
        records.append(Record(row, parts))
    return records


def test_real_jinju_chooses_stronger_nonconflicting_candidate_automatically():
    selector = select(frozen_jinju_records())
    short = selector.records["GCAN-26e832edfbbd3413438ba65a"]
    long = selector.records["GCAN-d77d0aa4cf68a62558eb51bc"]
    assert short.row["rna_cds_status"] == "SUPPORTED"
    assert short.eligible
    assert short.decision == "REJECT_INCOMPATIBLE_GRAPH"
    assert not short.selected
    assert long.row["rna_cds_status"] == "CONFLICT"
    assert long.row["rna_cds_stop_conflict"] == "1"
    assert long.rna_cds_conflict
    assert long.decision == "REJECT_INCOMPATIBLE_GRAPH"
    assert not long.selected
    old = selector.records["GCAN-f25132743f098767e8652d6a"]
    assert old.row["rna_cds_status"] == "UNRESOLVED"
    assert old.eligible
    assert old.decision == "SELECTED"
    assert old.selected
    assert old.key in long.blockers
    assert old.key in short.blockers


def test_real_jinju_original_decisions_are_reproduced_without_boundary_ranking():
    fixture = load_private_json("jinju_rna_boundary.json")
    legacy = deepcopy(POLICY)
    legacy.pop("rna_cds_boundary_ranking")
    selector = select(frozen_jinju_records(), legacy)
    for ident, expected in fixture["original_selection"].items():
        assert selector.records[ident].decision == expected["decision"]
        assert selector.records[ident].selected == bool(int(expected["selected"]))


def test_real_jinju_rankings_and_unrelated_loci_survive_checkpoint_restoration():
    records = frozen_jinju_records()
    original = select(records)
    audits = {k: audit_row(r, POLICY["policy_id"]) for k, r in original.records.items()}
    restored = select([])
    restored.restore(records, audits)
    restored.update([r for r in records if r.key == "GCAN-d77d0aa4cf68a62558eb51bc"])
    assert {k: audit_row(r, POLICY["policy_id"]) for k, r in restored.records.items()} == audits


def ranked_model(key, aa, *, full_rna=False, conflict=False, strong_class=True):
    end = 100 + 3 * (aa + 1) - 1
    record = make_record(
        key, ((100, end, 0),), protein_aa=aa,
        protein_tier=4, protein_unique=strong_class,
        protein_precision=1, terminals=2 if strong_class else 1,
    )
    c = CandidateRecord(key, key, "chr1", "+", record.cds, ("augustus",),
                        3 * (aa + 1), "complete")
    rnas = [transcript(((90, end + 10),))] if full_rna else []
    evidence = {}
    if conflict:
        rnas.append(transcript(((90, end - 1), (end + 50, end + 100)), name="conflict"))
        evidence = junctions(end, end + 49)
    record.row.update({k: format_value(v) for k, v in rna_cds_metrics(c, rnas, evidence).items()})
    return record


@pytest.mark.parametrize("long_conflicts", [True, False])
def test_conflict_penalty_can_favor_short_or_long_without_a_length_exception(long_conflicts):
    short = ranked_model("short", 84, conflict=not long_conflicts)
    long = ranked_model("long", 434, conflict=long_conflicts)
    expected = "short" if long_conflicts else "long"
    for records in ([short, long], [long, short]):
        assert [r.key for r in select(records).selected] == [expected]


@pytest.mark.parametrize("long_supported", [True, False])
def test_full_rna_breaks_a_tie_within_the_same_evidence_class(long_supported):
    records = [
        ranked_model("short", 84, full_rna=not long_supported),
        ranked_model("long", 434, full_rna=long_supported),
    ]
    assert [r.key for r in select(records).selected] == [
        "long" if long_supported else "short"
    ]


def test_complete_short_rna_does_not_override_a_stronger_protein_evidence_class():
    short = ranked_model("short", 84, full_rna=True, strong_class=False)
    long = ranked_model("long", 434)
    assert [r.key for r in select([short, long]).selected] == ["long"]


def test_all_conflicting_candidates_still_return_the_best_eligible_model_with_a_flag():
    short = ranked_model("short", 84, conflict=True)
    long = ranked_model("long", 434, conflict=True)
    result = select([short, long])
    assert [r.key for r in result.selected] == ["long"]
    row = audit_row(result.selected[0], POLICY["policy_id"])
    assert row["rna_cds_conflict_penalty"] == 1
    assert json.loads(row["selection_rank_json"])[0] == 1
    assert "SELECTED_WITH_UNRESOLVED_RNA_CDS_CONFLICT" in row["reasons"]


def test_changed_competitor_evidence_replaces_a_winner_after_checkpoint_restore():
    records = [ranked_model("short", 84, conflict=True), ranked_model("long", 434)]
    original = select(records)
    audits = {k: audit_row(r, POLICY["policy_id"]) for k, r in original.records.items()}
    restored = select([])
    restored.restore(records, audits)
    restored.update([
        ranked_model("short", 84, full_rna=True), ranked_model("long", 434, conflict=True),
    ])
    assert [r.key for r in restored.selected] == ["short"]
