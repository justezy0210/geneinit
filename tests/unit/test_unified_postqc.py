import csv
import gzip
import hashlib
import json
from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest

from geneinit.evidence.components import (
    CandidateRecord, JunctionAggregate, RnaTranscript, score_candidates, write_evidence_table,
)
from geneinit.evidence.junctions import AGGREGATE_COLUMNS
from geneinit.models.candidate import CdsPart, audit_orf
from geneinit.postqc.core import PostQcError, run_postqc
from geneinit.postqc.unified import RepairEvidence
from geneinit.selection.unified import UnifiedSelector
from geneinit.utr.core import load_raw_gff3, prepare_selected_coding_models
from geneinit.utr.grouping import group_coding_models
from geneinit.release.core import build_structural_release


ROOT = Path(__file__).resolve().parents[2]


def table(path, columns, rows):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def make_inputs(tmp_path, policy_name="unified_v3.json"):
    genome = list("A" * 2000)
    for start, seq in [
        (100, "ATG"), (510, "TAA"), (700, "ATG"), (1197, "TAA"),
        (400, "GT"), (448, "AG"), (763, "GT"), (898, "AG"),
    ]:
        genome[start-1:start-1+len(seq)] = seq
    genome = {"chr1": "".join(genome)}
    candidates = [
        CandidateRecord("left", "left", "chr1", "+",
                        (CdsPart(100, 399, 0), CdsPart(450, 512, 0)),
                        ("augustus",), 363, "complete"),
        CandidateRecord("right", "right", "chr1", "+",
                        (CdsPart(700, 762, 0), CdsPart(900, 1199, 0)),
                        ("augustus",), 363, "complete"),
    ]
    junctions = {
        ("chr1", left, right, "+"): JunctionAggregate(10, 2, 1, 1)
        for left, right in ((400, 449), (763, 899), (400, 899))
    }
    evidence_rows = score_candidates(candidates, [], [], [], junctions, genome, {}, {}, 2)
    evidence = tmp_path / "evidence.tsv"
    write_evidence_table(evidence, evidence_rows)
    candidate_columns = (
        "candidate_id", "structure_id", "seqid", "strand", "start", "end",
        "cds_signature", "predictors", "source_gene_ids", "source_transcript_ids",
    )
    candidate_table = table(
        tmp_path / "candidates.tsv",
        candidate_columns,
        [{key: row.get(key, row["candidate_id"]) for key in candidate_columns}
         for row in evidence_rows],
    )
    proteins = tmp_path / "proteins.faa"
    proteins.write_text("".join(
        f">{candidate.candidate_id}\n"
        f"{audit_orf(candidate.seqid, candidate.strand, candidate.cds, genome).protein_sequence}\n"
        for candidate in candidates
    ))
    outputs = SimpleNamespace(**{
        key: str(tmp_path / f"{key}.tsv")
        for key in ("base_audit", "base_selected", "base_summary", "loci",
                    "edges", "audit", "selected", "summary", "complete")
    })
    policy = tmp_path / "policy.json"
    policy.write_text((ROOT / "policies" / policy_name).read_text())
    runpy.run_path(
        str(ROOT / "workflow/scripts/select_coding_models.py"),
        init_globals={"snakemake": SimpleNamespace(
            params=SimpleNamespace(repo_root=ROOT, assembly_id="ASM",
                                   sample_count=2, replicates_per_condition=2),
            input=SimpleNamespace(policy=policy, proteins=proteins,
                                  components=evidence, te_evidence=[]),
            output=outputs,
        )},
    )
    with Path(outputs.audit).open() as handle:
        audits = {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert all(row["selected"] == "1" for row in audits.values())
    coding, coding_map = tmp_path / "coding.gff3", tmp_path / "coding.tsv"
    prepare_selected_coding_models(
        candidate_table, Path(outputs.audit), ["chr1"], coding, coding_map,
    )
    grouped, grouping = tmp_path / "grouped.gff3", tmp_path / "grouping.tsv"
    group_coding_models(
        coding, coding_map, candidate_table, tmp_path / "unused_source_models.tsv",
        ["chr1"], grouped, grouping, tmp_path / "group_pairs.tsv",
        tmp_path / "group_summary.tsv", tmp_path / "group_complete.tsv",
        defer_grouping=True,
    )
    junction_path = tmp_path / "junctions.tsv.gz"
    with gzip.open(junction_path, "wt") as handle:
        handle.write("#" + "\t".join(AGGREGATE_COLUMNS) + "\n")
        for (seqid, start, end, strand), aggregate in sorted(junctions.items()):
            row = {column: 0 for column in AGGREGATE_COLUMNS}
            row.update(seqid=seqid, intron_start=start, intron_end=end, strand=strand,
                       total_hq_reads=10, sample_count_hq_ge1=2,
                       condition_count_ge2rep_hq_ge1=1, condition_count_ge2rep_hq_ge3=1)
            handle.write("\t".join(str(row[column]) for column in AGGREGATE_COLUMNS) + "\n")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + genome["chr1"] + "\n")
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    context = {
        "selection_complete": outputs.complete,
        "sample_count": 2, "replicates_per_condition": 2, "expected_protein_copies": 2,
        "paths": {key: str(empty) for key in (
            "merged_gtf", "direct_chains", "direct_hits", "nonhc_chains",
            "nonhc_hits", "nonhc_trace", "repeats",
        )},
    }
    for key, content in {
        "nonhc_trace": "1\tchr1\t1\t2000\n",
        "direct_chains": "chr1\tProtHint\tCDSpart\t1500\t1530\t.\t+\t0\tgrp=prot_MSTRG.1;\n",
        "nonhc_chains": "1\tProtHint\tCDSpart\t1600\t1630\t.\t+\t0\tgrp=prot_1_g;\n",
        "direct_hits": "MSTRG.1\tprot\n",
        "nonhc_hits": "1_g\tprot\n",
    }.items():
        path = tmp_path / f"{key}.input"
        path.write_text(content)
        context["paths"][key] = str(path)
    args = [
        grouped, grouping, candidate_table, evidence, Path(outputs.audit),
        junction_path, fasta, policy, ["chr1"],
        *[tmp_path / name for name in (
            "postqc.gff3", "transcript_decisions.tsv", "gene_decisions.tsv",
            "postqc_edges.tsv", "promotions.tsv", "reviews.tsv", "retirements.tsv",
            "postqc_summary.tsv", "postqc_complete.tsv", "iterations.tsv",
        )],
        None, "ASM", tmp_path / "te_classification.json", context,
    ]
    return args, genome, junctions


@pytest.mark.parametrize("policy_name", [
    "unified_v3.json", "unified_v4_single_exon50_boundary.json",
    "unified_v5_single_exon50_competition.json",
])
def test_workflow_selection_and_uniform_repair_reach_a_fixed_point(tmp_path, policy_name):
    args, genome, junctions = make_inputs(tmp_path, policy_name)
    run_postqc(*args)
    with (tmp_path / "postqc_summary.tsv").open() as handle:
        summary = {r["metric"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}
    assert summary["selection_engine"] == "unified_v1"
    assert summary["input_transcript_count"] == "2"
    assert summary["output_transcript_count"] == "1"
    assert summary["retired_transcript_count"] == "2"
    assert summary["new_candidate_evaluations"] == "1"
    assert summary["cached_candidate_evaluations"] == "2"
    assert summary["converged_iteration"] == "2"
    _, models = load_raw_gff3(tmp_path / "postqc.gff3")
    assert len(models) == 1
    model = next(iter(models.values()))
    parts = tuple(CdsPart(f.start, f.end, int(f.phase)) for f in model.cds)
    assert audit_orf("chr1", "+", parts, genome).status == "complete"
    scorer = RepairEvidence(transcripts=[], chains=[], repeats=[], retained_targets={},
                            seed_targets={}, genome=genome, junctions=junctions,
                            expected_protein_copies=2)
    scored = scorer.score(model)
    assert scored.row["predictor_tier"] == "0"
    assert scored.row["predictors"] == ""
    assert scored.row["protein_tier"] == "0"
    assert scored.row["rna_tier"] == "4"
    assert scored.row["rna_cds_boundary_version"] == "2"
    assert scored.row["rna_cds_status"] == "UNRESOLVED"
    if policy_name != "unified_v3.json":
        # A newly reconstructed model must use the same boundary checks as an
        # original candidate, including an intron crossing its terminal base.
        rna = RnaTranscript("conflicting_repair_rna", "chr1", "+", (
            (parts[0].start, parts[0].end),
            (parts[-1].start, parts[-1].end - 1),
            (parts[-1].end + 101, parts[-1].end + 150),
        ))
        conflicting_junctions = {
            **junctions,
            ("chr1", parts[-1].end, parts[-1].end + 100, "+"):
                JunctionAggregate(10, 2, 1, 1),
        }
        scorer = RepairEvidence(
            transcripts=[rna], chains=[], repeats=[], retained_targets={},
            seed_targets={}, genome=genome, junctions=conflicting_junctions,
            expected_protein_copies=2,
        )
        repaired = scorer.score(model)
        assert repaired.row["rna_cds_status"] == "CONFLICT"
        assert repaired.row["rna_cds_stop_conflict"] == "1"
        selector = UnifiedSelector(
            "ASM", json.loads(args[7].read_text())["selector_core"],
            sample_count=2, replicates_per_condition=2,
        )
        selector.update([repaired])
        assert selector.selected
        assert selector.records[repaired.key].rna_cds_conflict
        assert "SELECTED_WITH_UNRESOLVED_RNA_CDS_CONFLICT" in selector.records[repaired.key].reasons
    release = tmp_path / "release"
    release.mkdir()
    build_structural_release(
        tmp_path / "postqc.gff3", tmp_path / "transcript_decisions.tsv",
        tmp_path / "postqc_edges.tsv", tmp_path / "reviews.tsv",
        args[6], args[7], ["chr1"], "ASM", "TEST",
        *[release / name for name in (
            "annotation.gff3", "transcripts.fna", "cds.fna", "proteins.faa",
            "longest_proteins.faa", "gene_id_map.tsv", "transcript_id_map.tsv",
            "candidate_lineage.tsv", "structural_decisions.tsv", "review_flags.tsv",
            "structural_qc.tsv", "structural_qc.json", "release_status.json",
            "release_manifest.json", "checksums.sha256", "policy.json",
            "release.complete.tsv",
        )],
        candidates_path=args[2], selection_audit_path=args[4],
        generated_te_path=tmp_path / "te_classification.json",
    )
    assert (release / "proteins.faa").read_text().count(">") == 1
    assert json.loads((release / "policy.json").read_text())["status"] == "EXPERIMENTAL_NOT_VALIDATED"
    for line in (release / "checksums.sha256").read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        assert hashlib.sha256((release / filename).read_bytes()).hexdigest() == expected


def test_policy_change_cannot_reuse_old_selection_checkpoint(tmp_path):
    args, _, _ = make_inputs(tmp_path)
    policy_path = args[7]
    policy = json.loads(policy_path.read_text())
    policy["selector_core"]["hard_qc"]["min_protein_aa"] = 80
    policy_path.write_text(json.dumps(policy))
    with pytest.raises(PostQcError, match="SHA256 mismatch"):
        run_postqc(*args)


def test_repair_from_paths_uses_sample_rna_without_inheriting_parent_scores(tmp_path):
    args, genome, junctions = make_inputs(tmp_path, "unified_v5_single_exon50_competition.json")
    run_postqc(*args)
    _, models = load_raw_gff3(tmp_path / "postqc.gff3")
    model = next(iter(models.values()))
    exons = [(f.start, f.end) for f in sorted(model.cds, key=lambda f: f.start)]
    sample = tmp_path / "sample.gtf"
    sample.write_text("".join(
        f'chr1\ttest\texon\t{a}\t{b}\t.\t+\t.\ttranscript_id "repair_sample";\n'
        for a, b in exons))
    paths = dict(args[-1]["paths"])
    before = RepairEvidence.from_paths(
        paths, genome=genome, junctions=junctions, expected_protein_copies=2).score(model)
    paths["sample_gtfs"] = [str(sample)]
    after = RepairEvidence.from_paths(
        paths, genome=genome, junctions=junctions, expected_protein_copies=2).score(model)
    assert before.row["rna_cds_status"] == "UNRESOLVED"
    assert after.row["rna_cds_status"] == "SUPPORTED"
    assert "repair_sample" in after.row["rna_cds_support_transcripts_json"]
    assert before.row["rna_tier"] == after.row["rna_tier"]
    assert after.row["predictor_tier"] == "0"
    selector = UnifiedSelector(
        "ASM", json.loads(args[7].read_text())["selector_core"],
        sample_count=2, replicates_per_condition=2)
    selector.update([after])
    assert selector.selected
