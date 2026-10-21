rule run_structural_postqc:
    input:
        code=source_inputs(workflow.basedir, "postqc"),
        gff3=GROUPED_GFF3,
        te_evidence=lambda wildcards: [config["te_evidence"]] if config.get("te_evidence") else [],
        grouping_mapping=GROUPING_MAPPING,
        grouping_complete=GROUPING_COMPLETE,
        candidates=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        component_evidence=COMPONENT_EVIDENCE,
        selection_audit=f"{SELECTION_DIR}/{ASSEMBLY_ID}.selection_audit.tsv",
        junctions=JUNCTION_AGGREGATE,
        junction_index=f"{JUNCTION_AGGREGATE}.csi",
        genome=lambda wildcards: config["genome"],
        policy=lambda wildcards: config["policy"],
        selection_complete=SELECTION_COMPLETE,
        merged_gtf=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
        sample_gtfs=expand(
            f"{TRANSCRIPTS_DIR}/samples/{{sample}}/{ASSEMBLY_ID}.{{sample}}.stringtie.gtf",
            sample=SAMPLES,
        ),
        direct_chains=f"{GENEMARK_RESULT}/selector_evidence/direct_top_chains.gff",
        direct_hits=f"{GENEMARK_RESULT}/selector_evidence/direct_protein_hits.tsv",
        nonhc_chains=f"{GENEMARK_RESULT}/selector_evidence/nonhc_top_chains.gff",
        nonhc_hits=f"{GENEMARK_RESULT}/selector_evidence/nonhc_protein_hits.tsv",
        nonhc_trace=f"{GENEMARK_RESULT}/selector_evidence/nonhc.trace",
        repeats=lambda wildcards: config["repeat_gff"],
    output:
        gff3=POSTQC_GFF3,
        te_classification=f"{POSTQC_DIR}/{ASSEMBLY_ID}.te_classification.json",
        transcripts=f"{POSTQC_DIR}/{ASSEMBLY_ID}.transcript_decisions.tsv",
        genes=f"{POSTQC_DIR}/{ASSEMBLY_ID}.gene_decisions.tsv",
        edges=f"{POSTQC_DIR}/{ASSEMBLY_ID}.coding_edge_audit.tsv",
        promotions=f"{POSTQC_DIR}/{ASSEMBLY_ID}.bridge_promotions.tsv",
        bridge_reviews=f"{POSTQC_DIR}/{ASSEMBLY_ID}.bridge_reviews.tsv",
        retirements=f"{POSTQC_DIR}/{ASSEMBLY_ID}.connector_retirements.tsv",
        iterations=f"{POSTQC_DIR}/{ASSEMBLY_ID}.iterations.tsv",
        summary=f"{POSTQC_DIR}/{ASSEMBLY_ID}.postqc.summary.tsv",
        complete=POSTQC_COMPLETE,
        unified_selection=str(Path(POSTQC_GFF3).with_suffix(".unified_selection.tsv")),
        unified_rounds=str(Path(POSTQC_GFF3).with_suffix(".unified_rounds.tsv")),
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
        seqids=REFERENCE_SEQUENCE_IDS,
        sample_count=len(SAMPLES),
        replicates_per_condition=REPLICATES_PER_CONDITION,
        expected_protein_copies=int(config["expected_protein_copies"]),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/iterative_postqc/structural.log",
    conda:
        "../../envs/evidence.yaml"
    script:
        "../scripts/run_structural_postqc.py"


rule structural_postqc:
    input:
        POSTQC_COMPLETE
