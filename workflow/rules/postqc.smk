rule run_structural_postqc:
    input:
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
    params:
        repo_root=str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
        seqids=REFERENCE_SEQUENCE_IDS,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/iterative_postqc/structural.log",
    conda:
        "../../envs/evidence.yaml"
    script:
        "../scripts/run_structural_postqc.py"


rule structural_postqc:
    input:
        POSTQC_COMPLETE
