rule prepare_selected_coding_models:
    input:
        code=source_inputs(workflow.basedir, "coding_prepare"),
        candidates=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        selection_audit=f"{SELECTION_DIR}/{ASSEMBLY_ID}.selection_audit.tsv",
        selection_complete=SELECTION_COMPLETE,
    output:
        gff3=CODING_GFF3,
        mapping=CODING_MAPPING,
        summary=f"{CODING_DIR}/{ASSEMBLY_ID}.selected_coding.summary.tsv",
        complete=CODING_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        seqids=REFERENCE_SEQUENCE_IDS,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/coding/prepare.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/prepare_selected_coding_models.py"


rule group_coding_models:
    input:
        code=source_inputs(workflow.basedir, "coding_group"),
        gff3=CODING_GFF3,
        mapping=CODING_MAPPING,
        coding_complete=CODING_COMPLETE,
        candidates=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        source_models=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.source_models.tsv",
        policy=lambda wildcards: config["policy"],
    output:
        gff3=GROUPED_GFF3,
        mapping=GROUPING_MAPPING,
        pairs=f"{GROUPING_DIR}/{ASSEMBLY_ID}.provisional_genes.pair_audit.tsv",
        summary=f"{GROUPING_DIR}/{ASSEMBLY_ID}.provisional_genes.summary.tsv",
        complete=GROUPING_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        seqids=REFERENCE_SEQUENCE_IDS,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/coding/grouping.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/group_coding_models.py"


rule coding_models:
    input:
        CODING_COMPLETE


rule gene_grouping:
    input:
        GROUPING_COMPLETE
