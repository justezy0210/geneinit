rule select_coding_models:
    input:
        code=source_inputs(workflow.basedir, "selection"),
        preflight=PREFLIGHT_COMPLETE,
        components=COMPONENT_EVIDENCE,
        proteins=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.proteins.faa",
        te_evidence=lambda wildcards: [config["te_evidence"]] if config.get("te_evidence") else [],
        component_complete=COMPONENT_COMPLETE,
        policy=lambda wildcards: config["policy"],
    output:
        base_audit=f"{SELECTION_DIR}/balanced_base/{ASSEMBLY_ID}.balanced_base_audit.tsv",
        base_selected=f"{SELECTION_DIR}/balanced_base/{ASSEMBLY_ID}.balanced_base_selected_candidate_ids.txt",
        base_summary=f"{SELECTION_DIR}/balanced_base/{ASSEMBLY_ID}.balanced_base.summary.tsv",
        loci=f"{SELECTION_DIR}/{ASSEMBLY_ID}.candidate_loci.tsv",
        edges=f"{SELECTION_DIR}/{ASSEMBLY_ID}.conflict_edges.tsv",
        audit=f"{SELECTION_DIR}/{ASSEMBLY_ID}.selection_audit.tsv",
        selected=SELECTED_KEYS,
        summary=f"{SELECTION_DIR}/{ASSEMBLY_ID}.selection.summary.tsv",
        complete=SELECTION_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
        sample_count=len(SAMPLES),
        evidence_gate_version="confirmed_te_only_v2",
        replicates_per_condition=REPLICATES_PER_CONDITION,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/selection/core.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/select_coding_models.py"


rule coding_selection:
    input:
        SELECTION_COMPLETE
