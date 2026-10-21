rule build_structural_release:
    input:
        gff3=UTR_GFF3,
        te_evidence=lambda wildcards: [config["te_evidence"]] if config.get("te_evidence") else [],
        candidates=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        selection_audit=f"{SELECTION_DIR}/{ASSEMBLY_ID}.selection_audit.tsv",
        generated_te=f"{POSTQC_DIR}/{ASSEMBLY_ID}.te_classification.json",
        transcript_decisions=f"{POSTQC_DIR}/{ASSEMBLY_ID}.transcript_decisions.tsv",
        edge_audit=f"{POSTQC_DIR}/{ASSEMBLY_ID}.coding_edge_audit.tsv",
        bridge_reviews=f"{POSTQC_DIR}/{ASSEMBLY_ID}.bridge_reviews.tsv",
        postqc_complete=POSTQC_COMPLETE,
        utr_complete=UTR_COMPLETE,
        genome=lambda wildcards: config["genome"],
        policy=lambda wildcards: config["policy"],
    output:
        gff3=f"{RELEASE_STAGING_DIR}/annotation.gff3",
        te_gff3=f"{RELEASE_STAGING_DIR}/te_annotation.gff3",
        te_cds=f"{RELEASE_STAGING_DIR}/te_cds.fna",
        te_proteins=f"{RELEASE_STAGING_DIR}/te_proteins.faa",
        te_audit=f"{RELEASE_STAGING_DIR}/te_candidate_audit.tsv",
        te_summary=f"{RELEASE_STAGING_DIR}/te_summary.json",
        transcripts_fna=f"{RELEASE_STAGING_DIR}/transcripts.fna",
        cds_fna=f"{RELEASE_STAGING_DIR}/cds.fna",
        proteins_faa=f"{RELEASE_STAGING_DIR}/proteins.faa",
        longest_proteins_faa=f"{RELEASE_STAGING_DIR}/longest_proteins.faa",
        gene_id_map=f"{RELEASE_STAGING_DIR}/gene_id_map.tsv",
        transcript_id_map=f"{RELEASE_STAGING_DIR}/transcript_id_map.tsv",
        candidate_lineage=f"{RELEASE_STAGING_DIR}/candidate_lineage.tsv",
        structural_decisions=f"{RELEASE_STAGING_DIR}/structural_decisions.tsv",
        review_flags=f"{RELEASE_STAGING_DIR}/review_flags.tsv",
        qc_tsv=f"{RELEASE_STAGING_DIR}/structural_qc.tsv",
        qc_json=f"{RELEASE_STAGING_DIR}/structural_qc.json",
        release_status=f"{RELEASE_STAGING_DIR}/release_status.json",
        manifest=f"{RELEASE_STAGING_DIR}/release_manifest.json",
        checksums=f"{RELEASE_STAGING_DIR}/checksums.sha256",
        policy=f"{RELEASE_STAGING_DIR}/policy.json",
        complete=RELEASE_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        seqids=REFERENCE_SEQUENCE_IDS,
        assembly_id=ASSEMBLY_ID,
        gene_prefix=str(config["gene_prefix"]),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/release/structural.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/build_structural_release.py"


rule structural_release:
    input:
        RELEASE_COMPLETE
