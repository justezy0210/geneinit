rule build_candidate_universe:
    input:
        code=source_inputs(workflow.basedir, "models"),
        prediction=PREDICTION_COMPLETE,
        genome=lambda wildcards: config["genome"],
        augustus=f"{AUGUSTUS_RESULT}/augustus.hints.gtf",
        genemark=f"{GENEMARK_RESULT}/genemark.gtf",
        masking_rescue=([MASKING_RESCUE_GTF] if MASKING_RESCUE_ENABLED else []),
    output:
        sources=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.source_models.tsv",
        universe=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        gff3=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.gff3",
        cds=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.cds.fna",
        protein=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.proteins.faa",
        summary=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.summary.tsv",
        complete=CANDIDATE_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/models/candidate_universe.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/build_candidate_universe.py"


rule candidate_universe:
    input:
        CANDIDATE_COMPLETE
