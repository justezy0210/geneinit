rule predictor_container_preflight:
    input:
        code=source_inputs(workflow.basedir, "predictor_container"),
        preflight=INPUT_PREFLIGHT_COMPLETE,
        sif=lambda wildcards: _deployment_value("predictor_sif"),
    output:
        versions=f"{PREDICTION_DIR}/00_container/software_versions.txt",
        complete=PREDICTOR_PREFLIGHT,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        deployment_snapshot=lambda wildcards: json.dumps(config.get("deployment", {}), sort_keys=True),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/container_preflight.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/predictor_container_preflight.py"


rule genemark_etp:
    input:
        code=source_inputs(workflow.basedir, "genemark"),
        container=PREDICTOR_PREFLIGHT,
        genome=lambda wildcards: config["genome"],
        proteins=lambda wildcards: config["protein_fasta"],
        merged_gtf=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
        bams=SAMPLE_BAMS,
        csis=SAMPLE_CSIS,
    output:
        work=directory(GENEMARK_WORK),
        genemark=f"{GENEMARK_RESULT}/genemark.gtf",
        supported=f"{GENEMARK_RESULT}/genemark_supported.gtf",
        rna_hints=f"{GENEMARK_RESULT}/rnaseq_hints.gff",
        protein_hints=f"{GENEMARK_RESULT}/gmst_protein_hints.gff",
        protein_evidence=f"{GENEMARK_RESULT}/gmst_protein_evidence.gff",
        direct_top_chains=f"{GENEMARK_RESULT}/selector_evidence/direct_top_chains.gff",
        direct_protein_hits=f"{GENEMARK_RESULT}/selector_evidence/direct_protein_hits.tsv",
        nonhc_top_chains=f"{GENEMARK_RESULT}/selector_evidence/nonhc_top_chains.gff",
        nonhc_protein_hits=f"{GENEMARK_RESULT}/selector_evidence/nonhc_protein_hits.tsv",
        nonhc_trace=f"{GENEMARK_RESULT}/selector_evidence/nonhc.trace",
        training=f"{GENEMARK_RESULT}/training.gtf",
        model=f"{GENEMARK_RESULT}/initial_genemark_model.mod",
        complete=GENEMARK_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        samples=SAMPLES,
        assembly_id=ASSEMBLY_ID,
        species=ASSEMBLY_ID,
        proc=GENEMARK_PROC,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/genemark_etp.log",
    threads: workflow.cores
    conda:
        "../../envs/predictor-host.yaml"
    script:
        "../scripts/run_genemark_etp.py"


rule build_full_etp_hints:
    input:
        code=source_inputs(workflow.basedir, "prediction_hints"),
        container=PREDICTOR_PREFLIGHT,
        genemark_complete=GENEMARK_COMPLETE,
        work=GENEMARK_WORK,
        genome=lambda wildcards: config["genome"],
    output:
        hints=ETP_HINTS,
        complete=ETP_HINTS_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        proc=GENEMARK_PROC,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/etp_hints.log",
    conda:
        "../../envs/predictor-host.yaml"
    script:
        "../scripts/build_etp_hints.py"


rule prepare_augustus_config:
    input:
        code=source_inputs(workflow.basedir, "augustus_config"),
        container=PREDICTOR_PREFLIGHT,
    output:
        config=directory(AUGUSTUS_CONFIG),
        complete=f"{PREDICTION_DIR}/03_augustus/config.complete.tsv",
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/augustus_config.log",
    conda:
        "../../envs/predictor-host.yaml"
    script:
        "../scripts/prepare_augustus_config.py"


rule augustus_direct:
    input:
        code=source_inputs(workflow.basedir, "augustus"),
        prediction_resume=([AUGUSTUS_PREDICTION_RESUME] if AUGUSTUS_PREDICTION_RESUME else []),
        optimization_start=([AUGUSTUS_OPTIMIZATION_START] if AUGUSTUS_OPTIMIZATION_START else []),
        training_snapshot=([AUGUSTUS_TRAINING_SNAPSHOT["genbank"]] if AUGUSTUS_TRAINING_SNAPSHOT else []),
        selection_reference=([AUGUSTUS_TRAINING_SNAPSHOT["selection"]] if AUGUSTUS_TRAINING_SNAPSHOT else []),
        container=PREDICTOR_PREFLIGHT,
        genemark_complete=GENEMARK_COMPLETE,
        hints_complete=ETP_HINTS_COMPLETE,
        augustus_config_complete=f"{PREDICTION_DIR}/03_augustus/config.complete.tsv",
        genome=lambda wildcards: config["genome"],
        training=f"{GENEMARK_RESULT}/training.gtf",
        hints=ETP_HINTS,
        extrinsic_config=AUGUSTUS_EXTRINSIC_CONFIG,
        augustus_config=AUGUSTUS_CONFIG,
    output:
        work=directory(AUGUSTUS_WORK),
        gff3=f"{AUGUSTUS_RESULT}/augustus.hints.gff3",
        gtf=f"{AUGUSTUS_RESULT}/augustus.hints.gtf",
        protein=f"{AUGUSTUS_RESULT}/augustus.hints.aa",
        codingseq=f"{AUGUSTUS_RESULT}/augustus.hints.codingseq",
        trained=directory(f"{AUGUSTUS_RESULT}/trained_species"),
        complete=AUGUSTUS_COMPLETE,
    params:
        assembly_id=ASSEMBLY_ID,
        species=AUGUSTUS_SPECIES,
        rounds=AUGUSTUS_ROUNDS,
        kfold=AUGUSTUS_KFOLD,
        training_flank=AUGUSTUS_TRAINING_FLANK,
        training_selection=AUGUSTUS_TRAINING,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/augustus_etp.log",
    threads: workflow.cores
    conda:
        "../../envs/predictor-host.yaml"
    script:
        "../scripts/run_augustus_etp.py"


rule prediction:
    input:
        code=source_inputs(workflow.basedir, "prediction_complete"),
        genemark=GENEMARK_COMPLETE,
        hints=ETP_HINTS_COMPLETE,
        augustus=AUGUSTUS_COMPLETE,
        masking_rescue=([MASKING_RESCUE_COMPLETE] if MASKING_RESCUE_ENABLED else []),
    output:
        complete=PREDICTION_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/complete.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/prediction_complete.py"


rule run_masking_rescue:
    input:
        code=source_inputs(workflow.basedir, "masking_rescue"),
        augustus_complete=AUGUSTUS_COMPLETE,
        genemark_complete=GENEMARK_COMPLETE,
        genome=lambda wildcards: config["genome"],
        repeats=lambda wildcards: config["repeat_gff"],
        augustus=f"{AUGUSTUS_RESULT}/augustus.hints.gtf",
        genemark=f"{GENEMARK_RESULT}/genemark.gtf",
        transcripts=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
        junctions=JUNCTION_AGGREGATE,
        direct_chains=f"{GENEMARK_RESULT}/selector_evidence/direct_top_chains.gff",
        nonhc_chains=f"{GENEMARK_RESULT}/selector_evidence/nonhc_top_chains.gff",
        trace=f"{GENEMARK_RESULT}/selector_evidence/nonhc.trace",
        hints=ETP_HINTS,
        extrinsic=AUGUSTUS_EXTRINSIC_CONFIG,
        augustus_config=AUGUSTUS_CONFIG,
    output:
        work=directory(f"{MASKING_RESCUE_DIR}/work"),
        gtf=MASKING_RESCUE_GTF,
        targets=f"{MASKING_RESCUE_DIR}/targets.json",
        summary=f"{MASKING_RESCUE_DIR}/summary.json",
        complete=MASKING_RESCUE_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        species=AUGUSTUS_SPECIES,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/prediction/masking_rescue.log",
    threads: workflow.cores
    conda:
        "../../envs/predictor-host.yaml"
    script:
        "../scripts/run_masking_rescue.py"


rule masking_rescue:
    input:
        MASKING_RESCUE_COMPLETE
