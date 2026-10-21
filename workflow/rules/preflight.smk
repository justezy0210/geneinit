import json


rule preflight:
    input:
        genome=lambda wildcards: config["genome"],
        repeat_gff=lambda wildcards: config["repeat_gff"],
        protein_fasta=lambda wildcards: config["protein_fasta"],
        samples=lambda wildcards: config["samples"],
        policy=lambda wildcards: config["policy"],
        te_evidence=lambda wildcards: [config["te_evidence"]] if config.get("te_evidence") else [],
        rna_inputs=_sample_input_files,
        rna_reference=[str(RNA_REFERENCE_FASTA)] if RNA_REFERENCE_FASTA else [],
        prediction_files=lambda wildcards: [
            value for value in (
                PREDICTION_CONFIG.get("augustus_prediction_resume"),
                PREDICTION_CONFIG.get("augustus_optimization_start_parameters"),
                *AUGUSTUS_TRAINING_SNAPSHOT.values(),
            ) if value
        ],
    output:
        normalized_config=f"{PREFLIGHT_DIR}/normalized_config.json",
        normalized_samples=f"{PREFLIGHT_DIR}/normalized_samples.tsv",
        input_manifest=f"{PREFLIGHT_DIR}/input_manifest.tsv",
        summary=f"{PREFLIGHT_DIR}/preflight.json",
        complete=PREFLIGHT_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        config_snapshot=lambda wildcards: json.dumps(config, sort_keys=True),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/preflight.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/preflight.py"
