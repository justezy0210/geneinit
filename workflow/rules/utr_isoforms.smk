UTR_CODING_SHARDS = [f"{UTR_PREP_DIR}/coding/{shard}.coding.gtf" for shard in UTR_SHARDS]
UTR_RNA_SHARDS = [f"{UTR_PREP_DIR}/rna/{shard}.rna.gtf" for shard in UTR_SHARDS]
UTR_RANKED_ALTERNATIVES = [
    f"{UTR_DIR}/01_ranked/{shard}/alternatives.gff" for shard in UTR_SHARDS
]
UTR_NORMALIZED_GFF3 = [
    f"{UTR_DIR}/02_refined/{shard}/normalized.gff3" for shard in UTR_SHARDS
]
UTR_NORMALIZED_MAPPINGS = [
    f"{UTR_DIR}/02_refined/{shard}/normalized.mapping.tsv" for shard in UTR_SHARDS
]


rule prepare_postqc_utr_shards:
    input:
        code=source_inputs(workflow.basedir, "utr_prepare"),
        gff3=POSTQC_GFF3,
        transcript_decisions=f"{POSTQC_DIR}/{ASSEMBLY_ID}.transcript_decisions.tsv",
        postqc_complete=POSTQC_COMPLETE,
        candidates=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        selection_audit=f"{SELECTION_DIR}/{ASSEMBLY_ID}.selection_audit.tsv",
        merged_rna=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
    output:
        coding=UTR_CODING_SHARDS,
        rna=UTR_RNA_SHARDS,
        mapping=UTR_PREPARATION_MAPPING,
        tasks=UTR_TASKS,
        complete=UTR_PREP_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
        seqids=REFERENCE_SEQUENCE_IDS,
        shards=UTR_SHARDS,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/final_utr/prepare.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/prepare_postqc_utr_shards.py"


rule rank_utr_isoforms:
    input:
        coding=f"{UTR_PREP_DIR}/coding/{{shard}}.coding.gtf",
        rna=f"{UTR_PREP_DIR}/rna/{{shard}}.rna.gtf",
        preparation=UTR_PREP_COMPLETE,
        bams=SAMPLE_BAMS,
        indices=SAMPLE_CSIS,
    output:
        alternatives=f"{UTR_DIR}/01_ranked/{{shard}}/alternatives.gff",
        status=f"{UTR_DIR}/01_ranked/{{shard}}/status.tsv",
    params:
        protocols=[f"{record.layout}:{record.strandedness}" for record in SAMPLE_RECORDS],
        paired=[record.ingenannot_paired for record in SAMPLE_RECORDS],
        stranded=[record.ingenannot_stranded for record in SAMPLE_RECORDS],
        protocol_modes=[record.ingenannot_protocol_mode for record in SAMPLE_RECORDS],
    wildcard_constraints:
        shard="|".join(UTR_SHARDS),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/final_utr/rank/{{shard}}.log",
    threads: 4
    conda:
        "../../envs/utr.yaml"
    script:
        "../scripts/rank_utr_isoforms.py"


rule refine_utr_shard:
    input:
        code=source_inputs(workflow.basedir, "utr_refine"),
        rna=f"{UTR_PREP_DIR}/rna/{{shard}}.rna.gtf",
        coding=f"{UTR_PREP_DIR}/coding/{{shard}}.coding.gtf",
        alternatives=f"{UTR_DIR}/01_ranked/{{shard}}/alternatives.gff",
        rank_status=f"{UTR_DIR}/01_ranked/{{shard}}/status.tsv",
        mapping=UTR_PREPARATION_MAPPING,
    output:
        raw=f"{UTR_DIR}/02_refined/{{shard}}/ingenannot.raw.gff3",
        gff3=f"{UTR_DIR}/02_refined/{{shard}}/normalized.gff3",
        mapping=f"{UTR_DIR}/02_refined/{{shard}}/normalized.mapping.tsv",
        summary=f"{UTR_DIR}/02_refined/{{shard}}/normalized.summary.tsv",
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
    wildcard_constraints:
        shard="|".join(UTR_SHARDS),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/final_utr/refine/{{shard}}.log",
    conda:
        "../../envs/utr.yaml"
    script:
        "../scripts/refine_utr_shard.py"


rule aggregate_utr_shards:
    input:
        code=source_inputs(workflow.basedir, "utr_aggregate"),
        gff3s=UTR_NORMALIZED_GFF3,
        mappings=UTR_NORMALIZED_MAPPINGS,
        preparation=UTR_PREPARATION_MAPPING,
        tasks=UTR_TASKS,
    output:
        gff3=UTR_GFF3,
        mapping=UTR_MAPPING,
        summary=f"{UTR_DIR}/{ASSEMBLY_ID}.utr_models.summary.tsv",
        complete=UTR_COMPLETE,
    params:
        repo_root=lambda wildcards: str(workflow.basedir),
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/final_utr/aggregate.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/aggregate_utr_shards.py"


rule utr_isoforms:
    input:
        UTR_COMPLETE
