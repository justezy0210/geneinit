JUNCTION_TABLES = [
    f"{JUNCTION_DIR}/samples/{SAMPLE_BY_ID[sample].layout}/{sample}/"
    f"{sample}.junctions.tsv.gz"
    for sample in SAMPLES
]
JUNCTION_INDICES = [f"{path}.csi" for path in JUNCTION_TABLES]


rule extract_sample_junctions:
    input:
        bam=_sample_bam,
        csi=_sample_csi,
    output:
        table=f"{JUNCTION_DIR}/samples/{{layout}}/{{sample}}/{{sample}}.junctions.tsv.gz",
        index=f"{JUNCTION_DIR}/samples/{{layout}}/{{sample}}/{{sample}}.junctions.tsv.gz.csi",
        summary=f"{JUNCTION_DIR}/samples/{{layout}}/{{sample}}/{{sample}}.junctions.summary.tsv",
    wildcard_constraints:
        layout="PE|SE",
        sample="|".join(map(re.escape, SAMPLES)),
    params:
        repo_root=str(workflow.basedir),
        sample=lambda wildcards: _sample_record(wildcards).sample,
        condition=lambda wildcards: _sample_record(wildcards).condition,
        layout=lambda wildcards: _sample_record(wildcards).layout,
        strandedness=lambda wildcards: _sample_record(wildcards).strandedness,
        min_mapq=JUNCTION_MIN_MAPQ,
        min_anchor=JUNCTION_MIN_ANCHOR,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/evidence/junctions/{{layout}}/{{sample}}.log",
    threads: 2
    conda:
        "../../envs/evidence.yaml"
    script:
        "../scripts/extract_junctions.py"


rule aggregate_junctions:
    input:
        genome=lambda wildcards: config["genome"],
        tables=JUNCTION_TABLES,
        indices=JUNCTION_INDICES,
    output:
        table=JUNCTION_AGGREGATE,
        index=f"{JUNCTION_AGGREGATE}.csi",
        summary=f"{JUNCTION_DIR}/aggregate/{ASSEMBLY_ID}.hq_junctions.summary.tsv",
        complete=JUNCTION_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
        samples=SAMPLES,
        conditions=[SAMPLE_BY_ID[sample].condition for sample in SAMPLES],
        min_mapq=JUNCTION_MIN_MAPQ,
        min_anchor=JUNCTION_MIN_ANCHOR,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/evidence/junctions/aggregate.log",
    conda:
        "../../envs/evidence.yaml"
    script:
        "../scripts/aggregate_junctions.py"


rule junction_evidence:
    input:
        JUNCTION_COMPLETE


rule score_component_evidence:
    input:
        candidates=f"{CANDIDATE_DIR}/{ASSEMBLY_ID}.candidate_universe.tsv",
        candidate_complete=CANDIDATE_COMPLETE,
        genome=lambda wildcards: config["genome"],
        merged_gtf=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
        junctions=JUNCTION_AGGREGATE,
        junction_index=f"{JUNCTION_AGGREGATE}.csi",
        junction_complete=JUNCTION_COMPLETE,
        direct_chains=f"{GENEMARK_RESULT}/selector_evidence/direct_top_chains.gff",
        direct_hits=f"{GENEMARK_RESULT}/selector_evidence/direct_protein_hits.tsv",
        nonhc_chains=f"{GENEMARK_RESULT}/selector_evidence/nonhc_top_chains.gff",
        nonhc_hits=f"{GENEMARK_RESULT}/selector_evidence/nonhc_protein_hits.tsv",
        nonhc_trace=f"{GENEMARK_RESULT}/selector_evidence/nonhc.trace",
        repeats=lambda wildcards: config["repeat_gff"],
    output:
        table=COMPONENT_EVIDENCE,
        summary=f"{COMPONENT_DIR}/{ASSEMBLY_ID}.component_evidence.summary.tsv",
        complete=COMPONENT_COMPLETE,
    params:
        repo_root=str(workflow.basedir),
        assembly_id=ASSEMBLY_ID,
        expected_protein_copies=int(config["expected_protein_copies"]),
        repeat_scope="classified_te_v1",
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/evidence/components.log",
    conda:
        "../../envs/preflight.yaml"
    script:
        "../scripts/score_component_evidence.py"


rule component_evidence:
    input:
        COMPONENT_COMPLETE
