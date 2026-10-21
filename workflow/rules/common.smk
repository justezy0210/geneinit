from pathlib import Path
from collections import Counter
import re

from geneinit.contracts.samples import load_samples
from geneinit.contracts.workflow import source_inputs
from geneinit.contracts.preflight import input_preflight_config
from geneinit.contracts.alignments import check_staging_paths, rna_input_config


def _safe_runtime_name(value, field):
    value = str(value)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError(f"{field} must match [A-Za-z0-9_.-]+: {value!r}")
    return value


def _sample_input_files(_wildcards):
    paths = []
    for record in SAMPLE_RECORDS:
        inputs = (record.bam, record.csi) if RNA_INPUT_MODE == "bam" else (record.read1, record.read2)
        paths.extend(str(path) for path in inputs if path is not None)
    return paths


def _reference_sequence_ids(path):
    sequence_ids = []
    seen = set()
    with Path(path).open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            seqid = raw[1:].strip().split(maxsplit=1)[0]
            if not seqid or seqid in seen:
                raise ValueError(f"duplicate/empty reference sequence ID at line {number}: {seqid!r}")
            seen.add(seqid)
            sequence_ids.append(seqid)
    if not sequence_ids:
        raise ValueError(f"reference FASTA has no sequences: {path}")
    return sequence_ids


EXECUTION_ID = _safe_runtime_name(config.get("execution_id", "development"), "execution_id")
EXECUTION_ROOT = str(config.get("execution_root", "runs")).rstrip("/")
if not EXECUTION_ROOT or EXECUTION_ROOT == "/":
    raise ValueError("execution_root must be a specific directory, not filesystem root")

PREFLIGHT_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/00_preflight"
PREFLIGHT_COMPLETE = f"{PREFLIGHT_DIR}/preflight.complete.tsv"
INPUT_PREFLIGHT_DIR = f"{PREFLIGHT_DIR}/inputs"
INPUT_PREFLIGHT_COMPLETE = f"{INPUT_PREFLIGHT_DIR}/preflight.complete.tsv"
INPUT_VALIDATION = f"{INPUT_PREFLIGHT_DIR}/input_validation.json"
REFERENCE_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/01_reference"
READS_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/02_reads"
ALIGNMENT_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/03_alignment"
TRANSCRIPTS_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/04_transcripts"
PREDICTION_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/05_prediction"
EVIDENCE_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/06_evidence"
ASSEMBLY_ID = str(config["assembly_id"])
INDEX_PREFIX = f"{REFERENCE_DIR}/hisat2/{ASSEMBLY_ID}"
HISAT2_INDEX = [f"{INDEX_PREFIX}.{part}.ht2" for part in range(1, 9)]

RNA_INPUT_MODE, RNA_REFERENCE_FASTA = rna_input_config(config)
SAMPLE_RECORDS = load_samples(Path(config["samples"]), input_mode=RNA_INPUT_MODE)
SAMPLE_BY_ID = {record.sample: record for record in SAMPLE_RECORDS}
SAMPLES = sorted(SAMPLE_BY_ID)
PE_SAMPLES = sorted(record.sample for record in SAMPLE_RECORDS if record.layout == "PE")
SE_SAMPLES = sorted(record.sample for record in SAMPLE_RECORDS if record.layout == "SE")
SAMPLE_BAMS = [
    f"{ALIGNMENT_DIR}/{SAMPLE_BY_ID[sample].layout}/{sample}/"
    f"{ASSEMBLY_ID}.{sample}.sorted.bam"
    for sample in SAMPLES
]
SAMPLE_CSIS = [f"{bam}.csi" for bam in SAMPLE_BAMS]
RNA_IMPORT_AUDITS = [f"{bam}.import.json" for bam in SAMPLE_BAMS] if RNA_INPUT_MODE == "bam" else []
if RNA_INPUT_MODE == "bam":
    imported_sources = [path for sample in SAMPLE_RECORDS for path in (sample.bam, sample.csi)]
    for sample, bam, csi in zip(SAMPLE_RECORDS, SAMPLE_BAMS, SAMPLE_CSIS):
        check_staging_paths(sample, Path(bam), Path(csi), Path(f"{bam}.import.json"),
                            other_sources=imported_sources)
CONDITION_SAMPLE_COUNTS = Counter(record.condition for record in SAMPLE_RECORDS)
REPLICATES_PER_CONDITION = min(CONDITION_SAMPLE_COUNTS.values())

PREDICTION_CONFIG = config.get("prediction", {})
MASKING_RESCUE_ENABLED = bool(PREDICTION_CONFIG.get("masking_rescue", False))
MASKING_RESCUE_DIR = f"{PREDICTION_DIR}/04_masking_rescue"
MASKING_RESCUE_GTF = f"{MASKING_RESCUE_DIR}/augustus.masking_rescue.gtf"
MASKING_RESCUE_COMPLETE = f"{MASKING_RESCUE_DIR}/masking_rescue.complete.tsv"
GENEMARK_PROC = str(PREDICTION_CONFIG.get("genemark_proc", "protein_evidence"))
AUGUSTUS_SPECIES = str(
    PREDICTION_CONFIG.get("augustus_species", f"{ASSEMBLY_ID}_etp")
)
AUGUSTUS_PREDICTION_RESUME = PREDICTION_CONFIG.get("augustus_prediction_resume")
AUGUSTUS_ROUNDS = int(PREDICTION_CONFIG.get("augustus_rounds", 5))
AUGUSTUS_KFOLD = int(PREDICTION_CONFIG.get("augustus_kfold", 8))
AUGUSTUS_OPTIMIZATION_START = PREDICTION_CONFIG.get("augustus_optimization_start_parameters")
AUGUSTUS_TRAINING_FLANK = int(PREDICTION_CONFIG.get("augustus_training_flank", 1000))
AUGUSTUS_TRAINING = dict(PREDICTION_CONFIG.get("augustus_training", {}))
AUGUSTUS_TRAINING_SNAPSHOT = dict(PREDICTION_CONFIG.get("augustus_training_snapshot", {}))
AUGUSTUS_EXTRINSIC_CONFIG = str(
    Path(workflow.basedir) / "config/augustus/geneinit_etp.cfg"
)

EVIDENCE_CONFIG = config.get("evidence", {})
JUNCTION_MIN_MAPQ = int(EVIDENCE_CONFIG.get("junction_min_mapq", 10))
JUNCTION_MIN_ANCHOR = int(EVIDENCE_CONFIG.get("junction_min_anchor", 8))
JUNCTION_DIR = f"{EVIDENCE_DIR}/01_junctions"
JUNCTION_AGGREGATE = f"{JUNCTION_DIR}/aggregate/{ASSEMBLY_ID}.hq_junctions.tsv.gz"
JUNCTION_COMPLETE = f"{JUNCTION_DIR}/aggregate/junctions.complete.tsv"
CANDIDATE_DIR = f"{EVIDENCE_DIR}/02_candidates"
CANDIDATE_COMPLETE = f"{CANDIDATE_DIR}/candidate_universe.complete.tsv"
COMPONENT_DIR = f"{EVIDENCE_DIR}/03_components"
COMPONENT_EVIDENCE = f"{COMPONENT_DIR}/{ASSEMBLY_ID}.component_evidence.tsv"
COMPONENT_COMPLETE = f"{COMPONENT_DIR}/component_evidence.complete.tsv"
SELECTION_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/07_selection"
SELECTED_KEYS = f"{SELECTION_DIR}/{ASSEMBLY_ID}.selected_candidate_ids.txt"
SELECTION_COMPLETE = f"{SELECTION_DIR}/selection.complete.tsv"

REFERENCE_SEQUENCE_IDS = _reference_sequence_ids(config["genome"])
UTR_SHARDS = [f"{index:06d}" for index in range(1, len(REFERENCE_SEQUENCE_IDS) + 1)]
UTR_SHARD_TO_SEQID = dict(zip(UTR_SHARDS, REFERENCE_SEQUENCE_IDS, strict=True))
CODING_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/08_coding_models"
CODING_GFF3 = f"{CODING_DIR}/{ASSEMBLY_ID}.selected_coding.gff3"
CODING_MAPPING = f"{CODING_DIR}/{ASSEMBLY_ID}.selected_coding.mapping.tsv"
CODING_COMPLETE = f"{CODING_DIR}/selected_coding.complete.tsv"
GROUPING_DIR = f"{CODING_DIR}/01_grouped"
GROUPED_GFF3 = f"{GROUPING_DIR}/{ASSEMBLY_ID}.provisional_genes.gff3"
GROUPING_MAPPING = f"{GROUPING_DIR}/{ASSEMBLY_ID}.provisional_genes.mapping.tsv"
GROUPING_COMPLETE = f"{GROUPING_DIR}/gene_grouping.complete.tsv"
POSTQC_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/09_iterative_postqc"
POSTQC_GFF3 = f"{POSTQC_DIR}/{ASSEMBLY_ID}.postqc.coding.gff3"
POSTQC_COMPLETE = f"{POSTQC_DIR}/postqc.complete.tsv"
UTR_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/10_final_utr"
UTR_PREP_DIR = f"{UTR_DIR}/00_prepared"
UTR_PREPARATION_MAPPING = f"{UTR_PREP_DIR}/{ASSEMBLY_ID}.postqc_coding.mapping.tsv"
UTR_TASKS = f"{UTR_PREP_DIR}/{ASSEMBLY_ID}.utr_tasks.tsv"
UTR_PREP_COMPLETE = f"{UTR_PREP_DIR}/preparation.complete.tsv"
UTR_GFF3 = f"{UTR_DIR}/{ASSEMBLY_ID}.final_utr_models.gff3"
UTR_MAPPING = f"{UTR_DIR}/{ASSEMBLY_ID}.final_utr_models.mapping.tsv"
UTR_COMPLETE = f"{UTR_DIR}/utr_isoforms.complete.tsv"
RELEASE_STAGING_DIR = f"{EXECUTION_ROOT}/{EXECUTION_ID}/11_release.staging/{ASSEMBLY_ID}"
RELEASE_COMPLETE = f"{RELEASE_STAGING_DIR}/release.complete.tsv"

PREDICTOR_PREFLIGHT = f"{PREDICTION_DIR}/00_container/preflight.complete.tsv"
GENEMARK_WORK = f"{PREDICTION_DIR}/01_genemark/work"
GENEMARK_RESULT = f"{PREDICTION_DIR}/01_genemark/results"
GENEMARK_COMPLETE = f"{GENEMARK_RESULT}/genemark.complete.tsv"
ETP_HINTS = f"{PREDICTION_DIR}/02_hints/{ASSEMBLY_ID}.etp_hints.gff"
ETP_HINTS_COMPLETE = f"{PREDICTION_DIR}/02_hints/etp_hints.complete.tsv"
AUGUSTUS_CONFIG = f"{PREDICTION_DIR}/03_augustus/config"
AUGUSTUS_WORK = f"{PREDICTION_DIR}/03_augustus/work"
AUGUSTUS_RESULT = f"{PREDICTION_DIR}/03_augustus/results"
AUGUSTUS_COMPLETE = f"{AUGUSTUS_RESULT}/augustus.complete.tsv"
PREDICTION_COMPLETE = f"{PREDICTION_DIR}/prediction.complete.tsv"


def _deployment_value(key):
    deployment = config.get("deployment")
    if not isinstance(deployment, dict) or key not in deployment:
        raise ValueError(
            f"target requires config deployment.{key}; see config/config.example.yaml"
        )
    return deployment[key]


def _sample_record(wildcards):
    record = SAMPLE_BY_ID[wildcards.sample]
    if hasattr(wildcards, "layout") and wildcards.layout != record.layout:
        raise ValueError(
            f"sample {record.sample} has layout {record.layout}, not {wildcards.layout}"
        )
    return record


def _clean_read1(wildcards):
    record = _sample_record(wildcards)
    if record.layout == "PE":
        return f"{READS_DIR}/PE/{record.sample}/{record.sample}.clean.R1.fastq.gz"
    return f"{READS_DIR}/SE/{record.sample}/{record.sample}.clean.fastq.gz"


def _sample_bam(wildcards):
    record = _sample_record(wildcards)
    return (
        f"{ALIGNMENT_DIR}/{record.layout}/{record.sample}/"
        f"{ASSEMBLY_ID}.{record.sample}.sorted.bam"
    )


def _sample_csi(wildcards):
    return _sample_bam(wildcards) + ".csi"


def _hisat2_strand_option(wildcards):
    value = _sample_record(wildcards).hisat2_strandness
    return f"--rna-strandness {value}" if value else ""


def _stringtie_strand_option(wildcards):
    return _sample_record(wildcards).stringtie_flag or ""
