"""Stage evidence and run GeneMark-ETP without its unstranded StringTie fallback."""

from pathlib import Path
import os
import shutil
import sys
import traceback

import yaml


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.deployment import (  # noqa: E402
    apptainer_command,
    load_predictor_deployment,
)
from geneinit.predictor.common import (  # noqa: E402
    PredictorError,
    atomic_copy,
    require_nonempty,
    run_logged,
    sha256_file,
    write_completion,
)
from geneinit.predictor.commands import genemark_etp_command  # noqa: E402


def symlink_checked(source: Path, target: Path) -> None:
    source = require_nonempty(source, "BAM or CSI").resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.symlink_to(source)


log = Path(snakemake.log[0])
try:
    deployment = load_predictor_deployment(dict(snakemake.config))
    genome = Path(snakemake.input.genome).resolve()
    proteins = Path(snakemake.input.proteins).resolve()
    merged_gtf = Path(snakemake.input.merged_gtf).resolve()
    bams = [Path(value).resolve() for value in snakemake.input.bams]
    csis = [Path(value).resolve() for value in snakemake.input.csis]
    sample_ids = list(snakemake.params.samples)
    if not (len(sample_ids) == len(bams) == len(csis)):
        raise PredictorError("sample, BAM, and CSI counts differ")

    work = Path(snakemake.output.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    staged_gtf = work / "rnaseq/stringtie/transcripts_merged.gff"
    atomic_copy(merged_gtf, staged_gtf)
    if staged_gtf.is_symlink():
        raise PredictorError("staged StringTie GTF must be a regular file")

    bam_dir = work / "input_bam"
    bam_dir.mkdir(parents=True, exist_ok=True)
    for sample, bam, csi in zip(sample_ids, bams, csis, strict=True):
        symlink_checked(bam, bam_dir / f"{sample}.bam")
        symlink_checked(csi, bam_dir / f"{sample}.bam.csi")

    cfg = work / "geneinit.gmetp.yaml"
    cfg.write_text(
        yaml.safe_dump(
            {
                "species": str(snakemake.params.species),
                "genome_path": str(genome),
                "rnaseq_sets": sample_ids,
                "protdb_path": str(proteins),
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    inner_command = genemark_etp_command(
        config=cfg,
        workdir=work,
        bam_dir=bam_dir,
        proc=str(snakemake.params.proc),
        threads=snakemake.threads,
    )
    command = apptainer_command(
        deployment,
        inner_command,
        (work, genome, proteins, *bams, *csis),
    )
    if "--force" in command:
        raise PredictorError("GeneMark-ETP --force is forbidden")
    run_logged(command, log)

    proc = str(snakemake.params.proc)
    internal = {
        "genemark": work / "genemark.gtf",
        "supported": work / "genemark_supported.gtf",
        "rna_hints": work / "rnaseq/hints/hintsfile_merged.gff",
        "protein_hints": work / f"rnaseq/hints/{proc}/prothint/prothint.gff",
        "protein_evidence": work / f"rnaseq/hints/{proc}/prothint/evidence.gff",
        "direct_top_chains": work / f"rnaseq/hints/{proc}/prothint/top_chains.gff",
        "direct_protein_hits": work / f"rnaseq/hints/{proc}/prothint/diamond/diamond.out",
        "nonhc_top_chains": work / f"{proc}/nonhc/prothint/top_chains.gff",
        "nonhc_protein_hits": work / f"{proc}/nonhc/prothint/diamond/diamond.out",
        "nonhc_trace": work / f"{proc}/nonhc/nonhc.trace",
        "training": work / f"{proc}/model/training.gtf",
        "model": work / f"{proc}/model/output.mod",
    }
    for label, source in internal.items():
        require_nonempty(source, label)
        atomic_copy(source, Path(getattr(snakemake.output, label)))

    write_completion(
        Path(snakemake.output.complete),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("assembly_id", snakemake.params.assembly_id),
            ("container_sha256", deployment.sif_sha256),
            ("config_sha256", sha256_file(cfg)),
            ("staged_stringtie_sha256", sha256_file(staged_gtf)),
            ("genemark_gtf_sha256", sha256_file(Path(snakemake.output.genemark))),
            ("genemark_supported_gtf_sha256", sha256_file(Path(snakemake.output.supported))),
            ("direct_top_chains_sha256", sha256_file(Path(snakemake.output.direct_top_chains))),
            ("direct_protein_hits_sha256", sha256_file(Path(snakemake.output.direct_protein_hits))),
            ("nonhc_top_chains_sha256", sha256_file(Path(snakemake.output.nonhc_top_chains))),
            ("nonhc_protein_hits_sha256", sha256_file(Path(snakemake.output.nonhc_protein_hits))),
            ("nonhc_trace_sha256", sha256_file(Path(snakemake.output.nonhc_trace))),
            ("force_used", "false"),
        ),
    )
except Exception:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("\n" + traceback.format_exc())
    raise
