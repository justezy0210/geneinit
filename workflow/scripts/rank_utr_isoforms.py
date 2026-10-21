from pathlib import Path
import os
import shutil
import subprocess


def contract_bool(value):
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"invalid InGenAnnot boolean contract: {value!r}")
    return normalized == "true"


protocols = [str(value) for value in snakemake.params.protocols]
paired_flags = [contract_bool(value) for value in snakemake.params.paired]
stranded_flags = [contract_bool(value) for value in snakemake.params.stranded]
protocol_modes = [str(value) for value in snakemake.params.protocol_modes]
if not (
    len(snakemake.input.bams)
    == len(protocols)
    == len(paired_flags)
    == len(stranded_flags)
    == len(protocol_modes)
):
    raise ValueError("InGenAnnot BAM protocol metadata length mismatch")

rna_gtf = Path(snakemake.input.rna)
coding_gtf = Path(snakemake.input.coding)
alternatives = Path(snakemake.output.alternatives)
status = Path(snakemake.output.status)
log = Path(snakemake.log[0])
alternatives.parent.mkdir(parents=True, exist_ok=True)
status.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)


def feature_count(path: Path, feature_types: set[str]) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(
            fields[2] in feature_types
            for raw in handle
            if raw.strip() and not raw.startswith("#")
            for fields in [raw.rstrip("\n").split("\t")]
            if len(fields) == 9
        )


coding_count = feature_count(coding_gtf, {"transcript", "mRNA"})
rna_count = feature_count(rna_gtf, {"transcript", "mRNA"})
if not coding_count or not rna_count:
    alternatives.write_text("", encoding="utf-8")
    state = "NO_CODING" if not coding_count else "NO_RNA"
    log.write_text(f"InGenAnnot ranking skipped: {state}\n", encoding="utf-8")
else:
    fof = alternatives.parent / "bam.fof"
    with fof.open("w", encoding="utf-8") as handle:
        for bam, paired, stranded in zip(
            snakemake.input.bams, paired_flags, stranded_flags, strict=True
        ):
            # Keep an explicitly supplied BAM path instead of resolving a
            # symlink: a deployment may intentionally pair that path with a
            # coordinate-matched local .csi/.bai view.
            handle.write(
                f"{Path(bam).absolute()}\t{paired}\t{stranded}\n"
            )
    prefix = alternatives.parent / "ranked_isoforms"
    environment = dict(os.environ)
    environment["PYTHONNOUSERSITE"] = "1"
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[variable] = "1"
    mpl_dir = alternatives.parent / ".mplconfig"
    mpl_dir.mkdir(exist_ok=True)
    environment["MPLCONFIGDIR"] = str(mpl_dir.resolve())
    command = [
        "ingenannot", "-v", "2", "-p", str(snakemake.threads),
        "isoform_ranking", str(rna_gtf), "-f", str(fof),
        "--prefix", str(prefix),
    ]
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=environment, check=True)
    produced = Path(f"{prefix}.alternatives.gff")
    if produced.is_file() and produced.stat().st_size:
        shutil.move(str(produced), alternatives)
        state = "HAS_RANKED_ALTERNATIVES"
    else:
        alternatives.write_text("", encoding="utf-8")
        state = "NO_RANKED_ALTERNATIVE"

status.write_text(
    "metric\tvalue\n"
    f"status\t{state}\n"
    f"coding_models\t{coding_count}\n"
    f"rna_transcripts\t{rna_count}\n"
    f"protocols\t{','.join(protocols)}\n"
    f"protocol_modes\t{','.join(protocol_modes)}\n",
    encoding="utf-8",
)
