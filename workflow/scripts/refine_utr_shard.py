import csv
from pathlib import Path
import os
import subprocess
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.utr.core import (  # noqa: E402
    POSTQC_PREPARATION_COLUMNS,
    UTR_MAPPING_COLUMNS,
    normalize_refinement,
    restore_postqc_gene_hierarchy,
    sha256_file,
    write_normalized_shard,
)


from geneinit.utr.ranked import validate_ranked_rna  # noqa: E402


coding = Path(snakemake.input.coding)
alternatives = Path(snakemake.input.alternatives)
preparation = Path(snakemake.input.mapping)
output_gff3 = Path(snakemake.output.gff3)
output_mapping = Path(snakemake.output.mapping)
summary = Path(snakemake.output.summary)
raw_gff3 = Path(snakemake.output.raw)
log = Path(snakemake.log[0])
for path in (output_gff3, output_mapping, summary, raw_gff3, log):
    path.parent.mkdir(parents=True, exist_ok=True)

with preparation.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = set(POSTQC_PREPARATION_COLUMNS) - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"preparation mapping lacks columns: {sorted(missing)}")
    preparation_rows = {row["candidate_id"]: row for row in reader}


def coding_model_count() -> int:
    with coding.open(encoding="utf-8") as handle:
        return sum(
            fields[2] in {"transcript", "mRNA"}
            for raw in handle
            if raw.strip() and not raw.startswith("#")
            for fields in [raw.rstrip("\n").split("\t")]
            if len(fields) == 9
        )


rank_validation = {"exon_strands_restored": 0, "unstranded_transcripts_excluded": 0}
count = coding_model_count()
if count == 0:
    output_gff3.write_text("##gff-version 3\n", encoding="utf-8")
    with output_mapping.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(UTR_MAPPING_COLUMNS)
    raw_gff3.write_text("", encoding="utf-8")
    log.write_text("InGenAnnot UTR refinement skipped: NO_CODING\n", encoding="utf-8")
    updated = 0
else:
    raw_input = None
    validated = alternatives
    if alternatives.stat().st_size:
        validated = raw_gff3.parent / "alternatives.validated.gtf"
        rank_validation = validate_ranked_rna(
            alternatives, Path(snakemake.input.rna), validated,
            raw_gff3.parent / "ranked_strand_repairs.tsv",
        )
    if validated.stat().st_size:
        environment = dict(os.environ)
        environment["PYTHONNOUSERSITE"] = "1"
        mpl_dir = raw_gff3.parent / ".mplconfig"
        mpl_dir.mkdir(exist_ok=True)
        environment["MPLCONFIGDIR"] = str(mpl_dir.resolve())
        command = [
            "ingenannot", "-v", "2", "utr_refine", str(coding),
            str(validated), str(raw_gff3), "--utr_mode", "rank",
        ]
        with log.open("w", encoding="utf-8") as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=environment, check=True)
        raw_input = raw_gff3
    else:
        raw_gff3.write_text("", encoding="utf-8")
        log.write_text("InGenAnnot UTR refinement skipped: NO_STRAND_VALIDATED_ALTERNATIVE\n", encoding="utf-8")
    order, models = normalize_refinement(coding, raw_input)
    restore_postqc_gene_hierarchy(order, models, preparation_rows)
    rows = write_normalized_shard(output_gff3, output_mapping, order, models, preparation_rows)
    updated = sum(row["utr_refine_evidence"] not in {"", ".", "NA", "NONE"} for row in rows)

summary.write_text(
    "metric\tvalue\n"
    "status\tPASS_MIGRATION_PENDING\n"
    f"coding_model_count\t{count}\n"
    f"ranked_exon_strands_restored_from_source\t{rank_validation['exon_strands_restored']}\n"
    f"unknown_strand_rna_transcripts_excluded\t{rank_validation['unstranded_transcripts_excluded']}\n"
    f"utr_updated_model_count\t{updated}\n"
    f"coding_gtf_sha256\t{sha256_file(coding)}\n"
    f"normalized_gff3_sha256\t{sha256_file(output_gff3)}\n"
    f"mapping_sha256\t{sha256_file(output_mapping)}\n",
    encoding="utf-8",
)
