"""Build AUGUSTUS hints directly from GeneMark-ETP evidence."""

from pathlib import Path
import subprocess
import sys
import tempfile
import traceback


sys.path.insert(0, str(Path(snakemake.params.repo_root)))

from geneinit.contracts.deployment import (  # noqa: E402
    apptainer_command,
    load_predictor_deployment,
)
from geneinit.predictor.common import (  # noqa: E402
    PredictorError,
    fasta_lengths,
    require_nonempty,
    sha256_file,
    validate_gff_coordinates,
    write_completion,
)


def feature_rows(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.open(encoding="utf-8") if line.strip()]


def source_count(rows: list[str], wanted: str) -> int:
    count = 0
    for row in rows:
        attributes = row.split("\t", 8)[8].split(";")
        if any(value.strip() in {f"src={wanted}", f"source={wanted}"} for value in attributes):
            count += 1
    return count


def fields(row: str) -> list[str]:
    values = row.split("\t")
    if len(values) != 9:
        raise PredictorError(f"hint row does not have 9 columns: {row[:120]}")
    return values


log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
try:
    deployment = load_predictor_deployment(dict(snakemake.config))
    work = Path(snakemake.input.work).resolve()
    genome = Path(snakemake.input.genome).resolve()
    proc = str(snakemake.params.proc)
    rna_hints = require_nonempty(work / "rnaseq/hints/hintsfile_merged.gff", "RNA hints")
    hc_hints = require_nonempty(
        work / f"rnaseq/hints/{proc}/prothint/prothint_augustus.gff",
        "HC protein hints",
    )
    nonhc_hints = require_nonempty(
        work / f"{proc}/nonhc/prothint/prothint_augustus.gff",
        "non-HC protein hints",
    )
    nonhc_trace = require_nonempty(
        work / f"{proc}/nonhc/nonhc.trace", "non-HC trace"
    )

    output = Path(snakemake.output.hints).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".geneinit_hints_", dir=output.parent) as raw_tmp:
        temporary = Path(raw_tmp)
        converted_log = temporary / "format_back.stderr.log"
        command = apptainer_command(
            deployment,
            (
                "/usr/bin/perl",
                "/opt/ETP/bin/format_back.pl",
                str(nonhc_hints),
                str(nonhc_trace),
            ),
            (work,),
            pwd=work,
        )
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        converted_log.write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise PredictorError(f"GeneMark format_back.pl failed; inspect {converted_log}")
        converted_rows = [row for row in completed.stdout.splitlines() if row.strip()]
        if not converted_rows:
            raise PredictorError("GeneMark format_back.pl produced no non-HC hints")
        rna_rows = feature_rows(rna_hints)
        hc_rows = feature_rows(hc_hints)
        raw_rows = converted_rows + rna_rows + hc_rows + rna_rows
        if len(raw_rows) != len(converted_rows) + len(rna_rows) * 2 + len(hc_rows):
            raise PredictorError("ETP hint assembly lost an evidence row")
        if source_count(raw_rows, "E") != len(rna_rows) * 2:
            raise PredictorError("full hints do not contain exactly two RNA-hint copies")

        no_merge: list[str] = []
        to_merge: list[str] = []
        for row in raw_rows:
            attributes = fields(row)[8]
            if "src=C" in row and ("grp=" in attributes or "group=" in attributes):
                no_merge.append(row)
            else:
                to_merge.append(row)
        to_merge.sort(
            key=lambda row: (
                fields(row)[0],
                fields(row)[2],
                int(fields(row)[4]),
                int(fields(row)[3]),
            )
        )
        merge_input = temporary / "hints.to_merge.sorted.gff"
        merge_input.write_text("\n".join(to_merge) + "\n", encoding="utf-8")
        join_command = apptainer_command(
            deployment,
            ("/usr/bin/perl", "/opt/Augustus/scripts/join_mult_hints.pl"),
            (temporary,),
        )
        with merge_input.open("r", encoding="utf-8") as stdin:
            joined = subprocess.run(
                join_command,
                stdin=stdin,
                text=True,
                capture_output=True,
                check=False,
            )
        if joined.returncode != 0:
            raise PredictorError(f"join_mult_hints.pl failed: {joined.stderr.strip()}")
        joined_rows = [row for row in joined.stdout.splitlines() if row.strip()]
        final_rows = joined_rows + no_merge
        final_rows.sort(
            key=lambda row: (
                fields(row)[0],
                int(fields(row)[3]),
                int(fields(row)[4]),
                fields(row)[2],
                fields(row)[6],
                fields(row)[7],
                fields(row)[8],
            )
        )
        candidate = temporary / "candidate.gff"
        candidate.write_text("\n".join(final_rows) + "\n", encoding="utf-8")
        validate_gff_coordinates(candidate, fasta_lengths(genome))
        candidate.replace(output)

    write_completion(
        Path(snakemake.output.complete),
        (
            ("schema_version", 1),
            ("status", "PASS"),
            ("method", "GeneMark format_back + RNA/HC assembly + AUGUSTUS join_mult_hints + deterministic sort"),
            ("rna_copy_factor", 2),
            ("converted_nonhc_rows", len(converted_rows)),
            ("raw_hint_rows", len(raw_rows)),
            ("joined_hint_rows", len(final_rows)),
            ("src_E_rows", source_count(final_rows, "E")),
            ("src_C_rows", source_count(final_rows, "C")),
            ("src_P_rows", source_count(final_rows, "P")),
            ("src_M_rows", source_count(final_rows, "M")),
            ("container_sha256", deployment.sif_sha256),
            ("hints_sha256", sha256_file(output)),
        ),
    )
except Exception:
    log.write_text(traceback.format_exc(), encoding="utf-8")
    raise
else:
    log.write_text("full ETP hints PASS\n", encoding="utf-8")
