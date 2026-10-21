"""Versioned, taxon-neutral predictor command contracts."""

from __future__ import annotations

from pathlib import Path


def genemark_etp_command(
    *,
    config: Path,
    workdir: Path,
    bam_dir: Path,
    proc: str,
    threads: int,
) -> tuple[str, ...]:
    """Return GeneMark-ETP arguments; no taxon mode is assumed.

    GeneMark-ETP launches Perl helpers through ``#!/usr/bin/env perl``. Keep
    the image's existing PATH for external tools, but put the system Perl
    first because that interpreter owns GeneMark-ETP's CPAN dependencies.
    """

    return (
        "/bin/bash",
        "-c",
        'export PATH="/usr/bin:/bin:${PATH}"; exec /usr/bin/perl "$@"',
        "geneinit-gmetp",
        "/opt/ETP/bin/gmetp.pl",
        "--cfg",
        str(config),
        "--workdir",
        str(workdir),
        "--bam",
        str(bam_dir),
        "--proc",
        proc,
        "--softmask",
        "--gc",
        "0",
        "--extend",
        "1",
        "--cores",
        str(threads),
        "--clean",
        "--verbose",
    )


def augustus_new_species_command(
    *,
    species: str,
    augustus_config: Path,
) -> tuple[str, ...]:
    """Create a new eukaryotic AUGUSTUS species parameter set."""

    return (
        "/usr/bin/perl",
        "/opt/Augustus/scripts/new_species.pl",
        f"--species={species}",
        f"--AUGUSTUS_CONFIG_PATH={augustus_config}",
        "--silent",
    )


def augustus_training_genbank_command(
    *, training_gtf: Path, genome: Path, flank: int, output: Path
) -> tuple[str, ...]:
    """Convert curated training genes to AUGUSTUS' compact GenBank format."""

    return (
        "/usr/bin/perl",
        "/opt/Augustus/scripts/gff2gbSmallDNA.pl",
        str(training_gtf),
        str(genome),
        str(flank),
        str(output),
    )


def augustus_etraining_command(
    *, species: str, training_genbank: Path, augustus_config: Path
) -> tuple[str, ...]:
    """Train AUGUSTUS parameters for one species."""

    return (
        "/opt/Augustus/bin/etraining",
        f"--species={species}",
        f"--AUGUSTUS_CONFIG_PATH={augustus_config}",
        "--stopCodonExcludedFromCDS=false",
        str(training_genbank),
    )


def augustus_optimize_command(
    *,
    species: str,
    training_genbank: Path,
    augustus_config: Path,
    threads: int,
    rounds: int,
    kfold: int,
    onlytrain: Path | None = None,
) -> tuple[str, ...]:
    """Optimize AUGUSTUS meta-parameters directly."""

    return (
        "/usr/bin/perl",
        "/opt/Augustus/scripts/optimize_augustus.pl",
        f"--species={species}",
        f"--AUGUSTUS_CONFIG_PATH={augustus_config}",
        f"--rounds={rounds}",
        f"--cpus={threads}",
        f"--kfold={kfold}",
        "--aug_exec_dir=/opt/Augustus/bin/",
        *((f"--onlytrain={onlytrain}",) if onlytrain is not None else ()),
        str(training_genbank),
    )


def augustus_prediction_command(
    *,
    species: str,
    genome: Path,
    hints: Path,
    extrinsic_config: Path,
    augustus_config: Path,
    output: Path,
) -> tuple[str, ...]:
    """Predict genes with trained AUGUSTUS and ETP hints directly."""

    return (
        "/opt/Augustus/bin/augustus",
        f"--species={species}",
        f"--AUGUSTUS_CONFIG_PATH={augustus_config}",
        f"--extrinsicCfgFile={extrinsic_config}",
        f"--hintsfile={hints}",
        "--alternatives-from-evidence=true",
        "--UTR=off",
        "--exonnames=on",
        "--protein=on",
        "--codingseq=on",
        "--gff3=on",
        "--softmasking=1",
        "--allow_hinted_splicesites=gcag,atac",
        f"--outfile={output}",
        str(genome),
    )
