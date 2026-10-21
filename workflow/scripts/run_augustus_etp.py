"""Train and run AUGUSTUS directly from GeneMark-ETP training genes and hints."""

from collections import defaultdict
from pathlib import Path
import json
import shutil
import sys
import traceback


sys.path.insert(0, str(Path(snakemake.input.extrinsic_config).resolve().parents[2]))

from geneinit.contracts.deployment import (  # noqa: E402
    apptainer_command,
    load_predictor_deployment,
)
from geneinit.models.candidate import (  # noqa: E402
    CandidateModel,
    parse_predictor_gtf,
    read_fasta,
    transcript_order,
)
from geneinit.predictor.common import (  # noqa: E402
    PredictorError,
    fasta_lengths,
    require_nonempty,
    run_logged,
    sha256_file,
    validate_gff_coordinates,
    write_completion,
)
from geneinit.predictor.commands import (  # noqa: E402
    augustus_etraining_command,
    augustus_new_species_command,
    augustus_optimize_command,
    augustus_prediction_command,
    augustus_training_genbank_command,
)
from geneinit.predictor.parallel import predict_parallel  # noqa: E402
from geneinit.predictor.resume import load_prediction_resume, model_hashes  # noqa: E402
from geneinit.predictor.training import (  # noqa: E402
    evaluation_metrics,
    read_training_records,
    require_clean_training,
    require_matching_selection,
    split_training,
    training_errors,
    write_genbank,
    write_proteins,
)


def safe_attribute(value: str, label: str) -> str:
    if not value or any(character in value for character in "\t\n;,=\""):
        raise PredictorError(f"unsafe AUGUSTUS {label}: {value!r}")
    return value


def write_fasta(path: Path, models: list[CandidateModel], kind: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for model in models:
            identifier = safe_attribute(model.source_transcript_id, "transcript ID")
            sequence = (
                model.orf.protein_sequence
                if kind == "protein"
                else model.orf.coding_sequence
            )
            handle.write(f">{identifier}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")
    temporary.replace(path)


def write_normalized_predictions(
    gff3: Path, gtf: Path, models: list[CandidateModel]
) -> tuple[int, int]:
    by_gene: dict[str, list[CandidateModel]] = defaultdict(list)
    for model in models:
        by_gene[safe_attribute(model.source_gene_id, "gene ID")].append(model)

    gff_tmp = gff3.with_name(f".{gff3.name}.tmp")
    gtf_tmp = gtf.with_name(f".{gtf.name}.tmp")
    ordered_genes = sorted(
        by_gene.items(),
        key=lambda item: (
            item[1][0].seqid,
            min(model.start for model in item[1]),
            max(model.end for model in item[1]),
            item[0],
        ),
    )
    with gff_tmp.open("w", encoding="utf-8") as gff_handle, gtf_tmp.open(
        "w", encoding="utf-8"
    ) as gtf_handle:
        gff_handle.write("##gff-version 3\n")
        for gene_id, members in ordered_genes:
            seqids = {model.seqid for model in members}
            strands = {model.strand for model in members}
            if len(seqids) != 1 or len(strands) != 1:
                raise PredictorError(f"AUGUSTUS gene ID reused across loci: {gene_id}")
            seqid = members[0].seqid
            strand = members[0].strand
            gene_start = min(model.start for model in members)
            gene_end = max(model.end for model in members)
            gff_handle.write(
                f"{seqid}\tAUGUSTUS\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t."
                f"\tID={gene_id}\n"
            )
            for model in sorted(members, key=lambda item: (item.start, item.end, item.source_transcript_id)):
                transcript_id = safe_attribute(
                    model.source_transcript_id, "transcript ID"
                )
                gff_handle.write(
                    f"{seqid}\tAUGUSTUS\tmRNA\t{model.start}\t{model.end}\t.\t{strand}\t."
                    f"\tID={transcript_id};Parent={gene_id}\n"
                )
                gtf_attributes = (
                    f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                )
                gtf_handle.write(
                    f"{seqid}\tAUGUSTUS\ttranscript\t{model.start}\t{model.end}\t.\t"
                    f"{strand}\t.\t{gtf_attributes}\n"
                )
                exon_parts = model.exon_parts or tuple(
                    (part.start, part.end) for part in model.cds_parts
                )
                for index, (start, end) in enumerate(
                    sorted(exon_parts, reverse=strand == "-"), 1
                ):
                    exon_id = f"{transcript_id}.exon{index}"
                    gff_handle.write(
                        f"{seqid}\tAUGUSTUS\texon\t{start}\t{end}\t.\t{strand}\t."
                        f"\tID={exon_id};Parent={transcript_id}\n"
                    )
                    gtf_handle.write(
                        f"{seqid}\tAUGUSTUS\texon\t{start}\t{end}\t.\t{strand}\t."
                        f"\t{gtf_attributes}\n"
                    )
                for index, part in enumerate(
                    transcript_order(model.cds_parts, strand), 1
                ):
                    cds_id = f"{transcript_id}.cds{index}"
                    gff_handle.write(
                        f"{seqid}\tAUGUSTUS\tCDS\t{part.start}\t{part.end}\t.\t{strand}\t"
                        f"{part.phase}\tID={cds_id};Parent={transcript_id}\n"
                    )
                    gtf_handle.write(
                        f"{seqid}\tAUGUSTUS\tCDS\t{part.start}\t{part.end}\t.\t{strand}\t"
                        f"{part.phase}\t{gtf_attributes}\n"
                    )
    gff_tmp.replace(gff3)
    gtf_tmp.replace(gtf)
    return len(ordered_genes), len(models)


def genbank_locus_count(path: Path) -> int:
    return sum(1 for line in path.open(encoding="utf-8") if line.startswith("LOCUS"))


log = Path(snakemake.log[0])
try:
    deployment = load_predictor_deployment(dict(snakemake.config))
    genome = Path(snakemake.input.genome).resolve()
    training_gtf = Path(snakemake.input.training).resolve()
    hints = Path(snakemake.input.hints).resolve()
    extrinsic_config = Path(snakemake.input.extrinsic_config).resolve()
    config_dir = Path(snakemake.input.augustus_config).resolve()
    work = Path(snakemake.output.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    stage_logs = work / "logs"
    stage_logs.mkdir(parents=True, exist_ok=True)

    species = str(snakemake.params.species)
    species_dir = config_dir / f"species/{species}"
    resume = None
    if snakemake.input.prediction_resume:
        resume = load_prediction_resume(
            Path(snakemake.input.prediction_resume[0]), genome=genome, hints=hints,
            training_gtf=training_gtf, extrinsic=extrinsic_config, species=species,
            species_dir=species_dir, work=work, settings={
                "rounds": int(snakemake.params.rounds), "kfold": int(snakemake.params.kfold),
                "training_flank": int(snakemake.params.training_flank),
                "training_selection": dict(snakemake.params.training_selection),
                "genetic_code": int(snakemake.config.get("genetic_code", 1)),
            },
        )
        selection = json.loads((work / "training/selection.json").read_text())
        evaluation = json.loads((work / "training/heldout_evaluation.json").read_text())
        initial_evaluation, final_evaluation = evaluation["initial"], evaluation["optimized"]
        training_loci = selection["splits"]["onlytrain"]["loci"] + selection["splits"]["optimize"]["loci"]
        kfold, optimize_threads = resume["kfold"], resume["optimization_cpus"]
    else:
        if species_dir.exists():
            shutil.rmtree(species_dir)

        commands = [
            (
                "new_species",
                augustus_new_species_command(
                    species=species, augustus_config=config_dir
                ),
                (config_dir,),
            )
        ]
        for name, inner_command, binds in commands:
            command = apptainer_command(deployment, inner_command, binds, pwd=work)
            if any("fungus" in value.lower() or "braker" in value.lower() for value in command):
                raise PredictorError("taxon-specific or BRAKER option leaked into direct AUGUSTUS")
            run_logged(command, stage_logs / f"{name}.log")

        if int(snakemake.config.get("genetic_code", 1)) != 1:
            raise PredictorError("AUGUSTUS training ORF selection currently requires genetic_code=1")
        raw_training_gb = work / "training.raw.gb"
        if snakemake.input.training_snapshot:
            snapshot = Path(snakemake.input.training_snapshot[0]).resolve()
            require_nonempty(snapshot, "reviewed training GenBank")
            shutil.copyfile(snapshot, raw_training_gb)
            (stage_logs / "training_snapshot.json").write_text(json.dumps({
                "genbank": str(snapshot),
                "selection_reference": str(Path(snakemake.input.selection_reference[0]).resolve()),
                "conversion_reused": True,
            }, indent=2) + "\n", encoding="utf-8")
        else:
            training_command = augustus_training_genbank_command(
                training_gtf=training_gtf,
                genome=genome,
                flank=int(snakemake.params.training_flank),
                output=raw_training_gb,
            )
            run_logged(
                apptainer_command(
                    deployment,
                    training_command,
                    (work, genome, training_gtf),
                    pwd=work,
                ),
                stage_logs / "gff2gbSmallDNA.log",
            )
        require_nonempty(raw_training_gb, "AUGUSTUS training GenBank")
        records = read_training_records(raw_training_gb)
        training_dir = work / "training"
        training_dir.mkdir(exist_ok=True)
        quality_gb = training_dir / "quality.gb"
        eligible = [record for record in records if not record.reason]
        if len(eligible) < 8:
            raise PredictorError("fewer than 8 complete ORFs with canonical splice sites")
        write_genbank(quality_gb, eligible)
        screen_log = stage_logs / "etraining.screen.log"
        screen = augustus_etraining_command(
            species=species, training_genbank=quality_gb, augustus_config=config_dir
        )
        run_logged(apptainer_command(deployment, screen, (work, config_dir), pwd=work), screen_log)
        excluded = training_errors(screen_log)
        if set(excluded) - {record.locus for record in eligible}:
            raise PredictorError(f"unrecognized locus in etraining screening log: {screen_log}")
        if any("end in stop codon" in message for message in excluded.values()):
            raise PredictorError(f"stop-codon errors remain after training normalization; inspect {screen_log}")
        eligible = [record for record in eligible if record.locus not in excluded]
        if len(eligible) < 8:
            raise PredictorError("fewer than 8 loci survive etraining structure screening")
        proteins = training_dir / "eligible.faa"
        write_proteins(proteins, eligible)
        protein_db = training_dir / "eligible"
        protein_hits = training_dir / "similarity.tsv"
        training_options = dict(snakemake.params.training_selection)
        allowed_options = {"max_loci", "optimize_loci", "holdout_loci", "seed", "identity", "coverage"}
        if set(training_options) - allowed_options:
            raise PredictorError(f"unknown AUGUSTUS training options: {set(training_options) - allowed_options}")
        diamond = "/opt/ETP/tools/diamond"
        for name, command in (
            ("training_makedb", (diamond, "makedb", "--in", str(proteins), "--db", str(protein_db),
                                "--threads", str(snakemake.threads))),
            ("training_similarity", (
                diamond, "blastp", "--query", str(proteins), "--db", str(protein_db),
                "--out", str(protein_hits), "--threads", str(snakemake.threads),
                "--sensitive", "--max-target-seqs", "0", "--evalue", "1e-5",
                "--id", str(training_options.get("identity", 80)),
                "--query-cover", str(training_options.get("coverage", 80)),
                "--subject-cover", str(training_options.get("coverage", 80)),
                "--outfmt", "6", "qseqid", "sseqid", "pident", "qstart", "qend", "qlen",
                "sstart", "send", "slen",
            )),
        ):
            run_logged(apptainer_command(deployment, command, (work,), pwd=work), stage_logs / f"{name}.log")
        selection = split_training(records, protein_hits, training_dir, excluded, **training_options)
        if snakemake.input.selection_reference:
            require_matching_selection(
                training_dir / "selection.tsv", Path(snakemake.input.selection_reference[0])
            )
        training_gb = training_dir / "train.gb"
        optimize_gb = training_dir / "optimize.gb"
        onlytrain_gb = training_dir / "onlytrain.gb"
        holdout_gb = training_dir / "holdout.gb"
        training_loci = genbank_locus_count(training_gb)
        if training_loci < 8:
            raise PredictorError(
                f"AUGUSTUS optimization requires at least 8 non-overlapping training loci; observed={training_loci}"
            )
        kfold = min(int(snakemake.params.kfold), selection["splits"]["optimize"]["loci"])
        optimize_threads = min(max(1, snakemake.threads), kfold)

        # Discard screening parameters: no held-out locus may influence the model.
        shutil.rmtree(species_dir)
        run_logged(apptainer_command(deployment, augustus_new_species_command(
            species=species, augustus_config=config_dir
        ), (config_dir,), pwd=work), stage_logs / "new_species.selected.log")
        etraining = augustus_etraining_command(
            species=species,
            training_genbank=training_gb,
            augustus_config=config_dir,
        )
        run_logged(
            apptainer_command(deployment, etraining, (work, config_dir), pwd=work),
            stage_logs / "etraining.initial.log",
        )
        require_clean_training(stage_logs / "etraining.initial.log")

        def evaluate_holdout(name):
            evaluation_log = stage_logs / f"heldout.{name}.log"
            command = (
                "/opt/Augustus/bin/augustus", f"--species={species}",
                f"--AUGUSTUS_CONFIG_PATH={config_dir}", "--UTR=off", "--softmasking=0",
                "--stopCodonExcludedFromCDS=false", str(holdout_gb),
            )
            run_logged(apptainer_command(deployment, command, (work, config_dir), pwd=work), evaluation_log)
            return evaluation_metrics(evaluation_log)

        initial_evaluation = evaluate_holdout("initial")
        (training_dir / "heldout_initial.json").write_text(
            json.dumps(initial_evaluation, indent=2) + "\n", encoding="utf-8"
        )
        # Preserve the original baseline evaluation. Previously optimized metadata
        # become the starting point for the new optimization only after that test.
        if snakemake.input.optimization_start:
            starting_parameters = Path(snakemake.input.optimization_start[0]).resolve()
            require_nonempty(starting_parameters, "starting optimization parameters")
            shutil.copyfile(starting_parameters, species_dir / f"{species}_parameters.cfg")
            (stage_logs / "optimization_start.json").write_text(json.dumps({
                "path": str(starting_parameters),
                "sha256": sha256_file(starting_parameters),
                "kfold": kfold,
                "cpus": optimize_threads,
            }, indent=2) + "\n", encoding="utf-8")
        optimize = augustus_optimize_command(
            species=species,
            training_genbank=optimize_gb,
            augustus_config=config_dir,
            threads=optimize_threads,
            rounds=int(snakemake.params.rounds),
            kfold=kfold,
            onlytrain=onlytrain_gb,
        )
        run_logged(
            apptainer_command(deployment, optimize, (work, config_dir), pwd=work),
            stage_logs / "optimize_augustus.log",
        )
        run_logged(
            apptainer_command(deployment, etraining, (work, config_dir), pwd=work),
            stage_logs / "etraining.final.log",
        )
        require_clean_training(stage_logs / "etraining.final.log")
        final_evaluation = evaluate_holdout("optimized")
        (training_dir / "heldout_evaluation.json").write_text(json.dumps({
            "reference": "independent GeneMark-ETP HC loci; not experimental ground truth",
            "holdout_loci": selection["splits"]["holdout"]["loci"],
            "initial": initial_evaluation, "optimized": final_evaluation,
            "f1_change": {level: final_evaluation[level]["f1"] - initial_evaluation[level]["f1"]
                          for level in initial_evaluation},
        }, indent=2) + "\n", encoding="utf-8")

    raw_gff3 = work / "augustus.direct.raw.gff3"
    def prediction_for(shard):
        prediction = augustus_prediction_command(
            species=species, genome=shard.fasta, hints=shard.hints,
            extrinsic_config=extrinsic_config, augustus_config=config_dir,
            output=shard.output,
        )
        if any("fungus" in value.lower() or "braker" in value.lower() for value in prediction):
            raise PredictorError("taxon-specific or BRAKER option leaked into direct AUGUSTUS")
        return apptainer_command(
            deployment,
            prediction,
            (work, extrinsic_config, config_dir),
            pwd=work,
        )

    prediction_summary = predict_parallel(
        genome=genome, hints=hints, work=work, output=raw_gff3,
        threads=snakemake.threads, command_for=prediction_for,
        wrap_command=lambda command, directory: apptainer_command(
            deployment, command, (directory,), pwd=work,
        ),
        log=stage_logs / "augustus_prediction.log",
        completed_sequences=resume["completed_sequences"] if resume else None,
    )
    require_nonempty(raw_gff3, "direct AUGUSTUS GFF3")
    validate_gff_coordinates(raw_gff3, fasta_lengths(genome))

    if resume and model_hashes(species_dir) != resume["model_sha256"]:
        raise PredictorError("trained model changed during resumed prediction")

    models = parse_predictor_gtf(raw_gff3, "augustus", read_fasta(genome))
    gff3 = Path(snakemake.output.gff3).resolve()
    gtf = Path(snakemake.output.gtf).resolve()
    gff3.parent.mkdir(parents=True, exist_ok=True)
    genes, transcripts = write_normalized_predictions(gff3, gtf, models)
    protein = Path(snakemake.output.protein).resolve()
    codingseq = Path(snakemake.output.codingseq).resolve()
    write_fasta(protein, models, "protein")
    write_fasta(codingseq, models, "codingseq")
    lengths = fasta_lengths(genome)
    validate_gff_coordinates(gff3, lengths)
    validate_gff_coordinates(gtf, lengths)

    require_nonempty(
        species_dir / f"{species}_parameters.cfg", "trained AUGUSTUS species"
    )
    trained = Path(snakemake.output.trained).resolve()
    if trained.exists():
        shutil.rmtree(trained)
    shutil.copytree(species_dir, trained)

    write_completion(
        Path(snakemake.output.complete),
        (
            ("schema_version", 2),
            ("status", "PASS"),
            ("assembly_id", snakemake.params.assembly_id),
            ("species", species),
            ("mode", "direct GeneMark-ETP training genes -> AUGUSTUS"),
            ("taxon_mode", "none"),
            ("braker_used", "false"),
            ("container_sha256", deployment.sif_sha256),
            ("rounds", snakemake.params.rounds),
            ("training_flank", snakemake.params.training_flank),
            ("training_loci", training_loci),
            ("training_reused", bool(resume)),
            ("training_container_sha256", resume["training_container_sha256"] if resume else deployment.sif_sha256),
            ("training_snapshot_reused", bool(snakemake.input.training_snapshot)),
            ("training_selection_method", selection["method"]),
            ("optimization_evaluation_loci", selection["splits"]["optimize"]["loci"]),
            ("onlytrain_loci", selection["splits"]["onlytrain"]["loci"]),
            ("holdout_loci", selection["splits"]["holdout"]["loci"]),
            ("initial_holdout_gene_f1", initial_evaluation["gene"]["f1"]),
            ("optimized_holdout_gene_f1", final_evaluation["gene"]["f1"]),
            ("kfold", kfold),
            ("optimization_cpus", optimize_threads),
            ("prediction_cpus", prediction_summary["maximum_active_workers"]),
            ("prediction_sequences", prediction_summary["sequence_count"]),
            ("prediction_reused_sequences", prediction_summary["reused_sequence_count"]),
            ("optimization_start_parameters", str(snakemake.input.optimization_start[0])
             if snakemake.input.optimization_start else ""),
            ("genes", genes),
            ("transcripts", transcripts),
            ("proteins", transcripts),
            ("gff3_sha256", sha256_file(gff3)),
            ("gtf_sha256", sha256_file(gtf)),
            ("protein_sha256", sha256_file(protein)),
            ("codingseq_sha256", sha256_file(codingseq)),
            ("extrinsic_config_sha256", sha256_file(extrinsic_config)),
        ),
    )
except Exception:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("\n" + traceback.format_exc())
    raise
else:
    log.write_text("direct AUGUSTUS training and prediction PASS\n", encoding="utf-8")
