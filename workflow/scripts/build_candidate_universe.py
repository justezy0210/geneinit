from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.models.candidate import (  # noqa: E402
    build_candidate_universe,
    fasta_records,
    parse_predictor_gtf,
    read_fasta,
    write_candidate_gff3,
    write_fasta,
    write_source_models,
    write_summary,
    write_universe,
)
from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402
from geneinit.predictor.masking import optional_augustus_models  # noqa: E402


genome = read_fasta(Path(snakemake.input.genome))
augustus_models = parse_predictor_gtf(
    Path(snakemake.input.augustus), "augustus", genome
)
genemark_models = parse_predictor_gtf(
    Path(snakemake.input.genemark), "genemark_etp", genome
)
rescue_models = []
for path in snakemake.input.masking_rescue:
    # Same predictor, distinct source IDs: never inflate predictor concordance.
    rescue_models.extend(optional_augustus_models(Path(path), genome))
models = sorted(
    augustus_models + genemark_models + rescue_models,
    key=lambda model: (
        model.seqid,
        model.start,
        model.end,
        model.strand,
        model.structure_id,
        model.predictor,
        model.source_transcript_id,
    ),
)
universe = build_candidate_universe(models)

write_source_models(Path(snakemake.output.sources), models)
write_universe(Path(snakemake.output.universe), universe)
write_candidate_gff3(Path(snakemake.output.gff3), universe)
write_fasta(Path(snakemake.output.cds), fasta_records(universe, "cds"))
write_fasta(Path(snakemake.output.protein), fasta_records(universe, "protein"))
write_summary(Path(snakemake.output.summary), models, universe)

write_completion(
    Path(snakemake.output.complete),
    (
        ("assembly_id", snakemake.params.assembly_id),
        ("source_model_count", len(models)),
        ("candidate_structure_count", len(universe)),
        ("exact_structure_duplicates_collapsed", len(models) - len(universe)),
        ("source_models_sha256", sha256_file(Path(snakemake.output.sources))),
        ("candidate_universe_sha256", sha256_file(Path(snakemake.output.universe))),
        ("candidate_gff3_sha256", sha256_file(Path(snakemake.output.gff3))),
        ("candidate_cds_sha256", sha256_file(Path(snakemake.output.cds))),
        ("candidate_protein_sha256", sha256_file(Path(snakemake.output.protein))),
    ),
)
