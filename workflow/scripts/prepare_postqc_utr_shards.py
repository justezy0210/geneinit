from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402
from geneinit.utr.core import prepare_postqc_utr_shards  # noqa: E402


shards = [str(value) for value in snakemake.params.shards]
seqids = [str(value) for value in snakemake.params.seqids]
if len(shards) != len(seqids) or len(set(shards)) != len(shards) or len(set(seqids)) != len(seqids):
    raise ValueError("UTR seqid/shard mapping is not one-to-one")
seqid_to_shard = dict(zip(seqids, shards, strict=True))
coding_paths = {
    shard: Path(path) for shard, path in zip(shards, snakemake.output.coding, strict=True)
}
rna_paths = {
    shard: Path(path) for shard, path in zip(shards, snakemake.output.rna, strict=True)
}

rows = prepare_postqc_utr_shards(
    Path(snakemake.input.gff3),
    Path(snakemake.input.transcript_decisions),
    Path(snakemake.input.candidates),
    Path(snakemake.input.selection_audit),
    Path(snakemake.input.merged_rna),
    seqid_to_shard,
    coding_paths,
    rna_paths,
    Path(snakemake.output.mapping),
    Path(snakemake.output.tasks),
)
write_completion(
    Path(snakemake.output.complete),
    (
        ("schema_version", 1),
        ("status", "PASS_MIGRATION_PENDING"),
        ("postqc_transcript_count", len(rows)),
        ("reference_sequence_count", len(seqids)),
        ("postqc_gff3_sha256", sha256_file(Path(snakemake.input.gff3))),
        ("transcript_decisions_sha256", sha256_file(Path(snakemake.input.transcript_decisions))),
        ("postqc_complete_sha256", sha256_file(Path(snakemake.input.postqc_complete))),
        ("preparation_mapping_sha256", sha256_file(Path(snakemake.output.mapping))),
        ("task_manifest_sha256", sha256_file(Path(snakemake.output.tasks))),
    ),
)
