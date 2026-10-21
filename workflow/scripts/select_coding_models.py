from collections import Counter
from contextlib import ExitStack, contextmanager
import csv
import json
from pathlib import Path
import sys


repo_root = Path(snakemake.params.repo_root)
sys.path.insert(0, str(repo_root))

from geneinit.evidence.te_classification import TeClassifier, confirmed_only, protein_sha256
from geneinit.models.candidate import read_fasta

from geneinit.evidence.repeats import REPEAT_SCOPE  # noqa: E402
from geneinit.evidence.junctions import write_metric_table  # noqa: E402
from geneinit.predictor.common import sha256_file, write_completion  # noqa: E402
from geneinit.selection.core import (  # noqa: E402
    AUDIT_COLUMNS,
    EDGE_COLUMNS,
    audit_row,
    build_loci_and_edges,
    edge_row,
    evaluate_and_select,
    iter_component_blocks,
)
from geneinit.selection.balanced import (  # noqa: E402
    BASE_AUDIT_COLUMNS,
    base_audit_row,
    select_balanced_base,
)


@contextmanager
def atomic_dict_writer(path: Path, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            yield writer
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_key_file(path: Path, keys):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial")
    try:
        temporary.write_text(
            "".join(f"{key}\n" for key in keys), encoding="utf-8"
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


policy_path = Path(snakemake.input.policy)
with policy_path.open(encoding="utf-8") as handle:
    policy_document = json.load(handle)
if policy_document.get("selector_schema") != "geneinit_selector_policy_v1":
    raise ValueError("policy does not contain geneinit_selector_policy_v1")
policy = policy_document["selector_core"]

neutral = confirmed_only(policy)
te_inputs = list(getattr(snakemake.input, "te_evidence", []))
classifier = TeClassifier(Path(te_inputs[0]) if te_inputs else None, str(snakemake.params.assembly_id)) if neutral else None
protein_hashes = {key: protein_sha256(seq) for key, seq in read_fasta(Path(snakemake.input.proteins)).items()} if neutral else {}
seen_contexts = set()
te_selected_count = 0
candidate_count = 0
locus_count = 0
edge_count = 0
eligible_count = 0
base_eligible_count = 0
base_supplement_count = 0
active_fusion_count = 0
base_statuses = Counter()
decisions = Counter()
routes = Counter()
base_selected_all: set[str] = set()
selected_order: list[str] = []

with ExitStack() as stack:
    base_audit_writer = stack.enter_context(
        atomic_dict_writer(Path(snakemake.output.base_audit), BASE_AUDIT_COLUMNS)
    )
    loci_writer = stack.enter_context(
        atomic_dict_writer(
            Path(snakemake.output.loci),
            ("candidate_id", "locus_id", "seqid", "strand", "start", "end"),
        )
    )
    edge_writer = stack.enter_context(
        atomic_dict_writer(Path(snakemake.output.edges), EDGE_COLUMNS)
    )
    audit_writer = stack.enter_context(
        atomic_dict_writer(Path(snakemake.output.audit), AUDIT_COLUMNS)
    )

    for records in iter_component_blocks(Path(snakemake.input.components)):
        if classifier is not None:
            for record in records:
                record.row.update(classifier.classify(record.row, protein_hashes[record.key]))
                seen_contexts.add(record.key)
        loci, edges = build_loci_and_edges(
            str(snakemake.params.assembly_id), records
        )
        base_selected, base_calls = select_balanced_base(loci, edges, policy)
        selected = evaluate_and_select(
            loci,
            edges,
            policy,
            sample_count=int(snakemake.params.sample_count),
            replicates_per_condition=int(
                snakemake.params.replicates_per_condition
            ),
            base_selected_keys=base_selected,
        )

        ordered_records = sorted(
            records,
            key=lambda item: (
                item.seqid, item.start, item.end, item.strand, item.key
            ),
        )
        for record in ordered_records:
            base_audit_writer.writerow(base_audit_row(base_calls[record.key]))
            loci_writer.writerow(
                {
                    "candidate_id": record.key,
                    "locus_id": record.locus_id,
                    "seqid": record.seqid,
                    "strand": record.strand,
                    "start": record.start,
                    "end": record.end,
                }
            )
            audit_writer.writerow(audit_row(record, str(policy["policy_id"])))

        for edge in edges:
            edge_writer.writerow(edge_row(edge))

        te_selected_count += sum(record.row.get("te_selected") == "1" for record in records)
        candidate_count += len(records)
        locus_count += len(loci)
        edge_count += len(edges)
        eligible_count += sum(record.eligible for record in records)
        base_eligible_count += sum(call.eligible for call in base_calls.values())
        base_supplement_count += sum(record.base_supplement for record in records)
        active_fusion_count += sum(record.active_fusion for record in records)
        base_statuses.update(call.status for call in base_calls.values())
        decisions.update(record.decision for record in records)
        routes.update(route for record in records for route in record.routes)
        base_selected_all.update(base_selected)
        selected_order.extend(record.key for record in selected)

if classifier is not None:
    unexpected = set(classifier.bundle.get("candidates", {})) - seen_contexts
    # GRNA contexts may describe reproducible post-QC generated models.
    unexpected = {key for key in unexpected if not key.startswith("GRNA-")}
    if unexpected:
        raise ValueError(f"TE evidence contains unknown candidates: {sorted(unexpected)[:10]}")
if len(selected_order) != len(set(selected_order)):
    raise ValueError("streamed selection produced duplicate selected candidates")
write_key_file(Path(snakemake.output.base_selected), sorted(base_selected_all))
write_key_file(Path(snakemake.output.selected), selected_order)

base_summary_rows = [
    ("assembly_id", snakemake.params.assembly_id),
    ("candidate_count", candidate_count),
    ("locus_count", locus_count),
    ("edge_count", edge_count),
    ("eligible_count", base_eligible_count),
    ("selected_count", len(base_selected_all)),
]
base_summary_rows.extend(
    (f"decision::{key}", value) for key, value in sorted(base_statuses.items())
)
write_metric_table(Path(snakemake.output.base_summary), base_summary_rows)

summary_rows = [
    ("assembly_id", snakemake.params.assembly_id),
    ("policy_id", policy["policy_id"]),
    ("policy_status", policy_document["status"]),
    ("candidate_count", candidate_count),
    ("locus_count", locus_count),
    ("edge_count", edge_count),
    ("eligible_count", eligible_count),
    ("selected_count", len(selected_order)),
    ("te_selected_count", te_selected_count),
    ("balanced_base_selected_count", len(base_selected_all)),
    ("base_supplement_count", base_supplement_count),
    ("active_fusion_count", active_fusion_count),
]
summary_rows.extend(
    (f"decision::{key}", value) for key, value in sorted(decisions.items())
)
summary_rows.extend(
    (f"route::{key}", value) for key, value in sorted(routes.items())
)
write_metric_table(Path(snakemake.output.summary), summary_rows)

write_completion(
    Path(snakemake.output.complete),
    (
        ("schema_version", 2),
        ("repeat_scope", REPEAT_SCOPE),
        ("rna_support_basis", "candidate_condition_junctions"),
        ("status", "PASS_MIGRATION_PENDING"),
        ("assembly_id", snakemake.params.assembly_id),
        ("policy_id", policy["policy_id"]),
        ("policy_sha256", sha256_file(policy_path)),
        ("candidate_count", candidate_count),
        ("selected_count", len(selected_order)),
        ("balanced_base_selected_count", len(base_selected_all)),
        ("base_supplement_count", base_supplement_count),
        ("base_audit_sha256", sha256_file(Path(snakemake.output.base_audit))),
        ("base_selected_keys_sha256", sha256_file(Path(snakemake.output.base_selected))),
        ("base_summary_sha256", sha256_file(Path(snakemake.output.base_summary))),
        ("loci_sha256", sha256_file(Path(snakemake.output.loci))),
        ("edges_sha256", sha256_file(Path(snakemake.output.edges))),
        ("audit_sha256", sha256_file(Path(snakemake.output.audit))),
        ("selected_keys_sha256", sha256_file(Path(snakemake.output.selected))),
        ("summary_sha256", sha256_file(Path(snakemake.output.summary))),
    ),
)
