# Confirmed-only TE policy (v2)

The default `policies/production_v2.json` now uses
`geneinit_te_confirmed_only_v2`. The previous policy is preserved as
`policies/production_v2_overlap_legacy.json` for comparisons. Migration status
remains pending biological evaluation.

The ordinary core, balanced-base and post-QC evidence routes still apply.
Classified TE overlap is reported but has no veto, rescue prerequisite or rank
weight. Unknown function, host/TE conflict and unassessed host evidence retain
GENERAL eligibility. GENERAL is a candidate pool, not proof of a non-TE gene.
Only confirmed TE is isolated before balanced-base selection, fusion detection
and core ranking. Core/base/bridge qualifications remain distinct.

## Evidence input

Set optional `te_evidence: /absolute/path/assembly.te_evidence.json` in the run
config. This release consumes completed evidence; it does **not** launch eggNOG,
DIAMOND, domain scans or intact-TE discovery, and ships **no biologically
calibrated TE-family signatures**. These analyses and calibration must be
performed upstream. Without a bundle, all candidates are NOT_ASSESSED and use
their existing non-TE selection requirements in GENERAL. An empty TE output in
that mode means no model was confirmed, not that no TE genes exist.

The JSON schema identifier is `geneinit_te_evidence_v2`. Required top-level
fields are `schema`, `assembly_id`, `signature_version`, `exception_version`.
Optional maps: `proteins`, `signatures`, `candidates`. The preflight locks the
bundle checksum; Snakemake tracks it through selection, post-QC and release.
All supplied search receipts, including otherwise unused proteins, are checked.
A configured missing file, FAILED search or mismatched identity is an error.

`proteins` is indexed by SHA256 of upper-case translated protein after removing
one terminal `*`. Identical proteins reuse the same results and classification
cache; loci are never merged by protein identity. Each entry may contain:

- `direct` and `host`: `status` is NOT_RUN, COMPLETE_HIT or COMPLETE_NO_HIT.
  Completed searches require `query_sha256`, `tool_version`, `database_sha256`,
  nonempty `options`, and `receipt` identifying the completed upstream run.
  `hits` is a list, nonempty only for COMPLETE_HIT. Receipt/DB metadata must
  describe the actual upstream result; this adapter does not reproduce searches.
- Direct hits: `kind: DIRECT_DOMAIN`, `family`, `domain`, `score`, `coverage`
  (fraction 0–1), `start_aa`. Direct TE_PROTEIN hits with `evalue`, `identity`
  (percent), `query_coverage` (fraction) and `alignment_aa` may add a suspect flag.
  TRANSFERRED_PFAM and eggNOG descriptions never satisfy direct domain criteria.
- Host hits: `kind: BUSCO`, `status: Complete|Duplicated`; or
  `kind: CELLULAR_SWISSPROT`, `accession`, `identity`, `query_coverage`,
  `subject_coverage`, `alignment_aa`, `evalue`; or a curated
  `CURATED_DOMESTICATED_GENE` match. These normalized kinds must be assigned from
  the same query and actual curated cellular/host evidence, not every Swiss-Prot
  match. Host support prevents TE separation.
- `eggnog_te_terms`: optional annotation terms; suspect flag only.

`signatures` maps TE family to an explicitly reviewed signature containing
`validated: true`, `validation_reference`, `domains` (list of `id`, `min_score`,
`min_coverage`), `ordered` (default true), and `intact_features` (required
family-specific feature names). Every domain must pass its thresholds in the
same protein and the required ordering. Unvalidated/missing families remain
GENERAL. No shared short domain or repeat-overlap threshold substitutes for a
family-specific signature.

`candidates` maps candidate ID to exact `seqid`, `strand`, CDS `start`, CDS
`end`, `cds_signature`, `protein_sha256`, and a `structures` list. Each intact
structure has `family`, `status: INTACT`, `seqid`, `strand`, `start`, `end`,
`source_sha256`, `structure_id` and `features`. This list is normalized from
actual intact-element output and must identify its original source, boundaries,
direction and verified family-specific features. Ordinary repeat-mask overlaps
are not intact-element evidence. All CDS bases must lie within the matching
intact element, whose family agrees with the passing domain architecture.

Optional review: `result: HOST|TE|UNRESOLVED`, `reason`, and
`evidence_snapshot_sha256` computed by
`geneinit.evidence.te_classification.review_snapshot(bundle, context)`.
The snapshot includes query evidence, signatures, versions, coordinates and
intact structures. Stale reviews are flagged and ignored. HOST/UNRESOLVED keep
GENERAL; TE cannot bypass missing intact/host/direct evidence. Resolve an
incorrect host hit in the versioned upstream evidence before a reviewed TE can
separate; an opinion alone is insufficient.

## Post-QC and outputs

New RNA bridges are classified using their newly translated CDS and exact
coordinates before joining the general graph. Their parents and bridge-only
route are recorded. New proteins without completed evidence are NOT_ASSESSED,
GENERAL; old classifications are never inherited. Deterministic signatures
suppress duplicates; the existing structural convergence limit remains five
iterations. No separate classification loop is needed for the fixed evidence
bundle, and uncertainty flags do not trigger another iteration.

The regular release GFF carries `candidate_pool`, `te_class`, `host_state` and
nonblocking `te_review_flags`. Additional release files:

- `te_annotation.gff3`, `te_cds.fna`, `te_proteins.faa`: only confirmed TE models
  passing their structural/evidence route and graph selection, with ORF checks.
- `te_candidate_audit.tsv`: each raw or accepted generated candidate once,
  identity, parents, origin/allowed route, evidence conditions, class, pool,
  final disposition, blocker and final IDs. Existing selector decisions remain
  in the upstream selection audit. Long-intron bridge proposals that fail the
  bridge route remain in the existing bridge review table.
- `te_summary.json`: raw/generated candidate counts, general and TE gene and
  transcript counts, final dispositions, and general annotation class counts.

Mixed intervals receive nonblocking MIXED_REVIEW cross-references. TE genes
must not simply be added to general genes when comparing biological locus
counts. Selection completion is not a new final annotation or BUSCO evaluation.

Family thresholds, source receipts, intact feature extraction and biological
precision require dataset-specific validation before a scientific claim of
confirmed TE removal. The unit-test TEST_TE signature is synthetic and must
never be used as a biological signature.
