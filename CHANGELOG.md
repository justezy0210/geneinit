# Changelog

## 0.0.1 — 2026-10-21

Initial versioned development release.

- Include the opt-in v10 selector, 50-aa single-exon minimum, corroborated RNA
  and protein ranking, and shared selection rules for original and repaired models.
- Separate input and policy validation, reuse evidence context, and track code
  dependencies so reruns start at the affected workflow stage.
- Exclude private regression data, containers, keys, generated sites, and local
  outputs from Git; skip only tests requiring unavailable private fixtures.
- Add an introductory README, execution reference, and explicit validation scope.
- Report the GeneInIt version with `geneinit-run --version` without launching
  Snakemake.

The software version does not promote a biological policy to production status.
Broader species accuracy and a fresh run through every external tool remain
outside the completed fixed-input validation scope.

### Included development changes

- Fix TE release grouping to use transcript IDs consistently and map components
  back to candidate IDs for output; cover two selected TE isoforms through
  complete structural release.
- Match an ordered TE signature against a valid increasing sequence of domain
  hits, including repeated domains and additional earlier hits.
- Reject FASTQ file aliases and byte-identical copies reused across samples.
- Refresh preflight audit records on configuration changes and include effective
  prediction/training and RNA thresholds in semantic hash schema version 2.
- Require checkpoint schema version 2 with input content hashes for AUGUSTUS
  resume; reject version 1 checkpoints rather than trusting file timestamps.
- Run all unit cases with pytest in CI, include pysam in the bootstrap lock,
  and fix the masking-rescue rule's repository-path lint.

- Add explicit RNA BAM + CSI input mode with library metadata, reference sequence
  validation, content locking, and read-only symlink staging. Rebuild StringTie
  and junction evidence through the existing downstream workflow; retain the
  default FASTQ route and current TE/RNA rescue policies.

- Create the standalone GeneInIt repository scaffold.
- Add versioned input schemas and a fail-closed preflight implementation.
- Add the initial Snakemake production profile and Python unit tests.
- Add PE/SE protocol contracts and Snakemake RNA evidence rules for fastp,
  HISAT2, samtools, and StringTie.
- Split the RNA toolchain into three resolved Conda environments because
  StringTie 3.0.3 requires `htslib <1.24` while samtools 1.24 requires
  `htslib 1.24`.
- Add a checksum- and license-gated predictor deployment contract.
- Add taxon-neutral GeneMark-ETP softmask execution with a protected staged
  StringTie GTF, direct GeneMark non-HC coordinate restoration, deterministic
  full ETP hint reconstruction, and direct AUGUSTUS training/prediction.
- Force GeneMark-ETP helper scripts to resolve the image's system Perl before
  Conda Perl, and fail container preflight unless the required CPAN modules and
  every ETP Perl entry point load successfully.
- Add reusable, path-independent predictor and all-in-one container recipes;
  keep the GeneMark key out of images and host homes during build validation,
  and lock the previous Perl and split-environment fixes with regression tests.
- Remove the fixed `--fungus` option and all BRAKER program/helper dependencies;
  lock only GeneMark-ETP and AUGUSTUS in the predictor deployment contract.
- Add a resolved predictor-host environment and a public predictor-DAG dry-run.
- Add policy v2 adaptive RNA-bridge intron thresholds, conservative fallback,
  and explicit review retention for non-promoted long-intron candidates.
- Preserve policy v1 for fixed-5-kb historical structural reproduction.
- Make `postqc.maximum_iterations` executable: structural corrections now run
  until a recorded no-change iteration, fail closed at the configured cap, and
  emit an iteration audit table.
- Reorder structural processing to coding-only grouping, iterative post-QC,
  one final UTR reconstruction for every active transcript (including promoted
  bridges), and release; preserve multi-isoform genes without duplicate GFF3
  gene rows.
- Replace whole-table candidate evidence and selection processing with exact,
  threshold-free coordinate-overlap block streaming; lockstep-stream post-QC
  candidate/evidence/audit inputs and retain only selected or bridge-eligible
  rows.
- Preserve active-fusion boundaries through post-QC, including exclusion from
  strong-RNA bridge anchors and rejection of synthetic bridges that span an
  active-fusion model; audit promotion-created final edges, make canonical
  ranking explicitly coding-evidence-only, and fail closed on UTR/release
  gene-hierarchy drift.
- Validate candidate IDs with an exact disk-backed uniqueness set across the
  complete streamed table, so non-adjacent duplicates fail closed without
  loading all identifiers into memory; remove suffix-list copies from pairwise
  post-QC and selection scans.
- Require synthetic strong-RNA bridges to link at least two actual released
  component roots under the same `assess_edge().linked()` rule used for
  predictor-supported bridges. Derive promotion hit counts from those links
  instead of assuming two exact-boundary anchors, preventing microexon
  overlaps below the coding-link threshold from being falsely promoted.
- Generalize final UTR BAM metadata to PE/SE and stranded/unstranded protocols;
  record InGenAnnot's orientation-neutral fallback for forward-stranded data.
- Make the default Snakemake target the final structural release; retain
  `preflight` as an explicit validation-only target.

## 2026-09-11 — restore strong-RNA classified-TE rescue

- Restore legacy HIGH/SEVERE TE rescue from replicated strong RNA without
  requiring protein tier >=2; report `STRONG_RECURRENT_RNA`.
- Preserve classified-TE CDS scope, candidate RNA calibration and thresholds,
  ORF/length QC, protein uniqueness, and graph/fusion/isoform checks.
- Retain the full unique protein and HIGH recurrent-RNA/unique-protein routes.

## 2026-09-11 — strong RNA and linked multi-hit protein TE rescue

- Rescue HIGH and SEVERE classified-TE models with replicated strong RNA,
  protein tier >=2, and a linked retained target despite ambiguous multi-hit
  protein status; report `STRONG_RNA_MODERATE_LINKED_MULTI_HIT_PROTEIN`.
- Preserve hard QC, RNA thresholds, the recorded protein multiplicity status,
  protein uniqueness routes, and graph/fusion/isoform checks. RNA-only models,
  unresolved targets, weak protein chains, and recurrent-only RNA do not gain
  the new rescue route.

## 2026-09-11 — classified TE and candidate RNA gates

- Use classified TE CDS unions for repeat gates; preserve total-repeat audit metrics.
- Evaluate replicated RNA support per candidate, independent of unrelated singleton conditions.
- Require independent protein support for RNA-based high/severe TE rescue and report reasons.
- Reject legacy component evidence; stamp repeat scope in workflow metadata.
- Include container patch sources in the all-in-one build context.

## 2026-09-13

- Apply confirmed-only TE policy across core/base selection, fusion graphs,
  bridge recovery and release. Preserve uncertainty as nonblocking annotation.
- Add versioned TE evidence import, candidate-specific intact context, separate
  TE GFF/FASTA and raw/generated candidate disposition audit.
- Preserve prior overlap policy for comparisons; no upstream searches or
  calibrated family signatures are implied by the code update.
