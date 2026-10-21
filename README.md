# GeneInIt

The current Jinju-data recommendation uses
`policies/unified_v10_single_exon50_corroborated_multisample.json` with one
selector for original and repaired candidates. The maintained workflow separates
input validation from downstream policy validation, shares evidence scoring, and
tracks imported code so reruns start at the affected stage. See
[incremental execution and validation](docs/development/incremental_execution.md).
The original production-v2 default remains available; select the v10 policy
explicitly. Its frozen experimental marker is retained, with completed validation
recorded outside the policy file.

GeneInIt is a reproducible Snakemake workflow for evidence-integrated eukaryotic
structural gene annotation. The intended production input is one nuclear,
soft-masked assembly, its coordinate-matched repeat GFF3, RNA-seq evidence
(raw reads or precomputed BAM + CSI), and an external protein evidence FASTA.

> Development status: preflight, RNA evidence, licensed predictor orchestration,
> candidate/evidence integration, coding-only provisional grouping, iterative
> transcript-first whole-genome structural post-QC, final InGenAnnot UTR
> reconstruction, and deterministic structural release staging are implemented
> and pass unit/lint/dry-run validation. GeneMark-ETP and AUGUSTUS are invoked
> directly without BRAKER or a fixed taxon mode. Candidate-heavy stages use the
> same threshold-free coordinate-block streaming path for every assembly. A real
> licensed predictor image is still required to validate a complete clean run.

## Core boundaries

- One assembly per run
- Input FASTA coordinates are the canonical output coordinates
- Structural annotation only
- No repeat discovery, functional annotation, BUSCO, OrthoFinder, synteny, or
  manual decision ledger in the core DAG
- No organism-name branch, fungus/plant mode, candidate-count cutoff, or
  agent/manual intervention in structural decisions
- No GeneMark key, private container, or private biological data is bundled

## Classified TE and candidate RNA evidence

Repeat gates use only known TE orders in the repeat GFF3, measured as the union
of overlapping CDS bases. `unknown`/unclassified repeats, RNA genes, simple/low
complexity repeats, and tentative orders such as `SINE?` do not contribute.
A known order is sufficient: `LTR/unknown` and `DNA/unknown` still count as TE.
Explicit `Classification`/`class` attributes take precedence over a library
`Name#class/family` suffix or recognized TE feature types. Arbitrary feature
names are never treated as TE classification. Supported orders are LTR, LINE,
SINE, DNA, RC, MITE, TIR, Helitron, Penelope, DIRS and Crypton. See the
[EDTA output description](https://github.com/oushujun/EDTA#outputs) for its GFF3
classification convention.

`repeat_scope=classified_te_v1` identifies the new component contract;
`repeat_*` metrics now describe classified TEs. `all_repeat_*` metrics retain
the total-repeat union and class breakdown for audit only. All selection,
balanced-base and post-QC repeat gates consume the classified fraction. The
input soft mask is not reconstructed by this change.

The default `production_v2.json` policy admits multi-exon candidates with a
complete ORF, passing structural QC, and recurrent RNA support for at least
50% of their intron connections. A recurrent connection requires at least two
HQ reads in total across at least two samples. The existing RNA tier threshold
also applies. `max_nonrecurrent_introns: null` removes the absolute limit on
missing connections; it does not change the 50% requirement. Single-exon
evidence routes remain unchanged. Graph conflicts, alternative-isoform evidence
and downstream structural QC still determine the final selection.

Strong RNA retains the per-connection requirement of at least three HQ reads
and at least one condition with two supporting replicates, with the same
policy fraction and the existing strong RNA tier requirement. Adding a
singleton condition cannot disable evidence in other replicated conditions.
The current junction contract supports `min_replicates_strong=2`; other values
fail explicitly. Legacy policies keep their 80% fraction and absolute cap.

Selection audits report `rna_recurrent_intron_count`,
`rna_recurrent_intron_fraction`, and `rna_nonrecurrent_introns_json` (1-based,
inclusive genomic coordinates). `RECURRENT_RNA_50` identifies the new route;
legacy 80% policies retain `RECURRENT_RNA_80`. Passing this gate retains the
predicted CDS structure and does not trim unsupported exons or establish that
every splice connection is correct.

The default confirmed-only policy reports classified TE overlap without using
it to reject a general candidate. Confirmed TE models remain in the separate
TE pool even when their RNA support passes 50%. Uncertain or conflicting TE
evidence remains in the general pool. Protein uniqueness is a mapping property,
not proof that a protein is non-TE. Audit rows retain both repeat fractions,
strong-intron counts, and the exact selection reason.

Old component tables must be regenerated with `score_component_evidence`
before selection or post-QC. The new Snakemake repeat-scope parameter invalidates
that stage when normal parameter rerun triggers are enabled; loaders also reject
legacy evidence rather than interpreting total-repeat fractions as TE fractions.
Reusing existing predictions does not require retraining AUGUSTUS, but selection
and its downstream outputs must be regenerated in a separate run/output area.
Single-exon routes and new-model reconstruction are unchanged.

## Precomputed RNA alignments

Set `rna.input_mode: bam`, supply the FASTA used for alignment as
`rna.reference_fasta`, and point `samples` to the BAM manifest described in
[docs/inputs.md](docs/inputs.md#precomputed-rna-bam--csi).
`geneinit-run --threads N --configfile run.yaml rna_import` validates and stages
BAM + CSI only. `rna_evidence` also rebuilds StringTie transcripts;
`junction_evidence` rebuilds junction support. Full runs reuse these alignments
for prediction and UTR reconstruction. No FASTQ files or HISAT2 index are needed
in BAM mode. The default `fastq` mode retains the existing read-processing path.

## Bootstrap

```bash
mamba env create \
  --prefix .conda/bootstrap \
  --file envs/bootstrap.yaml

conda run --prefix .conda/bootstrap \
  snakemake --profile profiles/local \
  --configfile config/config.yaml
```

The resolved RNA environments are intentionally separate:

```text
.conda/fastp
.conda/alignment   # HISAT2 + samtools
.conda/stringtie
.conda/utr         # InGenAnnot 0.0.20; executed with user-site disabled
```

Copy `config/config.example.yaml` to the ignored `config/config.yaml`, replace
all placeholders, and prepare a sample table following
`config/samples.schema.yaml`.

Running the bootstrap command without a target executes the full structural
workflow through `11_release.staging/`. Use the explicit `preflight` target
when only input and deployment validation is wanted.

The default `production_v2.json` policy learns RNA-bridge intron thresholds
from high-confidence introns in the current assembly. The fixed 5 kb policy is
retained as `production_v1.json` for historical structural-threshold
reproduction.

The external licensed predictor deployment must be declared with its local SIF
path, SHA-256, and explicit GeneMark license acknowledgement. GeneInIt checks
the checksum and required in-image programs before expensive work. It never
bundles the image or a GeneMark key.

An all-in-one image may run the same contracts without nested Apptainer. Set
`deployment.mode: in_container`, omit `apptainer_executable`, and run Snakemake
with `profiles/allinone`. The GeneMark key remains a runtime-only home-directory
file and is not embedded in the image.

Reusable, path-independent predictor and all-in-one build recipes are documented
in [`containers/README.md`](containers/README.md). Project directories should
call those scripts with environment variables rather than copy and modify the
Dockerfiles.

GeneInIt runs GeneMark-ETP with `/usr/bin` ahead of Conda environments in
`PATH`, and the reference container recipe pins every ETP Perl entry point to
`/usr/bin/perl`. This is required because GeneMark-ETP 1.02 includes helpers
with environment-based Perl shebangs while its CPAN dependencies belong to the
image's system Perl. Predictor preflight verifies those modules and
syntax-checks every ETP Perl script before expensive prediction begins.

AUGUSTUS training lowercases GenBank DNA (the prediction genome retains its
repeat mask), audits complete ORFs and splice sites, and screens `etraining`
messages even when its exit status is zero. DIAMOND groups proteins at 80%
identity and 80% coverage of both proteins; exact duplicates and overlapping
GenBank genomic intervals are also joined. One representative per connected
group is retained before applying the training size limit. Selection is seeded
and stratified by GC, CDS length, and exon count.

Training excludes models with a start or stop codon split across a splice
junction (first/last CDS exon shorter than three bases), because AUGUSTUS
etraining rejects these structures even if the spliced ORF is complete.
These exclusions do not remove the original GeneMark predictions.

`prediction.augustus_training` defaults to 6,000 total loci: 500 independent
holdout loci, 500 optimization evaluation loci, and 5,000 `--onlytrain` loci.
The holdout is excluded from all parameter training, including screening state,
and is evaluated before and after optimization. Final `etraining` uses the
5,500 non-holdout loci. For smaller datasets the total shrinks automatically;
datasets with no room for `onlytrain` require smaller explicit split sizes.
`work/training/selection.tsv` and `selection.json` record selection/exclusion
reasons and split statistics; `heldout_evaluation.json` records sensitivity,
specificity and F1 against GeneMark-ETP HC models, not experimental ground truth.
This selection currently requires the standard genetic code (`genetic_code: 1`).

To continue from a reviewed selection, set `prediction.augustus_training_snapshot`
with `genbank` pointing to its original converted GenBank and `selection` pointing
to its `selection.tsv`. The pipeline reuses that conversion, repeats quality and
similarity checks, and requires identical assignments before selected-data
training. This avoids changing models when the converter resolves tied,
overlapping annotations differently. Both files must come from the same genome
and GeneMark training input used for the current run. Omit this setting for a new
dataset; the normal conversion path remains available.

`prediction.augustus_prediction_resume` accepts checkpoint schema version 2.
The checkpoint must include SHA-256 hashes for the genome, hints, and training
GTF, in addition to their path, size, and modification time. Version 1 or
hash-free checkpoints are rejected before model files are copied. Recreate a
checkpoint only from verified original inputs, or start prediction again;
changing its version field alone does not verify old predictions.

Run GeneInIt with an explicit total CPU budget:

```bash
# Inside the all-in-one SIF
apptainer exec --cleanenv --bind /data:/data geneinit.sif \
  geneinit-run --threads 34 --configfile /data/config.yaml

# From a source checkout with the workflow dependencies installed
python -m geneinit.cli --threads 34 --configfile config/config.yaml
```

`--cores N` and `-c N` remain aliases for `--threads N`. A missing, nonpositive,
or conflicting budget fails before workflow execution. The budget is shared by
all concurrent Snakemake jobs, not granted separately to every job. On Linux,
the launcher also restricts the inherited CPU affinity to at most N available
logical CPUs; this bounds CPU use by tools with additional internal helper
threads. If a scheduler already grants fewer CPUs, that smaller allocation is
used. Internal I/O/helper threads can exist beyond N; N is the CPU execution
budget, not a limit on the number of operating-system thread objects.

Fastp, HISAT2, StringTie, GeneMark-ETP and AUGUSTUS prediction receive the workflow
allocation. UTR ranking retains up to four workers per sequence, so several
sequences can run within the same total budget. Sequential Python stages use
one CPU. Nested BLAS/OpenMP pools are restricted to one thread per worker.
HISAT2 and samtools share their allocation: alignment uses N-1 workers while
the sorter uses one. With N=1, alignment writes a temporary SAM before sorting,
so the two tools run sequentially; this requires temporary disk space.

AUGUSTUS prediction streams the genome and hints into one job per complete
reference sequence, starts at most N AUGUSTUS processes, and uses the bundled
`join_aug_pred.pl` to merge outputs in reference order and renumber IDs.
Softmasking, hints coordinates and model parameters are preserved. Splitting
inside a reference sequence is deliberately avoided, so long genes acquire no
new artificial boundaries; a one-sequence input therefore uses one prediction
worker. Worker failures stop pending/running prediction jobs and prevent final
output publication. `work/prediction_parallel.json` records sequence jobs and
the observed maximum number of active workers; `work/logs/augustus_prediction.log`
records progress. Changing N does not retrain the model or change the
optimization cross-validation setting by itself.

`prediction.augustus_kfold` defaults to 8. Optimization uses at most the smaller
of the workflow thread allocation and the fold count. To use up to 34 workers,
set this value to 34 and run `geneinit-run --threads 34`. Increasing folds also
increases repeated model training and changes the cross-validation partition;
it does not enlarge the optimization or independent holdout sets. A compatible
best parameter file can be supplied as
`prediction.augustus_optimization_start_parameters` to retain accepted metadata
from a previous run. It must come from the same species, selected data and
AUGUSTUS version. It is applied after the original initial holdout evaluation,
and the current optimization evaluates it again under the configured folds.

```bash
# RNA evidence only
conda run --prefix .conda/bootstrap \
  snakemake rna_evidence --profile profiles/local \
  --configfile config/config.yaml

# RNA -> GeneMark-ETP -> full ETP hints -> direct AUGUSTUS training/prediction
conda run --prefix .conda/bootstrap \
  snakemake prediction --profile profiles/local \
  --configfile config/config.yaml

# Candidate evidence -> balanced base -> unified coding-model selection
conda run --prefix .conda/bootstrap \
  snakemake coding_selection --profile profiles/local \
  --configfile config/config.yaml

# Materialize and provisionally group selected coding models (no UTR)
conda run --prefix .conda/bootstrap \
  snakemake gene_grouping --profile profiles/local \
  --configfile config/config.yaml

# Transcript-first merge/split, bridge recovery, and connector retirement;
# iterate to no structural change, with policy maximum_iterations as a hard cap
conda run --prefix .conda/bootstrap \
  snakemake structural_postqc --profile profiles/local \
  --configfile config/config.yaml

# One final InGenAnnot rank/refine pass over every active post-QC transcript
# with exact CDS/phase preservation gates
conda run --prefix .conda/bootstrap \
  snakemake utr_isoforms --profile profiles/local \
  --configfile config/config.yaml

# Coordinate-stable IDs, GFF3/FASTA, lineage, QC, checksums, and manifest
conda run --prefix .conda/bootstrap \
  snakemake structural_release --profile profiles/local \
  --configfile config/config.yaml
```

Until the bundled policy status is promoted from `MIGRATION_PENDING` to
`PRODUCTION`, this target writes to `11_release.staging/` and marks the result
`PASS_MIGRATION_PENDING`; it does not authorize a public production release.

InGenAnnot 0.0.20 directly represents paired/single-end and reverse-stranded or
unstranded BAMs. Forward-stranded libraries are counted orientation-neutrally
at the final UTR-ranking step because that InGenAnnot version has no forward-
strand switch; this fallback is recorded in every shard status. HISAT2,
StringTie, and junction evidence remain configured with the declared forward
protocol.

See [Scale-independent execution](docs/scalability.md) for the automatic block
contract and its invariants.

The target must precede `--configfile` with the pinned Snakemake CLI.

## Current verification

CI uses pytest so that function-based and parameterized tests run alongside the
unittest classes. Install the bootstrap lock to reproduce the CI dependencies,
including pysam for BAM-input tests:

```bash
conda-lock install --prefix .conda/bootstrap envs/locks/linux-64/bootstrap.conda-lock.yml
```

```bash
conda run --prefix .conda/bootstrap python -m pytest -q

XDG_CACHE_HOME=/tmp/geneinit_snakemake_cache \
  conda run --prefix .conda/bootstrap snakemake \
  --configfile tests/fixtures/preflight/config_prediction.yaml --lint
```

Rule-by-rule migration state is tracked in
`docs/development/migration_matrix.tsv`.
The fixed-input regression scope and results are documented in
`docs/development/golden_regression.md`.

### Confirmed-only TE selection

The default production-v2 policy retains structurally supported uncertain TE
candidates in the general pool and separates only confirmed TE evidence.
Repeat overlap is report-only. See [the evidence contract](docs/te-evidence.md)
for the optional completed-evidence input, separate TE outputs, and the explicit
limits of running without calibrated family signatures.
# Evidence-guided recovery in masked regions

Set `prediction.masking_rescue: true` to add a local AUGUSTUS pass before the
candidate universe is built. The existing trained parameters and ETP hints are
reused. Retained protein chains and recurrent junctions from the supplied RNA
samples identify incomplete structures overlapping classified TE annotations.
Only the implicated TE intervals are unmasked, with 20 kb of surrounding context.
BUSCO labels and external validation RNA are not inputs to this step.

Recovered, complete ORFs enter the usual evidence scoring, confirmed-TE
classification, selection, structural QC, UTR refinement and release stages.
They remain AUGUSTUS candidates, so the second pass cannot inflate independent
predictor concordance. The original repeat annotations are preserved. Target
coordinates, supporting evidence IDs and candidate counts are recorded under
`05_prediction/04_masking_rescue/`. Candidate recovery does not establish gene
function or resolve TE capture versus repeat misclassification.

The option is off for existing configurations that omit it. Completed predictors
can be reused in an isolated downstream run; never allow upstream rules to write
through links to a previous run. This recovery pass does not retrain GeneMark or
AUGUSTUS and is not a replacement for reviewing the resulting gene structures.
