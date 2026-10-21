# Unified selection v3 — experimental

This implementation consolidates admission and conflict resolution. It does
not claim identical membership to RNA50 v2, improved biological accuracy, or
a measured speedup. `production_v2.json` and the adopted annotations remain
unchanged. Select `policies/unified_v3.json` explicitly in a new run config.
The frozen policy marks the result `EXPERIMENTAL_NOT_VALIDATED`.

## Implemented workflow

1. Generate candidates and compute RNA/protein evidence as before.
2. Admit candidates using one structural gate and the union of core and
   supplemental evidence predicates; resolve conflicts once.
3. Inspect selected structures for RNA-supported repairs. Score each new
   structure against the original evidence and submit it to the same selector.
   Re-resolve the affected old and new competition loci, reusing unchanged
   evidence. Stop at a fixed point; fail at the existing iteration limit.
4. Assign biological gene containers after repair, then run the existing UTR,
   TE-release, FASTA, and structural-validation stages.

The former provisional grouping stage now writes one-transcript transport
containers for v3. It does not compute a second biological grouping graph.
Final grouping uses the existing coding-component definition.

## Deliberate decision changes

- The second balanced selector, its selected-set prerequisite, and
  grandfathered base-to-base conflicts are absent.
- Supplemental evidence predicates remain explicit admission routes.
  The weak single-exon support predicate participates in the common ranking;
  it no longer requires that the separate base selector found no main model.
- Core evidence retains its rank advantage. Additional isoform support is
  checked by one rule set, available to all candidates.
- Strict fusion precedence uses current RNA/protein evidence; it does not
  require a positive decision from a second selector. Fusion witnesses are
  evaluated from a fixed eligible snapshot.
- An existing rejected model cannot bypass selection merely by passing the
  old post-QC tier-2 bridge rule. It can re-enter when a changed local graph
  selects it under the same policy.
- A synthesized RNA model receives its own ORF, splice, junction, protein-chain
  and repeat metrics. No parental score or independent predictor vote is
  invented. Synthetic models never become new independent repair anchors.
- Known TE and general candidates remain separate competition pools.
- The 50/100-aa thresholds, ORF checks, RNA50 threshold, adaptive bridge
  lengths and final UTR contract remain in force. The 72-aa Doddi example
  therefore remains excluded until a separate length-policy change is tested.
- RNA hint copy weighting is unchanged; its biological effect is a separate
  question from this selector refactor.

## Cache and provenance

`UnifiedSelector.update` invalidates the complete old locus when a candidate
is removed or moved, and every existing locus connected to an upsert. It
resolves the resulting union, including previously rejected competitors.
Disconnected loci retain their decisions. A changed evidence row is
fingerprinted and reevaluated; unchanged rows retain their admission results.

The post-QC process restores original decisions from the completed audit,
without rerunning the initial selector. Policy, component-table and audit
SHA-256 checks and sample/replicate counts prevent incompatible restoration.
Changing global evidence sources still requires rebuilding their affected
component rows; this API is not a substitute for upstream invalidation.

Additional declared post-QC outputs:

- `*.unified_selection.tsv`: all candidate decisions in affected loci per
  local update.
- `*.unified_rounds.tsv`: new evaluations, reused evaluations, affected
  candidates and untouched candidates.

Legacy balanced-base output paths remain as empty compatibility files in v3,
with a `NOT_APPLICABLE` summary. Legacy runs emit empty unified-audit headers.
The ordinary audit, candidate IDs and release contracts remain available.
Existing iteration column `supported_bridge_promotions` counts restored
existing candidates in v3; the individual `status` values distinguish these
from synthetic RNA models. Counts of promotion/retirement events need not
equal the final active-model difference when later rounds replace a model.

## Validation and adoption

Unit and synthetic integration tests cover common hard QC, supplemental
admission, conflict resolution, TE isolation, deterministic fusion handling,
move/delete/evidence-change invalidation, checkpoint restoration, policy
mismatch rejection, a newly scored RNA repair, convergence, release generation,
and release checksum validation. External GeneMark/AUGUSTUS training and
InGenAnnot execution are not exercised by those tests.

The companion `replay_regression.py` replays frozen Jinju candidates through
selection and post-QC into a new directory. It records input/source hashes,
candidate-level changes, timings, and optional Doddi-candidate outcomes.
It does not rerun mapping or predictors, add UTRs, overwrite a release, build
a SIF, or establish biological superiority.

Adoption requires reviewing added/removed structures, the 6,170 formerly
supplemented models, the 390 old bridge recoveries and 196 synthetic bridges,
and known Doddi discrepancy loci. Full UTR/release validation on real data
and performance measurements follow that comparison. Final counts must be
measured rather than copied from v2.

## Controlled single-exon 50-aa comparison

`policies/unified_v3_single_exon50.json` lowers only
`selector_core.hard_qc.single_exon_min_protein_aa` from 100 to 50.
Policy identifiers and the note distinguish this experimental variant.
All other admission, ORF, graph, TE and post-QC rules are identical.
The original `unified_v3.json` remains the 100-aa control.

The companion `run_single_exon50.py` uses a completed 100-aa replay as its
comparison. It verifies matching upstream evidence and pipeline source,
then writes all 50–99-aa candidate decisions, newly retained short models,
excluded previous models and Doddi case changes. The common replay driver
accepts explicit `--policy` and `--comparison` arguments for this purpose.
Passing the new length floor does not bypass protein-evidence admission or
competition with another model. The 72-aa Doddi example can therefore still
fail its evidence route after passing structural QC.

The subsequent experimental boundary-ranking policy is
`policies/unified_v4_single_exon50_boundary.json`. It retains the 50-aa floor
and v3 admission routes, adds shared RNA/CDS boundary evidence, and uses it
to rank and automatically select competing candidates.
See [RNA/CDS boundary ranking](rna_cds_boundary.md)
for the evidence definitions, checkpoint behavior and validation limits.
