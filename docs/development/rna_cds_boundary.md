# RNA/CDS boundary ranking — experimental v4 rank1

Use `policies/unified_v4_single_exon50_boundary.json` explicitly in a new run.
The single-coding-exon minimum remains 50 aa. Admission evidence, ORF checks,
TE classification and fusion rules retain their v3 definitions. Boundary
evidence now informs automatic selection. The best eligible compatible model
set is returned with evidence and uncertainty recorded in the audit.

## Shared evidence calculation

`geneinit/evidence/rna_cds.py` compares each candidate against all overlapping
same-strand merged RNA transcripts. Both initial component scoring and
`RepairEvidence.score` call it through `best_rna_metrics`.

The comparison uses whole RNA introns, including those crossing the first or
last CDS boundary. It checks CDS/exon overlap, exact candidate intron inclusion,
and start/stop codon coordinates in transcription order, including split
codons. Introns entirely outside the CDS are UTR evidence, not CDS conflicts.
A complete compatible RNA contains the entire CDS chain, so the candidate's
existing genomic ORF audit supplies phase/start/stop validity. This calculation
does not discover new ORFs or demonstrate translation.

A conflict requires:

- An overlapping same-strand RNA transcript with some exonic CDS overlap.
- An RNA intron overlapping one or more candidate CDS bases.
- At least 3 HQ junction alignments, at least 2 supporting samples, and at least
  one condition with two supporting samples, using existing junction aggregates.

Missing, partial, opposite-strand and strand-unresolved RNA do not by themselves
produce a conflict. All overlapping isoforms are considered; the highest-ranked
transcript alone does not decide the outcome.

The component and selection audit tables preserve:

| Status | Meaning |
| --- | --- |
| SUPPORTED | At least one complete compatible RNA; no qualifying conflict |
| MIXED | Complete compatible RNA and qualifying conflict both exist |
| CONFLICT | Qualifying conflict without a complete compatible RNA |
| UNRESOLVED | Neither complete support nor qualifying conflict |

Transcript IDs, conflicting junction coordinates/counts, CDS overlap in those
introns, and start/stop conflict flags accompany the status.
These are structural observations, not independent truth labels.

## Automatic ranking without a preferred gene length

All candidates first pass the existing admission rules. Within each competition
pool, the graph uses this lexicographic priority (earlier items take precedence):

1. Candidates without a CONFLICT status precede candidates with that status.
   An actual recurrent contradiction is distinguished from absent/partial RNA.
2. Preserve the existing evidence-route class, including strong RNA, protein,
   predictor and complementary support routes.
3. Within the same route class, prefer complete RNA/CDS structural compatibility
   (SUPPORTED or MIXED) over missing/incomplete support.
4. Use the remaining existing rank fields, including independent protein
   metrics and predictor support. Length remains a late tie-breaker.

These are explicit heuristic priorities, not a fitted probability of accuracy.
No weighted sum, gene-specific exception or rule preferring a particular
protein length is introduced. Full RNA coverage of a short CDS does not
automatically outrank a better evidence-route class. Protein-chain coverage is
coverage of the retained alignment chain, not necessarily coverage of the
whole reference protein.

The graph selects compatible models normally; no candidate or whole locus is
held solely by the new rule. When all alternatives have a conflict, the best
eligible candidate can still be selected and receives the reason
`SELECTED_WITH_UNRESOLVED_RNA_CDS_CONFLICT`. Existing hard-QC and admission
failures are never bypassed.

The output audit records `rna_cds_conflict_penalty`,
`rna_cds_complete_support`, `selection_rank_json`, support/conflict transcript
IDs, junction counts, original evidence routes and competing blocker IDs.
Rank tuples sort ascending and are not confidence probabilities. Unselected
candidates remain in the universe and audit.

In the frozen Jinju fixture, 434 aa is now selected. The 301-aa model has a
recurrent terminal conflict; the 434-aa model has no detected contradiction
under this test and has a stronger existing route class than the 84-aa model.
The 84-aa model retains its complete RNA structural-support annotation but
loses competition with 434 aa. The unsupported connecting junction in 434 aa
remains recorded: automatic selection does not prove its biological correctness.

## Reuse and validation

New component evidence is mandatory: v4 fails on missing/stale boundary columns
instead of silently interpreting old tables as negative evidence. Component
completion schema is 3, with boundary evidence version 1. Existing v3 policies
remain usable and retain their selection rules.

Checkpoint restoration reconstructs the same ranking state. Local graph
updates rerank the affected locus, and structural repairs receive their own
fresh boundary metrics. The superseded review-only policy key is rejected
explicitly, and rank1 uses a new policy ID to prevent incompatible restoration.

Tests cover both strands, split codons, cross-boundary introns, UTR splices,
weak/no/partial RNA, compatible alternative isoforms, old-table rejection,
checkpoint restoration, and newly reconstructed models. Counterexamples check
that short or long models win according to their evidence, that a stronger
protein route can beat a complete short RNA, and that all-conflicting inputs
still produce a best eligible candidate with a flag. The frozen Jinju
fixture includes all 17 original candidates from the nine reviewed competition
loci, their 58 regional merged transcripts and 201 regional HQ junctions.
It reproduces the original local decisions before applying the new review rule.

## Limits and rollout

The automatic check uses the existing merged GTF and junction aggregates.
Sample-level RNA isoforms absent from the merge can change the rank later;
their absence from this check is not evidence that they do not exist.
An assembled transcript is not a measured full-length molecule. A shared
same-strand locus can contain alternative genes or isoforms.

This implementation ranks existing candidates and flags conflicts. It does not repair
the 301-aa endpoint, promote its incomplete 372-aa projection, classify the
eight TE-associated models as confirmed TE, or change the 150-aa host criterion.

The local fixture and synthetic workflow/release tests are regression checks.
A genome-wide before/after replay and inspection of changed loci are still needed
before adopting v4 as the production policy. Mapping, predictor training, full
UTR generation and the adopted annotations are outside these tests.
