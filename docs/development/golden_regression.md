# Golden regression baseline

> These counts describe the frozen M2/M27 result produced before UTR
> reconstruction was moved behind iterative post-QC. They remain a comparison
> baseline, not validation of the reordered DAG. A new full M2/M27 regression
> must confirm final coding membership and measure UTR recovery for promoted
> transcripts.

GeneInIt was previously regressed on the fixed M2 and M27 nuclear assemblies, existing
GeneMark-ETP/AUGUSTUS predictions, RNA evidence, and coordinate-matched repeat
annotation. This validates candidate construction through deterministic
structural release. It does not yet validate a clean execution of the licensed
predictor from raw reads.

## Final counts

| Metric | M2 | M27 |
| --- | ---: | ---: |
| Candidate structures | 26,027 | 26,240 |
| Selected coding models | 16,847 | 16,892 |
| Provisional genes | 13,287 | 13,278 |
| Final genes | 13,219 | 13,208 |
| Final transcripts/CDS/proteins | 16,973 | 17,043 |
| Supported predicted-model bridges | 104 | 132 |
| Synthetic RNA-junction bridges | 22 | 19 |
| Adaptive training introns | 69,250 | 69,082 |
| Q99.5 training intron | 629 bp | 671 bp |
| Q99.9 training intron | 3,551 bp | 5,007 bp |
| Held bridge: additional evidence | 27 | 20 |
| Held bridge: extreme | 22 | 27 |
| Release-blocking QC flags | 0 | 0 |

For both assemblies, GFF3 transcript counts equal transcript, CDS, and protein
FASTA counts. All recorded release checksums pass. ORF, hierarchy, and sequence
consistency failure counts are zero.

## Adaptive synthetic RNA bridge guard

A splice junction alone can connect otherwise separate genes across an
unusually long genomic interval. Policy v2 learns thresholds from unique
introns inside selected complete-ORF models whose exon/CDS boundaries match an
exact, replicated RNA junction. The automatic threshold is
`max(5 kb, Q99.5)`, capped at 50 kb. The review threshold is
`max(4 * automatic, Q99.9)`, capped at 500 kb. Fewer than 2,000 training
introns triggers a recorded 5/20 kb conservative fallback.

Both Ganoderma assemblies learned an effective 5/20 kb threshold. Their final
annotation GFF3 SHA-256 values were identical to the frozen fixed-threshold
results. The observed retained synthetic intron maxima were 4,656 bp for M2
and 4,607 bp for M27. Longer candidates are written to
`09_postqc/<assembly>.bridge_reviews.tsv` and propagated into release
`review_flags.tsv`; they do not merge gene components automatically. A complete
candidate model already supported by the predictor is assessed by the separate
predicted-model bridge policy.

The adaptive post-QC component and direct predictor command contract contain no
fixed fungal mode. A clean Cuscuta raw-input regression is still required before
plant production status is granted.

## Release state

The regression result is `PASS_MIGRATION_PENDING`. Production status remains
blocked until the licensed GeneMark-ETP/AUGUSTUS deployment has completed a clean
raw-input-to-release run and the frozen production policy is explicitly
approved.
