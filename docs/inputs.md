# Required inputs

GeneInIt accepts one nuclear eukaryotic assembly per run.

- Soft-masked genome FASTA
- Repeat GFF3 in the exact same coordinate system
- External protein evidence FASTA appropriate for the target taxon
- RNA sample table: raw FASTQ (default) or precomputed BAM + CSI
- Assembly representation and expected protein copy count
- Nuclear genetic code and public gene prefix
- Versioned structural policy; `production_v2.json` is the adaptive default

GeneInIt currently accepts nuclear genetic code `1` only. Other translation
tables fail preflight instead of being silently treated as standard code.

The `prediction` target additionally requires:

- readable external GeneMark-ETP/AUGUSTUS SIF
- exact SHA-256 of that SIF
- `genemark_license_acknowledged: true`
- an Apptainer executable
- an optional taxon-neutral `prediction.genemark_proc` run label (default:
  `protein_evidence`)

The image path is runtime-specific. Its SHA-256 is included in the semantic run
lock, while the private image and GeneMark key are never copied into the
repository.

The locked predictor stack is GeneMark-ETP 1.02 and AUGUSTUS 3.5.0. GeneInIt
invokes GeneMark-ETP directly, converts its non-HC evidence back to genome
coordinates with GeneMark's `format_back.pl`, and builds the complete ETP hint
set. It then calls AUGUSTUS' `new_species.pl`, `gff2gbSmallDNA.pl`, `etraining`,
`optimize_augustus.pl`, and `augustus` directly. BRAKER is neither executed nor
required in the container.

No `--fungus` option is added. Soft masking remains enabled because it describes
the input representation, not a taxon. Five optimization rounds and 1,000 bp
training flanks are defaults and may be changed under `prediction`. GeneMark
`--force` remains forbidden so the explicitly stranded, staged StringTie GTF is
retained.

Repeat GFF3 seqids must be a subset of genome seqids. A repeat-free contig does
not have to appear in the GFF3.

## RNA strandedness

| Layout | Allowed values |
|---|---|
| Paired-end | `RF`, `FR`, `unstranded` |
| Single-end | `R`, `F`, `unstranded` |

The precise HISAT2 and StringTie mapping is versioned in
`config/samples.schema.yaml` and is covered by protocol-specific tests before a
protocol is declared production-ready.

Assign each biological RNA library once. In FASTQ mode, a read file cannot be
reused in another sample, including through a symlink or hard link. Preflight
also rejects byte-identical FASTQ copies assigned to different samples using
the existing input hashes. Identical bytes in two distinct mate files within
one sample are permitted. File checks do not establish biological independence
or recognize every re-encoded copy; sample and replicate labels must describe
the actual libraries.

## Configuration identity and restarts

Preflight tracks the configuration contents and optional prediction evidence
files, so a configuration change refreshes its normalized configuration and
audit records on the next workflow run. The semantic configuration hash uses
schema version 2 and includes normalized prediction/training and RNA-junction
thresholds, as well as content hashes of optional resume, training-snapshot,
and optimization-start files. Omitted settings and explicit default settings
have the same identity. Runtime output paths remain excluded.

Hashes from semantic schema version 1 and 2 are not interchangeable, even when
the biological inputs are unchanged. Individual downstream rules retain their
own dependency and parameter checks. Keep previous release files when moving
an existing analysis to the updated software.

## Precomputed RNA BAM + CSI

Use one input mode for the entire run. Omitting `rna` is equivalent to
`rna.input_mode: fastq`; the existing seven-column FASTQ manifest is unchanged.
For BAM mode:

```yaml
rna:
  input_mode: bam
  reference_fasta: /absolute/path/reference_used_for_alignment.fa
samples: /absolute/path/rna_bam_samples.tsv
```

The tab-separated header must be exactly:

```text
sample	condition	replicate	layout	strandedness	bam	csi
```

See `config/samples.bam.schema.yaml` and `config/samples.bam.example.tsv`.
Paths are resolved relative to the working directory; absolute paths are
recommended. Each sample requires a coordinate-sorted, splice-aware RNA BAM
and its CSI index, including for small genomes. CRAM and BAI are not accepted.
CSI may have a separate source path. Assign each biological library once;
duplicate sample IDs, condition/replicate pairs, and BAM files are rejected.
Declare the original library's layout and strandedness using the table above.
Do not infer these labels from alignment rate or reuse one BAM as multiple
replicates. BAM tags needed by downstream transcript/junction tools (for example,
`XS` or `ts` for splice orientation) should be retained from the aligner.

`reference_fasta` declares which reference was used to create the BAMs. The
operator must establish this provenance from the mapping run. Preflight compares
its sequence IDs, order, and uppercase bases against `genome`; changes only in
soft-mask case, FASTA descriptions, or line wrapping are accepted. Actual base,
coordinate, or sequence-name changes require new alignments. This comparison
cannot prove the alignment history of a BAM whose header lacks reference hashes.
Reusing an alignment after a mask-only change retains its original mapping
choices; it does not remap reads under the new mask.

Preflight also checks BAM header/EOF integrity, the ordered reference dictionary,
CSI magic and indexed access (including beyond 512 Mb on long sequences), and
layout/sort order in up to 256 mapped primary records from the first 10,000
records. These are bounded checks, not an exhaustive validation of every record
or an independent strand inference. Whole-file SHA-256 locks cover both BAM
and CSI. The importer checks those locks again before staging.

Start in a new execution directory. Imported BAM and CSI are symlinked into
`03_alignment/{layout}/{sample}/`, with an `.import.json` audit beside each BAM.
Source files must remain available and unchanged, and accessible through container
binds, for the entire run. Source/output aliases are rejected before the DAG can
remove outputs. No original BAM is rewritten or reindexed. StringTie sample and
merged transcripts, junction support, and all subsequent evidence are generated
in the new run. No old GTF or junction table is imported. Alignment program
headers and declared library metadata are retained in the preflight audit.

BAM/CSI content and sample metadata contribute to the semantic run hash; moving
identical inputs to another path does not change it. Mask differences in the
actual prediction `genome` still change its whole-file content lock as usual.
