# GeneInIt

GeneInIt is being built as a one-assembly-per-run eukaryotic structural gene
annotation workflow. Its defining contract is that a soft-masked nuclear FASTA
is the only coordinate system used throughout one run.

The current development release reaches deterministic structural release
staging after transcript-first whole-genome post-QC. It has no organism-name
branches: biological choices come only from the versioned policy and declared
sample/assembly metadata. Coding identity is
structure-derived. Provisional grouping and post-QC operate on coding-only
models; after post-QC converges, one final UTR stage processes every active
transcript and fails if any CDS coordinate or phase changes. Post-QC ignores
provisional gene containers, rebuilds same-frame coding components, recovers
supported candidate/RNA-junction bridges, retires unsupported articulation
transcripts, and iterates until a no-change fixed point or the policy-defined
maximum of five rounds is reached. The release stage assigns stable
coordinate-ordered gene/transcript IDs, validates GFF3/CDS/protein consistency,
and writes lineage, QC, checksums, and a manifest. Candidate evidence and model
selection are streamed through exact coordinate-overlap blocks, while post-QC
lockstep-streams the full candidate/evidence/audit tables and retains only its
working set. It does not yet authorize a public production annotation. The
predictor stage runs GeneMark-ETP and AUGUSTUS directly without BRAKER or a
fixed taxon mode.

## Structural workflow

```text
preflight -> RNA evidence -> GeneMark-ETP/AUGUSTUS -> evidence integration
          -> coding selection -> coding-only provisional gene grouping
          -> iterative structural post-QC -> one final UTR reconstruction
          -> final release
```

See [Required inputs](inputs.md) for the current contract.
See [Scale-independent execution](scalability.md) for the streaming contract.
See [Golden regression](development/golden_regression.md) for validation scope
and structural release counts.
