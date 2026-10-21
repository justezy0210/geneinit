# Scale-independent execution

GeneInIt does not select a code path from the organism name, genome size,
candidate count, or a fungus/plant label. The same biological policy is applied
to every input assembly.

## Exact automatic partition

Candidate tables are coordinate ordered. Evidence scoring and coding-model
selection consume one span-overlap block at a time. A block closes only when
the next candidate is on another reference sequence or starts after the
maximum end coordinate of every candidate in the current block.

This boundary is exact rather than heuristic: an interval-overlap edge cannot
cross it. Opposite strands remain in the same conservative span block so the
canonical input order is preserved. Small or sparse inputs naturally create
small blocks; dense inputs create larger blocks. No candidate-count threshold
or organism-specific block size is used.

Post-QC reads candidate, component-evidence, and selection-audit tables once in
lockstep. It verifies candidate identity and coordinate order, retains evidence
for active transcripts, and retains unselected rows only when they satisfy the
versioned bridge policy. Bridge-to-transcript lookup uses an interval index.
Connector detection builds the graph once per iteration and computes
support-separating articulation nodes in linear graph time. At most one is
retired per connected component in a round, then the next round re-evaluates
the new graph to prevent cascading over-removal.

Every complete candidate-table scan also uses an exact temporary disk-backed
ID set. Duplicate identifiers therefore fail closed even when repeated rows are
far apart, without keeping millions of identifiers in Python memory. Pairwise
scans iterate by index or iterator combinations and do not allocate a copied
suffix list for every left-hand model.

## Result invariants

- Blocked and single-table execution use identical locus and edge definitions.
- No gene, coding edge, or bridge decision is allowed to cross an invalid
  sequence/strand axis.
- Stable locus and post-QC gene IDs are content-derived, not shard-derived.
- Candidate/evidence/audit order or membership drift fails the run.
- Every promoted supported or synthetic strong-RNA bridge must link at least
  two actual released component roots under the configured coding-link rule;
  promotion counts are derived from those verified links.
- Active-fusion models cannot serve as strong-RNA bridge anchors, and a
  synthetic strong-RNA bridge cannot overlap or span one; these boundaries
  remain hard through final post-QC.
- Final UTR restoration and release must reproduce the post-QC gene hierarchy
  exactly.

The peak memory of candidate-heavy evidence and selection stages is therefore
driven by the largest real span-overlap block rather than the total candidate
table. Selected transcript models, reference sequences, and external evidence
indexes are still held for the stages that require random access; these are
algorithmic working sets, not taxon-specific exceptions.
