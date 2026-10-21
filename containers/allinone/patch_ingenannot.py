"""Patch InGenAnnot 0.0.20 GTF exon identity to distinguish opposite strands."""
from pathlib import Path
import hashlib
import json
import ingenannot

path = Path(ingenannot.__file__).parent / "utils/gene_builder.py"
original_sha256 = "f125e364f8dd2291d85b8846a99ab840abbdd4c7b7c0a6d1c2347c23655ff385"
old = "return Exon('exon:{}-{}-{}'.format(feat.seqid,feat.start,feat.end)"
new = "return Exon('exon:{}-{}-{}-{}'.format(feat.seqid,feat.start,feat.end,feat.strand)"
text = path.read_text()
digest = hashlib.sha256(path.read_bytes()).hexdigest()
if digest == original_sha256:
    assert text.count(old) == 1
    path.write_text(text.replace(old, new))
else:
    # Verify idempotence against the exact pinned upstream file, not a substring.
    assert text.count(new) == 1
    assert hashlib.sha256(text.replace(new, old).encode()).hexdigest() == original_sha256
report = {"patch": "gtf-exon-identity-includes-strand", "path": str(path),
          "upstream_sha256": original_sha256,
          "patched_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
(path.parent / "geneinit_strand_patch.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report), flush=True)
