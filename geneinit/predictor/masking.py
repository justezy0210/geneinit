"""Generate additional AUGUSTUS candidates where evidence conflicts with TE masking.

This is candidate discovery, not a gene/TE classifier. Only the conflicting TE
intervals are unmasked, the trained model is reused, and recovered models retain
the AUGUSTUS predictor label so that rerunning it is not independent evidence.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from geneinit.evidence.metrics import IntervalIndex, introns_from_blocks, merge_intervals, overlap_bp
from geneinit.evidence.repeats import is_classified_te
from geneinit.models.candidate import parse_predictor_gtf
from geneinit.predictor.common import PredictorError
from geneinit.predictor.parallel import PredictionShard


METHOD = "evidence_guided_local_unmask_v1"
FLANK = 20_000


@dataclass(frozen=True)
class RescueTarget:
    seqid: str
    strand: str
    start: int
    end: int
    repeats: tuple[tuple[int, int], ...]
    evidence: tuple[str, ...]

    @property
    def identifier(self):
        payload = json.dumps([self.seqid, self.strand, self.start, self.end, self.repeats], separators=(",", ":"))
        return "GMR-" + hashlib.sha256(payload.encode()).hexdigest()[:16]


def index(items, stranded=False):
    return IntervalIndex(items, seqid=lambda x: x.seqid, start=lambda x: x.start,
                         end=lambda x: x.end, strand=(lambda x: x.strand) if stranded else None)


def discover_targets(models, transcripts, chains, repeats, junctions):
    """Require a missing RNA chain or missing protein structure inside classified TE.

    RNA route: >=2 introns independently observed in >=2 existing RNA samples
    each (>=3 HQ reads total), with a missing junction touching a TE. Protein
    route: a complete retained chain or >=2 introns and >=150 supported CDS bp.
    Neither route depends on BUSCO labels or a particular reference annotation.
    """
    raw_index = index(models, True)
    repeat_index = index([r for r in repeats if is_classified_te(r.classification)])
    proposed = []
    for source, entries in [("RNA", transcripts), ("PROTEIN", chains)]:
        for entry in entries:
            coding = entry.exons if source == "RNA" else entry.intervals
            if not coding:
                continue
            evidence_introns = set(entry.introns)
            if source == "RNA":
                evidence_introns = {
                    (a, b) for a, b in evidence_introns
                    if (j := junctions.get((entry.seqid, a, b, entry.strand))) is not None
                    and j.sample_count_hq_ge1 >= 2 and j.total_hq_reads >= 3
                }
                if len(evidence_introns) < 2:
                    continue
            elif not (entry.complete or (len(evidence_introns) >= 2 and sum(b-a+1 for a,b in coding) >= 150)):
                continue
            local = list(raw_index.query(entry.seqid, entry.start, entry.end, entry.strand))
            covered = []
            for model in local:
                parts = [(p.start, p.end) for p in model.cds_parts]
                introns = set(introns_from_blocks(parts))
                if evidence_introns <= introns:
                    if source == "RNA" or overlap_bp(parts, coding) >= .95 * sum(b-a+1 for a,b in coding):
                        covered.append(model)
            if covered:
                continue
            matching = [r for r in repeat_index.query(entry.seqid, entry.start, entry.end)
                        if overlap_bp(coding, [(r.start, r.end)])]
            if source == "RNA":
                # A UTR repeat must not trigger unmasking for an unrelated CDS defect.
                seen = max((set(introns_from_blocks([(p.start,p.end) for p in m.cds_parts])) for m in local),
                           key=lambda chain: len(chain & evidence_introns), default=set())
                missing = evidence_introns - seen
                matching = [r for r in matching if any(
                    r.start <= a-1 <= r.end or r.start <= b+1 <= r.end for a,b in missing)]
            if not matching:
                continue
            spans = merge_intervals((r.start, r.end) for r in matching)
            evidence_id = entry.transcript_id if source == "RNA" else entry.chain_id
            proposed.append(RescueTarget(entry.seqid, entry.strand,
                min(entry.start, spans[0][0]), max(entry.end, spans[-1][1]),
                spans, (source + ":" + evidence_id,)))
    # Merge biological target spans, not padded windows: flanks must not join
    # nearby but separate genes into chromosome-sized unmasked regions.
    result = []
    for target in sorted(proposed, key=lambda t: (t.seqid,t.strand,t.start,t.end)):
        if result and (result[-1].seqid,result[-1].strand)==(target.seqid,target.strand) and target.start<=result[-1].end:
            prev=result.pop()
            result.append(RescueTarget(target.seqid,target.strand,prev.start,max(prev.end,target.end),
                merge_intervals(prev.repeats+target.repeats),tuple(sorted(set(prev.evidence+target.evidence)))))
        else:
            result.append(target)
    return result


def masked_fasta(path):
    """Read the original case; model-auditing FASTA readers intentionally uppercase."""
    seqs = {}
    name = None
    with Path(path).open() as handle:
        for raw in handle:
            if raw.startswith(">"):
                name = raw[1:].split()[0]
                if name in seqs:
                    raise PredictorError("duplicate genome sequence ID")
                seqs[name] = []
            elif raw.strip():
                if name is None:
                    raise PredictorError("sequence before FASTA header")
                seqs[name].append(raw.strip())
    return {k: "".join(v) for k,v in seqs.items()}


def unmask_window(sequence, start, end, repeats):
    """Coordinates are 1-based inclusive; preserve every base and outside case."""
    if not 1 <= start <= end <= len(sequence):
        raise PredictorError("invalid rescue window")
    original = sequence[start-1:end]
    changed = list(original)
    for left,right in merge_intervals(repeats):
        if not start <= left <= right <= end:
            raise PredictorError("repeat outside rescue window")
        changed[left-start:right-start+1] = original[left-start:right-start+1].upper()
    changed = "".join(changed)
    assert original.upper() == changed.upper()
    return changed, sum(a != b for a,b in zip(original,changed))


def prepare_windows(targets, genome, hints, directory, flank=FLANK):
    directory.mkdir(parents=True, exist_ok=True)
    windows = []
    shards = []
    for t in targets:
        start=max(1,t.start-flank);end=min(len(genome[t.seqid]),t.end+flank)
        seq,changed=unmask_window(genome[t.seqid],start,end,t.repeats)
        if not changed:
            continue
        prefix=directory/t.identifier
        shard=PredictionShard(t.identifier,len(seq),prefix.with_suffix('.fa'),
            prefix.with_suffix('.hints.gff'),prefix.with_suffix('.gff3'),prefix.with_suffix('.log'))
        shard.fasta.write_text(f'>{t.seqid}\n{seq}\n')
        shard.hints.write_text('')
        windows.append({'target_id':t.identifier,'seqid':t.seqid,'strand':t.strand,'start':start,'end':end,
            'target_start':t.start,'target_end':t.end,'repeats':t.repeats,'evidence':t.evidence,'changed_case_bases':changed})
        shards.append(shard)
    # Index windows to stream the hints once without holding all hints in memory.
    by_bin=defaultdict(list)
    for i,w in enumerate(windows):
        for block in range((w['start']-1)//100_000,(w['end']-1)//100_000+1):
            by_bin[w['seqid'],block].append(i)
    buffers=defaultdict(list)
    buffered=0
    for raw in Path(hints).open():
        if raw.startswith('#') or not raw.strip():continue
        f=raw.rstrip().split('\t')
        if len(f)!=9:raise PredictorError('invalid hints record')
        a,b=int(f[3]),int(f[4])
        for i in by_bin.get((f[0],(a-1)//100_000),[]):
            w=windows[i]
            if a>=w['start'] and b<=w['end']:
                out=f.copy();out[3]=str(a-w['start']+1);out[4]=str(b-w['start']+1)
                buffers[i].append('\t'.join(out)+'\n')
                buffered+=1
        if buffered>=50_000:
            for i,lines in buffers.items():
                with shards[i].hints.open('a') as handle:handle.writelines(lines)
            buffers.clear()
            buffered=0
    for i,lines in buffers.items():
        with shards[i].hints.open('a') as handle:handle.writelines(lines)
    return windows,shards


def collect_predictions(windows, shards, genome, output):
    """Lift local coordinates and retain only central, complete novel structures."""
    seen=set();counts=defaultdict(int)
    with Path(output).open('w') as out:
        out.write('# '+METHOD+'; predictor=augustus; downstream evidence gates required\n')
        for w,shard in zip(windows,shards,strict=True):
            local={w['seqid']:genome[w['seqid']][w['start']-1:w['end']].upper()}
            for model in optional_augustus_models(shard.output,local):
                counts['local_models']+=1
                if model.strand!=w['strand'] or model.start<=1 or model.end>=shard.length:
                    counts['outside_target_or_boundary']+=1;continue
                if model.orf.status!='complete':
                    counts['incomplete_orf']+=1;continue
                offset=w['start']-1
                parts=[(p.start+offset,p.end+offset,p.phase) for p in model.cds_parts]
                if not overlap_bp([(a,b) for a,b,p in parts],w['repeats']):
                    counts['outside_target_or_boundary']+=1;continue
                signature=(w['seqid'],w['strand'],tuple(parts))
                if signature in seen:continue
                seen.add(signature)
                tid=w['target_id']+'.'+model.source_transcript_id
                gid=w['target_id']+'.'+model.source_gene_id
                attrs=f'gene_id "{gid}"; transcript_id "{tid}";'
                out.write(f'{w["seqid"]}\tAUGUSTUS\ttranscript\t{parts[0][0]}\t{parts[-1][1]}\t.\t{w["strand"]}\t.\t{attrs}\n')
                for a,b,phase in parts:
                    out.write(f'{w["seqid"]}\tAUGUSTUS\tCDS\t{a}\t{b}\t.\t{w["strand"]}\t{phase}\t{attrs}\n')
                    out.write(f'{w["seqid"]}\tAUGUSTUS\texon\t{a}\t{b}\t.\t{w["strand"]}\t.\t{attrs}\n')
                counts['candidate_models']+=1
    return dict(counts)


def optional_augustus_models(path, genome):
    """A successful local run may contain zero genes; malformed input still fails."""
    with Path(path).open() as handle:
        for number,raw in enumerate(handle,1):
            if raw.startswith('#') or not raw.strip():continue
            fields=raw.rstrip().split('\t')
            if len(fields)!=9:
                raise PredictorError(f'invalid rescue GTF record at {path}:{number}')
            if fields[2]=='CDS':
                return parse_predictor_gtf(Path(path),'augustus',genome)
    return []
