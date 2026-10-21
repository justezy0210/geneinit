"""Separate confirmed TE models and preserve all candidate dispositions."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

from geneinit.evidence.te_classification import TE_COLUMNS, model_identity, protein_sha256
from geneinit.models.candidate import read_fasta
from geneinit.postqc.core import Node, _candidate_model, build_edges, coding_components
from geneinit.utr.core import atomic_text, _write_table


def write_te_release(output_dir, candidates_path, audit_path, generated_path,
                     classifier, genome, general_models, gene_map, transcript_map,
                     lineage_rows, decisions, gene_prefix):
    from geneinit.release.core import _structural_validate, _write_gff_row, _write_fasta
    with audit_path.open() as handle:
        audits = {r['candidate_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    generated = json.loads(generated_path.read_text())['candidates'] if generated_path else {}
    overlap = set(audits) & set(generated)
    if overlap:
        raise ValueError(f"duplicate raw/generated candidates: {sorted(overlap)[:5]}")
    by_lineage = {r['candidate_id']: r for r in lineage_rows}
    current = {}
    for model in general_models.values():
        identity = {**audits.get(model.candidate_id, {}), **model_identity(model)}
        result = classifier.classify(identity, protein_sha256(_structural_validate(model, genome)[2]))
        if result['candidate_pool'] != 'GENERAL':
            raise ValueError(f"confirmed TE escaped into general release: {model.candidate_id}")
        previous = audits.get(model.candidate_id, generated.get(model.candidate_id))
        if previous is None or previous['protein_sha256'] != result['protein_sha256']:
            raise ValueError(f"release CDS lacks current classification: {model.candidate_id}")
        current[model.candidate_id] = result
    te_models = {}
    with candidates_path.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if audits[row['candidate_id']].get('te_selected') == '1':
                model = _candidate_model(row)
                te_models[model.candidate_id] = model
    for key, row in generated.items():
        if row['candidate_pool'] == 'TE':
            te_models[key] = _candidate_model(row)
    te_nodes = {model.transcript_id: Node(model=model, input_gene_id='.', source_status='TE_SELECTED',
                        selection_locus_id=audits.get(key, {}).get('locus_id', '.'),
                        routes=generated.get(key, {}).get('allowed_route', audits.get(key, {}).get('routes', '')),
                        external_support=True, rna_tier=0, protein_tier=0)
                for key, model in te_models.items()}
    components, _ = coding_components(te_nodes, build_edges(te_nodes), 3)
    # Coding edges use transcript IDs; release bookkeeping uses candidate IDs.
    components = {root: [te_nodes[tid].candidate_id for tid in members]
                  for root, members in components.items()}
    # The selector already resolved raw TE conflicts; identical generated CDS
    # are represented once, retaining the rejected candidate in the audit.
    sequences = {}
    te_ids = {}
    cross_refs = {}
    te_genes = 0
    signatures = set()
    flags = []
    with atomic_text(output_dir/'te_annotation.gff3') as handle:
        handle.write('##gff-version 3\n')
        for component in sorted(components.values(), key=lambda c: min((te_models[k].seqid, te_models[k].start, k) for k in c)):
            members = []
            for key in sorted(component):
                model = te_models[key]
                sig = (model.seqid, model.strand, tuple((f.start, f.end, f.phase) for f in model.cds))
                if sig in signatures:
                    continue
                signatures.add(sig)
                members.append(key)
            if not members:
                continue
            te_genes += 1
            gene = f'{gene_prefix}TE{te_genes:06d}'
            first = te_models[members[0]]
            _write_gff_row(handle, first.seqid, 'GeneInIt', 'gene',
                          min(te_models[k].start for k in members), max(te_models[k].end for k in members),
                          first.strand, '.', [('ID', gene), ('candidate_pool', 'TE')])
            for index, key in enumerate(members, 1):
                model = te_models[key]
                seqs = _structural_validate(model, genome)
                result = classifier.classify(model_identity(model), protein_sha256(seqs[2]))
                if result['candidate_pool'] != 'TE':
                    raise ValueError(f'TE export lacks confirmed evidence: {key}')
                current[key] = result
                tid = f'{gene}.t{index}'
                te_ids[key] = (gene, tid)
                sequences[tid] = seqs
                refs = sorted({gene_map[g.gene_id] for g in general_models.values()
                               if g.seqid == model.seqid and g.strand == model.strand
                               and g.start <= model.end and model.start <= g.end})
                cross_refs[key] = ','.join(refs)
                attrs = [('ID', tid), ('Parent', gene), ('candidate_id', key),
                         ('candidate_pool', 'TE'), ('te_class', result['effective_class'])]
                if refs:
                    attrs += [('general_gene_cross_references', ','.join(refs)), ('te_review_flags', 'MIXED_REVIEW')]
                    flags.extend(dict(flag_type='MIXED_REVIEW', severity='INFO', subject_a=tid,
                        subject_b=ref, seqid=model.seqid, start=model.start, end=model.end,
                        reason='TE and general transcripts share a genomic interval', release_blocking=0) for ref in refs)
                _write_gff_row(handle, model.seqid, 'GeneInIt', 'mRNA', model.start, model.end, model.strand, '.', attrs)
                for kind, features in [('exon', model.exons), ('CDS', model.cds)]:
                    for i, feature in enumerate(features, 1):
                        _write_gff_row(handle, model.seqid, 'GeneInIt', kind, feature.start, feature.end,
                            model.strand, feature.phase if kind == 'CDS' else '.',
                            [('ID', f'{tid}.{kind}{i}'), ('Parent', tid)])
    _write_fasta(output_dir/'te_proteins.faa', ((tid, seq[2]) for tid, seq in sequences.items()))
    _write_fasta(output_dir/'te_cds.fna', ((tid, seq[1]) for tid, seq in sequences.items()))
    columns = ('candidate_id', 'origin_stage', 'parent_candidate_ids', 'allowed_route',
               'seqid', 'strand', 'start', 'end', 'hard_qc_pass', 'eligible',
               'repeat_cds_fraction', 'all_repeat_cds_fraction', 'final_disposition',
               'blockers', 'final_gene_id', 'final_transcript_id', 'general_gene_cross_references') + TE_COLUMNS
    counts = Counter()
    class_counts = Counter()
    def rows():
        for key, old in list(audits.items()) + list(generated.items()):
            row = dict(old)
            row.update(current.get(key, {}))
            lineage = by_lineage.get(key, {})
            if key in te_ids:
                disposition = 'TE_SELECTED'
                gid, tid = te_ids[key]
            elif lineage.get('final_transcript_id', '.') != '.':
                disposition = 'GENERAL_SELECTED'
                gid, tid = lineage['final_gene_id'], lineage['final_transcript_id']
                class_counts[row['effective_class']] += 1
            else:
                gid = tid = '.'
                if lineage or old.get('selected') == '1':
                    disposition = 'POSTQC_REJECT'
                elif key in generated or old.get('eligible') == '1' or old.get('decision', '').startswith('REJECT_ACTIVE_FUSION'):
                    disposition = 'GRAPH_REJECT'
                elif old.get('hard_qc_pass') == '0':
                    disposition = 'QC_REJECT'
                else:
                    disposition = 'EVIDENCE_REJECT'
            row.update(candidate_id=key, origin_stage=old.get('origin_stage', 'CANDIDATE'),
                parent_candidate_ids=','.join(old.get('parent_candidate_ids', [])),
                allowed_route=old.get('allowed_route', old.get('routes', '')),
                final_disposition=disposition, final_gene_id=gid, final_transcript_id=tid,
                general_gene_cross_references=cross_refs.get(key, ''),
                te_selected=str(int(disposition == 'TE_SELECTED')))
            counts[disposition] += 1
            yield {col: row.get(col, '') for col in columns}
    _write_table(output_dir/'te_candidate_audit.tsv', columns, rows())
    summary = dict(policy_id=classifier.bundle.get('policy_id', 'geneinit_te_confirmed_only_v2'),
        evidence_sha256=classifier.digest, raw_candidate_count=len(audits), generated_candidate_count=len(generated),
        general_gene_count=len(gene_map), general_transcript_count=len(general_models),
        te_gene_count=te_genes, te_transcript_count=len(te_ids),
        dispositions=dict(counts), general_classes=dict(class_counts),
        note='TE counts are separate; cross-referenced loci must not be added as independent biological genes.')
    with atomic_text(output_dir/'te_summary.json') as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write('\n')
    return flags
