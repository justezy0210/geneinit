"""Behavior regressions for inclusion, identity, pool isolation and release."""
import copy
import csv
import json
from pathlib import Path
import runpy

import pytest

from geneinit.evidence.te_classification import TeClassifier, TeEvidenceError, POLICY_ID, protein_sha256
from geneinit.selection.core import build_loci_and_edges, evaluate_and_select, audit_row
from geneinit.selection.balanced import select_balanced_base, evaluate_base_record
from geneinit.postqc.core import _eligible_bridge

helpers = runpy.run_path(str(Path(__file__).with_name('test_selection_core.py')))
make_record = helpers['make_record']
POLICY = {**copy.deepcopy(helpers['POLICY']), 'te_policy_id': POLICY_ID}


def record(key='C1', fraction='1'):
    r = make_record(key, [(10, 99, 0), (200, 289, 0)])
    r.row['protein_phase_comparable_parts'] = '0'
    r.row['repeat_cds_fraction'] = fraction
    return r


def bundle(row, digest):
    provenance = dict(query_sha256=digest, tool_version='test-1', database_sha256='b'*64,
                      options={'mode': 'direct'}, receipt='fixture-run-1')
    return dict(schema='geneinit_te_evidence_v2', assembly_id='TEST', signature_version='test-1',
        exception_version='test-1', signatures={'TEST_TE': dict(validated=True, validation_reference='fixture only',
            ordered=True, domains=[dict(id='RT', min_score=30, min_coverage=.8)],
            intact_features=['LEFT_BOUNDARY', 'RIGHT_BOUNDARY'])},
        proteins={digest: dict(direct=dict(**provenance, status='COMPLETE_HIT', hits=[
            dict(kind='DIRECT_DOMAIN', family='TEST_TE', domain='RT', score=80, coverage=.9, start_aa=5)]),
            host=dict(**provenance, status='COMPLETE_NO_HIT', hits=[]))},
        candidates={row['candidate_id']: dict(**{k: row[k] for k in ('seqid','strand','start','end','cds_signature')},
            protein_sha256=digest, structures=[dict(family='TEST_TE', status='INTACT', seqid=row['seqid'],
                strand=row['strand'], start=int(row['start'])-5, end=int(row['end'])+5,
                source_sha256='c'*64, structure_id='T1', features=['LEFT_BOUNDARY','RIGHT_BOUNDARY'])])})


def classifier(tmp_path, data):
    path = tmp_path/'te.json'
    path.write_text(json.dumps(data))
    return TeClassifier(path, 'TEST')


def test_missing_and_weak_evidence_are_general(tmp_path):
    r = record()
    digest = protein_sha256('MAAA')
    assert TeClassifier().classify(r.row, digest)['candidate_pool'] == 'GENERAL'
    data = bundle(r.row, digest)
    data['proteins'][digest]['direct']['hits'][0]['kind'] = 'TRANSFERRED_PFAM'
    assert classifier(tmp_path, data).classify(r.row, digest)['candidate_pool'] == 'GENERAL'


@pytest.mark.parametrize('change, expected', [
    ('intact', 'TE_INCOMPLETE_EVIDENCE'), ('host_pending', 'HOST_SCREEN_PENDING'),
    ('host_conflict', 'HOST_TE_CONFLICT'), ('signature', 'TE_SUSPECT'),
])
def test_incomplete_and_conflicting_te_stay_general(tmp_path, change, expected):
    r = record(); digest = protein_sha256('MAAA'); data = bundle(r.row, digest)
    if change == 'intact': data['candidates']['C1']['structures'] = []
    if change == 'host_pending': data['proteins'][digest]['host'] = {'status':'NOT_RUN'}
    if change == 'signature': data['signatures']['TEST_TE']['validated'] = False
    if change == 'host_conflict':
        data['proteins'][digest]['host'].update(status='COMPLETE_HIT', hits=[{'kind':'BUSCO','status':'Complete'}])
    result = classifier(tmp_path,data).classify(r.row,digest)
    assert result['candidate_pool'] == 'GENERAL'
    assert result['effective_class'] == expected


def test_complete_te_separates_even_with_strong_rna_and_zero_overlap(tmp_path):
    r=record(fraction='0'); digest=protein_sha256('MAAA')
    result=classifier(tmp_path,bundle(r.row,digest)).classify(r.row,digest)
    assert result['candidate_pool'] == 'TE'
    r.row.update(result)
    loci,edges=build_loci_and_edges('TEST',[r])
    assert evaluate_and_select(loci,edges,POLICY,sample_count=4,replicates_per_condition=2) == []
    assert r.row['te_selected'] == '1'


@pytest.mark.parametrize('domains,hits,ordered,expected', [
    (['A', 'B'], [('B', 10), ('A', 50), ('B', 100)], True, 'TE'),
    (['A', 'B'], [('B', 100), ('A', 50), ('B', 10)], True, 'TE'),
    (['RT', 'RT'], [('RT', 5), ('RT', 40)], True, 'TE'),
    (['RT', 'RT'], [('RT', 5), ('RT', 5)], True, 'GENERAL'),
    (['A', 'B'], [('A', 50), ('B', 10)], True, 'GENERAL'),
    (['A', 'B'], [('A', 50), ('B', 50)], True, 'GENERAL'),
    (['A', 'B'], [('A', 50)], True, 'GENERAL'),
    (['A', 'B'], [('A', 50), ('B', 10)], False, 'TE'),
    (['RT', 'RT'], [('RT', 5)], False, 'TE'),
])
def test_domain_architecture_uses_a_valid_ordered_hit_chain(tmp_path, domains, hits, ordered, expected):
    r = record()
    digest = protein_sha256('M' + 'A' * 199)
    data = bundle(r.row, digest)
    signature = data['signatures']['TEST_TE']
    signature['ordered'] = ordered
    signature['domains'] = [dict(id=domain, min_score=30, min_coverage=.8) for domain in domains]
    data['proteins'][digest]['direct']['hits'] = [
        dict(kind='DIRECT_DOMAIN', family='TEST_TE', domain=domain,
             score=80, coverage=.9, start_aa=position)
        for domain, position in hits
    ]
    result = classifier(tmp_path, data).classify(r.row, digest)
    assert result['candidate_pool'] == expected
    assert result['te_architecture_pass'] == str(int(expected == 'TE'))
    # Domain presence and an ordered architecture remain separate audit fields.
    assert result['te_direct_pass'] == str(int(set(domains) <= {domain for domain, _ in hits}))


@pytest.mark.parametrize('change', ['hash','coordinate','failed','empty','receipt'])
def test_technical_failure_never_becomes_no_hit(tmp_path,change):
    r=record(); digest=protein_sha256('MAAA'); data=bundle(r.row,digest)
    if change=='hash': data['proteins'][digest]['direct']['query_sha256']='d'*64
    if change=='coordinate': data['candidates']['C1']['end']='291'
    if change=='failed': data['proteins'][digest]['host']={'status':'FAILED'}
    if change=='empty': data['proteins'][digest]['host']['status']='COMPLETE_HIT'
    if change=='receipt': del data['proteins'][digest]['host']['receipt']
    with pytest.raises(TeEvidenceError): classifier(tmp_path,data).classify(r.row,digest)


def test_same_protein_different_locus_does_not_share_intact_context(tmp_path):
    r=record(); digest=protein_sha256('MAAA'); c=classifier(tmp_path,bundle(r.row,digest))
    assert c.classify(r.row,digest)['candidate_pool']=='TE'
    other={**r.row,'candidate_id':'C2','start':'1010','end':'1289','cds_signature':'1010-1099:0,1200-1289:0'}
    assert c.classify(other,digest)['candidate_pool']=='GENERAL'
    assert len(c.cache)==1


def select(records):
    loci,edges=build_loci_and_edges('TEST',records)
    base,_=select_balanced_base(loci,edges,POLICY)
    return evaluate_and_select(loci,edges,POLICY,sample_count=4,replicates_per_condition=2,base_selected_keys=base)


def test_overlap_and_review_flags_cannot_change_ranking():
    def run(fraction):
        a,b=record('A',fraction),record('B','0')
        a.row['te_review_flags']='HOST_TE_CONFLICT,UNRESOLVED'
        return [r.key for r in select([a,b])]
    assert run('0') == run('1') == ['A']


def test_te_and_general_graphs_do_not_block_each_other(tmp_path):
    a,b=record('A'),record('B')
    digest=protein_sha256('MAAA')
    c=classifier(tmp_path,bundle(a.row,digest))
    a.row.update(c.classify(a.row,digest)); b.row.update(c.classify(b.row,digest))
    assert [r.key for r in select([a,b])] == ['B']
    assert a.row['te_selected']=='1'
    assert not b.blockers


def test_high_overlap_base_and_bridge_keep_non_te_requirements():
    r=record(); r.row['te_policy_id']=POLICY_ID
    assert evaluate_base_record(r).eligible
    evidence={**r.row,'candidate_id':'C1'}
    audit={'candidate_id':'C1','selected':'0','hard_qc_pass':'1','active_fusion':'0',
           'te_policy_id':POLICY_ID,'candidate_pool':'GENERAL'}
    bridge={'te_policy_id':POLICY_ID,'exclude_active_fusion':True,'minimum_rna_tier':2,'minimum_protein_tier':2}
    assert _eligible_bridge(r.row,evidence,audit,bridge)
    assert not _eligible_bridge(r.row,evidence,{**audit,'candidate_pool':'TE'},bridge)
    assert not _eligible_bridge(r.row,{**evidence,'orf_status':'partial'},audit,bridge)
    r.row['common_qc_pass']='0'
    assert select([r]) == []


def test_release_contains_uncertain_general_and_separate_confirmed_te(tmp_path):
    from geneinit.release.core import build_structural_release
    helpers=runpy.run_path(str(Path(__file__).with_name('test_release_core.py')))
    inputs=helpers['StructuralReleaseTests']()._write_inputs(tmp_path)
    r=record('C1'); r.row.update(start='10',end='99',cds_signature='10-99:0')
    # Identical proteins at different loci: only C2 has an intact TE context.
    digest=protein_sha256('M'+'A'*28)
    te={**r.row,'candidate_id':'C2','start':'110','end':'199','cds_signature':'110-199:0'}
    data=bundle(te,digest); c=classifier(tmp_path,data)
    sequence='C'*9+'ATG'+'GCC'*28+'TAA'+'C'*10+'ATG'+'GCC'*28+'TAA'+'C'*10
    inputs['genome'].write_text('>chr1\n'+sequence+'\n')
    candidates=tmp_path/'candidates.tsv'; audit=tmp_path/'selection.tsv'
    rows=[]
    for row in [r.row,te]:
        result=c.classify(row,digest)
        rows.append({**row,**result,'selected':'1' if row['candidate_id']=='C1' else '0',
                     'te_selected':'1' if row['candidate_id']=='C2' else '0','hard_qc_pass':'1','eligible':'1'})
    for path, values in [(candidates,[r.row,te]),(audit,rows)]:
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(values[0]),delimiter='\t');w.writeheader();w.writerows(values)
    gen=tmp_path/'generated.json';gen.write_text('{"candidates":{}}')
    inputs['policy'].write_text(json.dumps({'policy_id':POLICY_ID,'selector_core':POLICY,'status':'MIGRATION_PENDING'}))
    output=tmp_path/'release'
    names=['annotation.gff3','transcripts.fna','cds.fna','proteins.faa','longest_proteins.faa',
           'gene_id_map.tsv','transcript_id_map.tsv','candidate_lineage.tsv','structural_decisions.tsv','review_flags.tsv',
           'structural_qc.tsv','structural_qc.json','release_status.json','release_manifest.json','checksums.sha256','policy.json','release.complete.tsv']
    build_structural_release(inputs['gff3'],inputs['decisions'],inputs['edges'],inputs['bridge_reviews'],inputs['genome'],inputs['policy'],
        ['chr1'],'TEST','TESTG',*[output/name for name in names],tmp_path/'te.json',candidates,audit,gen)
    assert 'TE_INCOMPLETE_EVIDENCE' in (output/'annotation.gff3').read_text()
    assert 'CLASSIFIED_TE_OVERLAP' in (output/'annotation.gff3').read_text()
    assert 'candidate_id=C2' in (output/'te_annotation.gff3').read_text()
    summary=json.loads((output/'te_summary.json').read_text())
    assert summary['dispositions']=={'GENERAL_SELECTED':1,'TE_SELECTED':1}
    assert summary['te_gene_count']==1
    qc=json.loads((output/'structural_qc.json').read_text())
    assert qc['release_blocking_flag_count']==0


def test_release_groups_two_selected_complete_te_isoforms(tmp_path):
    from geneinit.models.candidate import audit_orf, read_fasta
    from geneinit.release.core import build_structural_release
    from geneinit.utr.core import load_raw_gff3

    release_helpers = runpy.run_path(str(Path(__file__).with_name('test_release_core.py')))
    inputs = release_helpers['StructuralReleaseTests']()._write_inputs(tmp_path)
    sequence = list(read_fasta(inputs['genome'])['chr1'] + 'C' * 520)
    sequence[209:299] = 'ATG' + 'GCC' * 29
    sequence[399:489] = 'GCC' * 30
    sequence[499:559] = 'GCC' * 20
    sequence[599:689] = 'GCC' * 29 + 'TAA'
    for start, end in [(300, 399), (490, 599), (460, 499), (560, 599)]:
        sequence[start - 1:start + 1] = 'GT'
        sequence[end - 2:end] = 'AG'
    genome = {'chr1': ''.join(sequence)}
    inputs['genome'].write_text('>chr1\n' + genome['chr1'] + '\n')
    general = make_record('C1', [(10, 99, 0)])
    alternatives = [
        make_record('GCAN-A', [(210, 299, 0), (400, 489, 0), (600, 689, 0)]),
        make_record('GCAN-B', [(210, 299, 0), (400, 459, 0), (500, 559, 0), (600, 689, 0)]),
    ]
    proteins, contexts, expected_sequences = {}, {}, {}
    for candidate in alternatives:
        orf = audit_orf('chr1', '+', candidate.cds, genome)
        assert orf.status == 'complete'
        protein = orf.protein_sequence.removesuffix('*')
        candidate.row['translated_protein_aa'] = str(len(protein))
        digest = protein_sha256(protein)
        data = bundle(candidate.row, digest)
        proteins.update(data['proteins'])
        contexts.update(data['candidates'])
        expected_sequences[candidate.key] = (orf.coding_sequence, protein)
    data['proteins'], data['candidates'] = proteins, contexts
    c = classifier(tmp_path, data)
    for candidate in [general, *alternatives]:
        orf = audit_orf('chr1', '+', candidate.cds, genome)
        candidate.row.update(c.classify(candidate.row, protein_sha256(orf.protein_sequence)))
    loci, edges = build_loci_and_edges('TEST', alternatives)
    assert [edge.relation for edge in edges] == ['ALT_SPLICE']
    assert evaluate_and_select(loci, edges, POLICY, sample_count=2, replicates_per_condition=2) == []
    assert all(candidate.row['te_selected'] == '1' for candidate in alternatives)

    candidates, audit = tmp_path/'candidates.tsv', tmp_path/'selection.tsv'
    general_audit = {**audit_row(general, POLICY['policy_id']),
                     'selected': '1', 'hard_qc_pass': '1', 'eligible': '1'}
    for path, rows in [
        (candidates, [candidate.row for candidate in [general, *alternatives]]),
        (audit, [general_audit] + [audit_row(candidate, POLICY['policy_id']) for candidate in alternatives]),
    ]:
        with path.open('w') as handle:
            columns = list(dict.fromkeys(key for row in rows for key in row))
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    generated = tmp_path/'generated.json'
    generated.write_text('{"candidates":{}}')
    inputs['policy'].write_text(json.dumps({'policy_id': POLICY_ID, 'selector_core': POLICY,
                                          'status': 'MIGRATION_PENDING'}))
    output = tmp_path/'release'
    names = ['annotation.gff3', 'transcripts.fna', 'cds.fna', 'proteins.faa', 'longest_proteins.faa',
             'gene_id_map.tsv', 'transcript_id_map.tsv', 'candidate_lineage.tsv', 'structural_decisions.tsv',
             'review_flags.tsv', 'structural_qc.tsv', 'structural_qc.json', 'release_status.json',
             'release_manifest.json', 'checksums.sha256', 'policy.json', 'release.complete.tsv']
    build_structural_release(
        inputs['gff3'], inputs['decisions'], inputs['edges'], inputs['bridge_reviews'],
        inputs['genome'], inputs['policy'], ['chr1'], 'TEST', 'TESTG',
        *[output/name for name in names], tmp_path/'te.json', candidates, audit, generated,
    )
    _, models = load_raw_gff3(output/'te_annotation.gff3')
    assert {model.candidate_id for model in models.values()} == {'GCAN-A', 'GCAN-B'}
    assert {model.gene_id for model in models.values()} == {'TESTGTE000001'}
    exported_cds = read_fasta(output/'te_cds.fna')
    exported_proteins = read_fasta(output/'te_proteins.faa')
    for transcript_id, model in models.items():
        cds, protein = expected_sequences[model.candidate_id]
        assert exported_cds[transcript_id] == cds
        assert exported_proteins[transcript_id] == protein
    with (output/'te_candidate_audit.tsv').open() as handle:
        final_audits = {row['candidate_id']: row for row in csv.DictReader(handle, delimiter='\t')}
    assert all(final_audits[candidate.key]['final_disposition'] == 'TE_SELECTED' for candidate in alternatives)
    assert all(final_audits[candidate.key]['final_gene_id'] == 'TESTGTE000001' for candidate in alternatives)
    summary = json.loads((output/'te_summary.json').read_text())
    assert summary['te_gene_count'] == 1 and summary['te_transcript_count'] == 2
    assert summary['dispositions'] == {'GENERAL_SELECTED': 1, 'TE_SELECTED': 2}
    assert 'PASS_MIGRATION_PENDING' in (output/'release.complete.tsv').read_text()


def test_reviews_are_advisory_and_invalidate_on_changed_evidence(tmp_path):
    from geneinit.evidence.te_classification import review_snapshot
    r=record(); digest=protein_sha256('MAAA'); data=bundle(r.row,digest)
    context=data['candidates']['C1']
    context['review']={'result':'UNRESOLVED','reason':'needs inspection',
                       'evidence_snapshot_sha256':review_snapshot(data,context)}
    assert classifier(tmp_path,data).classify(r.row,digest)['candidate_pool']=='GENERAL'
    context['structures'][0]['structure_id']='changed-evidence'
    result=classifier(tmp_path,data).classify(r.row,digest)
    assert 'STALE_REVIEW_IGNORED' in result['te_review_flags']
    assert result['review_result']=='NOT_REVIEWED'


def test_new_rna_bridge_is_classified_before_becoming_general_anchor():
    from geneinit.postqc.core import recover_strong_rna_bridges
    from geneinit.evidence.components import JunctionAggregate
    h=runpy.run_path(str(Path(__file__).with_name('test_postqc_core.py')))
    sequence=list('C'*300); sequence[99:102]='ATG';sequence[248:251]='TAA'
    junctions={('chr1',151,200,'+'):JunctionAggregate(total_hq_reads=5,sample_count_hq_ge1=2,
        condition_count_ge2rep_hq_ge1=1,condition_count_ge2rep_hq_ge3=1)}
    seen=[]
    def classify(model,left,right):
        seen.append((model.candidate_id,left,right));return False
    nodes={'left':h['node']('left','C1',100,150,gene_id='g1'),
           'right':h['node']('right','C2',201,251,gene_id='g2')}
    assert recover_strong_rna_bridges(nodes,junctions,{'chr1':''.join(sequence)},classify_generated=classify)==[]
    assert len(nodes)==2 and len(seen)==1 and seen[0][1:]==('C1','C2')
    assert seen[0][0].startswith('GRNA-')
