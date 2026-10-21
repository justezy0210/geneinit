from pathlib import Path
from types import SimpleNamespace as NS

import pytest

from geneinit.evidence.components import RnaTranscript, RepeatFeature, JunctionAggregate, ProteinChain
from geneinit.models.candidate import CdsPart, parse_predictor_gtf, build_candidate_universe
from geneinit.predictor.masking import discover_targets, unmask_window, prepare_windows, collect_predictions
from geneinit.predictor.masking import optional_augustus_models


def fixture(strand='+'):
    rna=RnaTranscript('rna1','chr',strand,((100,150),(201,250),(301,350)))
    raw=NS(seqid='chr',strand=strand,start=201,end=350,cds_parts=(CdsPart(201,250,0),CdsPart(301,350,1)))
    te=RepeatFeature('chr',90,360,'DNA/Helitron')
    junctions={('chr',a,b,strand):JunctionAggregate(10,3,1,1) for a,b in rna.introns}
    return rna,raw,te,junctions


@pytest.mark.parametrize('strand',['+','-'])
def test_missing_supported_chain_is_discovered_on_either_strand(strand):
    rna,raw,te,j=fixture(strand)
    targets=discover_targets([raw],[rna],[],[te],j)
    assert len(targets)==1 and targets[0].repeats==((90,360),)
    assert targets[0].evidence==('RNA:rna1',)
    raw.cds_parts=(CdsPart(100,150,0),*raw.cds_parts)
    raw.start=100
    assert discover_targets([raw],[rna],[],[te],j)==[]


def test_weak_rna_and_unrelated_repeat_do_not_trigger():
    rna,raw,te,j=fixture()
    weak={k:JunctionAggregate(1,1,0,0) for k in j}
    assert discover_targets([raw],[rna],[],[te],weak)==[]
    assert discover_targets([raw],[rna],[],[RepeatFeature('chr',320,340,'DNA/Helitron')],j)==[]
    assert discover_targets([raw],[rna],[],[RepeatFeature('chr',90,360,'Unknown')],j)==[]


def test_complete_protein_chain_does_not_require_expression():
    rna,raw,te,j=fixture()
    protein=ProteinChain('direct','chr','+','query','target',tuple(CdsPart(a,b,0) for a,b in rna.exons),rna.introns,((100,102),),((348,350),))
    assert len(discover_targets([raw],[],[protein],[te],{}))==1


def test_unmasking_preserves_sequence_and_outside_case():
    changed,n=unmask_window('aaCCttggAAaa',2,11,[(5,8)])
    assert changed=='aCCTTGGAAa' and n==4
    with pytest.raises(RuntimeError):unmask_window('aaaa',1,3,[(2,4)])


def test_window_coordinates_hints_and_augustus_identity(tmp_path):
    rna,raw,te,j=fixture()
    target=discover_targets([raw],[rna],[],[te],j)[0]
    genome={'chr':'a'*500}
    hints=tmp_path/'h.gff';hints.write_text('chr\tb2h\tintron\t151\t200\t9\t+\t.\tsrc=E\n')
    windows,shards=prepare_windows([target],genome,hints,tmp_path/'windows',flank=10)
    w=windows[0];assert (w['start'],w['end'])==(80,370)
    assert '\t72\t121\t' in shards[0].hints.read_text()
    # A complete local ORF is lifted back, with the original predictor identity.
    sequence=list(genome['chr']);sequence[99:108]=list('atgaaataa');genome['chr']=''.join(sequence)
    shards[0].output.write_text('chr\tAUGUSTUS\tCDS\t21\t29\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";\n')
    output=tmp_path/'out.gtf'
    counts=collect_predictions(windows,shards,genome,output)
    assert counts['candidate_models']==1
    models=parse_predictor_gtf(output,'augustus',{k:v.upper() for k,v in genome.items()})
    assert models[0].start==100 and models[0].end==108 and models[0].orf.status=='complete'
    assert models[0].predictor=='augustus'
    universe=build_candidate_universe(models+models)
    assert {m.predictor for m in universe[0].models}=={'augustus'}


def test_clipped_orf_is_not_rescued(tmp_path):
    rna,raw,te,j=fixture()
    target=discover_targets([raw],[rna],[],[te],j)[0]
    hints=tmp_path/'h.gff';hints.write_text('')
    genome={'chr':'a'*79+'ATGAAATAA'+'a'*412}
    windows,shards=prepare_windows([target],genome,hints,tmp_path/'windows',flank=10)
    shards[0].output.write_text('chr\tAUGUSTUS\tCDS\t1\t9\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";\n')
    counts=collect_predictions(windows,shards,genome,tmp_path/'out.gtf')
    assert counts.get('candidate_models',0)==0


def test_successful_window_without_genes_is_an_empty_contribution(tmp_path):
    rna,raw,te,j=fixture()
    target=discover_targets([raw],[rna],[],[te],j)[0]
    hints=tmp_path/'h.gff';hints.write_text('')
    genome={'chr':'a'*500}
    windows,shards=prepare_windows([target],genome,hints,tmp_path/'windows',flank=10)
    shards[0].output.write_text('##gff-version 3\n# no predicted genes\n')
    output=tmp_path/'out.gtf'
    assert collect_predictions(windows,shards,genome,output).get('candidate_models',0)==0
    assert optional_augustus_models(output,genome)==[]


def test_no_conflicts_is_an_empty_contribution(tmp_path):
    output=tmp_path/'out.gtf'
    assert collect_predictions([],[],{'chr':'ACGT'},output)=={}
    assert optional_augustus_models(output,{'chr':'ACGT'})==[]


def test_malformed_prediction_is_not_silently_ignored(tmp_path):
    output=tmp_path/'bad.gtf';output.write_text('not a GTF record\n')
    with pytest.raises(RuntimeError):optional_augustus_models(output,{'chr':'ACGT'})
