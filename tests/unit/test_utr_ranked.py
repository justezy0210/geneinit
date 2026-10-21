from pathlib import Path

import pytest

from geneinit.utr.core import UtrError
from geneinit.utr.ranked import validate_ranked_rna


def records(tid, strand):
    attrs = f'gene_id "g{tid}"; transcript_id "{tid}";'
    return (f'chr1\tRNA\ttranscript\t10\t30\t.\t{strand}\t.\t{attrs} rank "1";\n'
            f'chr1\tRNA\texon\t10\t30\t.\t{strand}\t.\t{attrs}\n')


def setup_inputs(tmp_path):
    rna, ranked = tmp_path/'rna.gtf', tmp_path/'ranked.gtf'
    rna.write_text(records('minus', '-') + records('plus', '+'))
    ranked.write_text(rna.read_text().replace('exon\t10\t30\t.\t+', 'exon\t10\t30\t.\t-'))
    return ranked, rna, tmp_path/'validated.gtf', tmp_path/'audit.tsv'


def test_source_proven_repair_preserves_ranks_ids_coordinates_and_input(tmp_path):
    ranked, rna, output, audit = setup_inputs(tmp_path)
    original = ranked.read_bytes()
    assert validate_ranked_rna(ranked, rna, output, audit)["exon_strands_restored"] == 1
    assert output.read_bytes() == rna.read_bytes()
    assert ranked.read_bytes() == original
    assert audit.read_text().splitlines()[1] == 'plus\tchr1\t10\t30\t-\t+\tRESTORE_SOURCE_EXON_STRAND'
    second = tmp_path/'second.gtf'
    assert validate_ranked_rna(output, rna, second, audit)["exon_strands_restored"] == 0
    assert second.read_bytes() == output.read_bytes()


@pytest.mark.parametrize('change', ['transcript_strand', 'exon_coordinate', 'sequence',
    'unknown_transcript', 'missing_exon', 'duplicate_exon', 'source_strand', 'source_coordinate'])
def test_refuse_unproven_changes(tmp_path, change):
    ranked, rna, output, audit = setup_inputs(tmp_path)
    text = ranked.read_text()
    if change == 'transcript_strand':
        text = text.replace('transcript\t10\t30\t.\t+', 'transcript\t10\t30\t.\t-')
    elif change == 'exon_coordinate':
        text = text.replace('exon\t10\t30', 'exon\t11\t30')
    elif change == 'sequence':
        text = text.replace('chr1', 'chr2')
    elif change == 'unknown_transcript':
        text = text.replace('"plus"', '"missing"')
    elif change == 'missing_exon':
        text = ''.join(text.splitlines(True)[:-1])
    elif change == 'duplicate_exon':
        text += text.splitlines(True)[-1]
    elif change == 'source_strand':
        rna.write_text(rna.read_text().replace('exon\t10\t30\t.\t+', 'exon\t10\t30\t.\t-'))
    else:
        rna.write_text(rna.read_text().replace('exon\t10\t30', 'exon\t11\t30'))
    ranked.write_text(text)
    with pytest.raises(UtrError):
        validate_ranked_rna(ranked, rna, output, audit)
    assert not output.exists()


def test_subset_of_ranked_transcripts_is_valid(tmp_path):
    ranked, rna, output, audit = setup_inputs(tmp_path)
    ranked.write_text(records('minus', '-'))
    assert validate_ranked_rna(ranked, rna, output, audit)["exon_strands_restored"] == 0
    assert output.read_bytes() == ranked.read_bytes()


def test_unknown_source_strand_is_not_silently_oriented(tmp_path):
    ranked, rna, output, audit = setup_inputs(tmp_path)
    rna.write_text(records('unknown', '.'))
    ranked.write_text(records('unknown', '+'))
    result = validate_ranked_rna(ranked, rna, output, audit)
    assert result == {'exon_strands_restored': 0, 'unstranded_transcripts_excluded': 1}
    assert output.read_text() == ''
    assert 'EXCLUDE_UNKNOWN_SOURCE_STRAND' in audit.read_text()
