"""Regression tests for observed failure modes and independent search oracles."""
import random
import json

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

import accession_sources as ac
from opencrispr_designer import (
    scan_opencrispr_sites, design_from_segments, screen_local_reference,
    reverse_complement, guide_format_variants, GuideRNA,
)
from local_screening import (
    TargetLocus, screen_reference, exact_target_loci, MAX_PANEL_BP,
)
from sequence_sources import parse_multifasta, manual_record, segments_from_ncbi_record
from validation import validate_opencrispr_guide
from workflow import guide_key, export_bundle, clear_design_state, reset_panel_if_changed

SPACER = 'ACGTACGTACGTACGTACGT'
LOCUS = SPACER + 'AGG'


def guide():
    return scan_opencrispr_sites(LOCUS)[0]


def test_unconfirmed_exact_is_never_silently_excluded():
    r = screen_reference(SPACER, {'other': LOCUS}, 3)
    assert r.total_hits == 1 and r.specificity_score == 50
    assert not r.intended_site_excluded
    assert validate_opencrispr_guide(guide(), panel_screen=r).status == 'REVIEW'
    with pytest.raises(ValueError, match='explicit intended_locus'):
        screen_local_reference(SPACER, {'other': LOCUS}, exclude_one_exact=True)


def test_cap_cannot_change_total_score_or_counts():
    panel = {str(i): LOCUS for i in range(301)}
    r = screen_reference(SPACER, panel, intended_locus=TargetLocus('0', 1, '+'), max_hits=1)
    full = screen_reference(SPACER, panel, intended_locus=TargetLocus('0', 1, '+'), max_hits=1000)
    assert r.total_hits == full.total_hits == 300
    assert r.mismatch_counts == full.mismatch_counts == (300, 0, 0, 0)
    assert r.specificity_score == full.specificity_score == .33
    assert len(r.hits) == 1 and r.truncated
    assert validate_opencrispr_guide(guide(), panel_screen=r).status == 'FAIL'


def test_confirmed_single_locus_and_reverse_coordinates():
    reverse = reverse_complement(LOCUS)
    panel = {'reverse': 'AAA' + reverse}
    r = screen_reference(SPACER, panel, intended_locus=TargetLocus('reverse', 7, '-'))
    assert r.total_hits == 0 and r.specificity_score == 100
    assert exact_target_loci(SPACER, panel) == [TargetLocus('reverse', 7, '-')]
    assert validate_opencrispr_guide(guide(), panel_screen=r).specificity_status.startswith('PASS')


@pytest.mark.parametrize('target', [TargetLocus('absent',1,'+'), TargetLocus('x',0,'+'),
                                    TargetLocus('x',1,'?'), TargetLocus('x',2,'+')])
def test_bad_on_target_coordinates_are_rejected(target):
    with pytest.raises(ValueError, match='locus'):
        screen_reference(SPACER, {'x': LOCUS}, intended_locus=target)


def test_coordinate_cannot_exclude_mismatched_site():
    with pytest.raises(ValueError, match='does not exactly match'):
        screen_reference(SPACER, {'x': 'T' + LOCUS[1:]}, intended_locus=TargetLocus('x',1,'+'))


@pytest.mark.parametrize('raw', ['>x\nACGT\n>x\nAAAA', '>x one\nACGT\n>x two\nAAAA',
                                 '>x\n>y\nACGT', '>\nACGT', 'ACGT\n>x\nACGT'])
def test_malformed_fasta_is_rejected(raw):
    with pytest.raises(ValueError):
        parse_multifasta(raw)


@pytest.mark.parametrize('panel', [{}, {'x':''}, {'':'ACGT'}, {'x':'A'*(MAX_PANEL_BP+1)}])
def test_empty_or_oversized_panel_is_rejected(panel):
    with pytest.raises(ValueError):
        screen_reference(SPACER, panel)


@pytest.mark.parametrize('spacer', ['A'*19, 'N'*20, 'A'*21])
def test_ambiguous_or_wrong_length_query_is_rejected(spacer):
    with pytest.raises(ValueError):
        screen_reference(spacer, {'x':LOCUS})
    with pytest.raises(ValueError):
        guide_format_variants(spacer)


@pytest.mark.parametrize('kwargs', [{'max_hits':0},{'max_hits':5001},{'max_mismatches':-1},{'max_mismatches':5}])
def test_invalid_search_settings(kwargs):
    with pytest.raises(ValueError):
        screen_reference(SPACER, {'x':LOCUS}, **kwargs)


def test_no_searchable_sites_produces_no_score():
    r = screen_reference(SPACER, {'x':'N'*40})
    assert r.scanned_sites == 0 and r.ambiguous_bases == 40
    assert r.specificity_score is None
    assert validate_opencrispr_guide(guide(), panel_screen=r).status == 'REVIEW'


def test_ambiguity_and_narrow_radius_require_review():
    r = screen_reference(SPACER, {'x':LOCUS, 'gap':'N'*40}, intended_locus=TargetLocus('x',1,'+'))
    assert r.ambiguous_bases == 40
    assert validate_opencrispr_guide(guide(), panel_screen=r).status == 'REVIEW'
    r = screen_reference(SPACER, {'x':LOCUS}, 0, TargetLocus('x',1,'+'))
    assert validate_opencrispr_guide(guide(), panel_screen=r).status == 'REVIEW'


def test_hit_order_prefers_exact_over_earlier_near_match():
    r = screen_reference(SPACER, {'a':'T'+LOCUS[1:], 'b':LOCUS}, max_hits=1)
    assert r.hits[0].contig == 'b' and r.total_hits == 2


def test_independent_oracle_for_both_strands_and_mismatches():
    rng = random.Random(18)
    # Oracle reverses the entire contig and searches forward, independently of
    # the implementation's CCN indexing on the original strand.
    seq = ''.join(rng.choice('ACGT') for _ in range(1000))
    seq += LOCUS + 'A'*7 + reverse_complement('T'+LOCUS[1:])
    expected = set()
    for strand, oriented in [('+',seq), ('-',str(Seq(seq).reverse_complement()))]:
        for i in range(len(oriented)-22):
            if oriented[i+21:i+23] != 'GG':
                continue
            target = oriented[i:i+20]
            mm = sum(a != b for a,b in zip(SPACER,target))
            if mm <= 4:
                start = i+1 if strand == '+' else len(seq)-(i+20)+1
                expected.add((start,strand,target,mm))
    r = screen_reference(SPACER, {'test':seq},4,max_hits=5000)
    actual = {(h.start,h.strand,h.spacer,h.mismatches) for h in r.hits}
    assert actual == expected


def record(seq='A'*100, record_id='NM_TEST.1'):
    rec = SeqRecord(Seq(seq),id=record_id)
    rec.annotations['molecule_type']='DNA'
    return rec


def test_short_annotated_exons_are_not_rejoined():
    r = record(LOCUS)
    r.features = [SeqFeature(FeatureLocation(0,23), type='CDS'),
                  SeqFeature(FeatureLocation(0,10), type='exon'),
                  SeqFeature(FeatureLocation(10,23), type='exon')]
    with pytest.raises(ValueError, match='short exons are not joined'):
        segments_from_ncbi_record(r)


def test_missing_transcript_exons_and_multiple_cds_are_blocked():
    r = record()
    r.features=[SeqFeature(FeatureLocation(0,100),type='CDS')]
    with pytest.raises(ValueError,match='exon boundaries'):
        segments_from_ncbi_record(r)
    r.features *= 2
    with pytest.raises(ValueError,match='Multiple CDS'):
        segments_from_ncbi_record(r)


def test_genomic_compound_cds_excludes_intron_and_preserves_reverse_orientation():
    r = record('A'*30 + LOCUS + 'A'*30, 'NC_TEST.1')
    r.features=[SeqFeature(CompoundLocation([FeatureLocation(0,30,strand=-1),
                                             FeatureLocation(53,83,strand=-1)]),type='CDS')]
    segments, _ = segments_from_ncbi_record(r)
    assert [len(s) for _,s in segments] == [30,30]
    assert not design_from_segments('x',segments)


def test_exon_cds_intersections_use_parts_not_bounding_span():
    r=record('A'*100,'NC_TEST.1')
    r.features=[SeqFeature(CompoundLocation([FeatureLocation(0,25),FeatureLocation(75,100)]),type='CDS'),
                SeqFeature(FeatureLocation(0,100),type='exon')]
    segments,_=segments_from_ncbi_record(r)
    assert [len(s) for _,s in segments] == [25,25]


def test_unannotated_small_genomic_record_allowed_with_context_warning():
    segments,warnings=segments_from_ncbi_record(record(LOCUS,'NC_TEST.1'))
    assert segments == [('accession_sequence', LOCUS)]
    assert any('No CDS' in x for x in warnings)


def test_ensembl_requested_version_is_not_silently_replaced(monkeypatch):
    class Response:
        def json(self): return {'version': 5}
    monkeypatch.setattr(ac,'_get',lambda *a,**k:Response())
    with pytest.raises(ValueError,match='version 4, received 5'):
        ac.fetch_ensembl_accession('ENST00000001.4')


def test_design_limits_and_duplicate_segments():
    for segments,count,score in [([('x',LOCUS)],0,0),([('x',LOCUS)],5,-1),
                                  ([('x',LOCUS),('x',LOCUS)],5,0),([('x','A'*500001)],5,0)]:
        with pytest.raises(ValueError): design_from_segments('x',segments,count,score)


def test_native_format_for_non_g_spacer_is_available():
    formats={x.label:x.sequence for x in guide_format_variants(SPACER)}
    assert formats == {'gX19':'G'+SPACER[1:], 'gX20':'G'+SPACER, 'X20':SPACER}


def test_export_and_state_do_not_reuse_old_panel_or_target():
    rec=manual_record('>x\n'+LOCUS)
    g=guide()
    other=manual_record('>y\n'+LOCUS)
    assert guide_key(rec,g) != guide_key(other,g)
    state={'oc_record':rec, 'oc_guides':[g], 'oc_local::old':123, 'theme':'dark'}
    clear_design_state(state)
    assert state == {'theme':'dark'}
    state={}
    reset_panel_if_changed(state,'panel',3,None,250)
    state['oc_screens']['x']=123
    reset_panel_if_changed(state,'panel',3,None,250)
    assert state['oc_screens']
    reset_panel_if_changed(state,'changed',3,None,250)
    assert state['oc_screens'] == {}
    screen=screen_reference(SPACER,{'x':LOCUS},intended_locus=TargetLocus('x',1,'+'))
    settings={'min_score':35,'min_spec_review':50,'max_guides':30}
    rows,bundle=export_bundle(rec,[g],{guide_key(rec,g):screen},settings)
    assert rows[0]['Off-target hits'] == 0
    data=json.loads(json.dumps(bundle))
    assert data['guides'][0]['panel_screen']['panel_sha256'] == screen.panel_sha256
    assert data['settings'] == settings
    assert data['guides'][0]['expression_formats']
    rows,_=export_bundle(rec,[g],{},settings)
    assert rows[0]['MIT specificity'] is None


def test_validator_rejects_cross_guide_screen():
    r=screen_reference('T'+SPACER[1:],{'x':LOCUS})
    with pytest.raises(ValueError,match='different targeting spacer'):
        validate_opencrispr_guide(guide(),panel_screen=r)


def test_validator_handles_invalid_symbols_as_fail():
    g=guide()
    g.spacer='-'*20
    assert validate_opencrispr_guide(g).status == 'FAIL'
