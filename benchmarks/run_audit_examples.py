"""Run deterministic software examples; no experimental efficacy is measured."""
from pathlib import Path
import sys
import json
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from opencrispr_designer import scan_opencrispr_sites, reverse_complement, guide_format_variants
from sequence_sources import parse_multifasta, segments_from_ncbi_record
from local_screening import screen_reference, TargetLocus
from workflow import reset_panel_if_changed

S='ACGTACGTACGTACGTACGT'
L=S+'AGG'
rows=[]
def add(name, actual, expected):
    assert actual == expected, (name,actual,expected)
    rows.append({'example':name,'observed':actual,'expected':expected,'status':'PASS'})
def rejected(name,fn,part):
    try: fn()
    except ValueError as exc:
        assert part in str(exc)
        add(name,'Rejected with actionable message','Rejected with actionable message')
    else: raise AssertionError(name)

r=screen_reference(S,{'other':L})
add('Unconfirmed exact match is retained',{'hits':r.total_hits,'score':r.specificity_score},{'hits':1,'score':50.0})
r=screen_reference(S,{'target':L},intended_locus=TargetLocus('target',1,'+'))
add('Confirmed intended coordinate is excluded',{'hits':r.total_hits,'score':r.specificity_score},{'hits':0,'score':100.0})
r=screen_reference(S,{'target':L,'copy':L},intended_locus=TargetLocus('target',1,'+'))
add('Second exact locus remains visible',r.total_hits,1)
panel={str(i):L for i in range(301)}
scores=[screen_reference(S,panel,intended_locus=TargetLocus('0',1,'+'),max_hits=cap).specificity_score for cap in (1,250,1000)]
add('Display limits 1 250 and 1000 preserve score',scores,[0.33,0.33,0.33])
rejected('Duplicate FASTA names',lambda:parse_multifasta('>same\n'+L+'\n>same\n'+L),'Duplicate FASTA')
r=screen_reference(S,{'unknown':'N'*40})
add('All ambiguous reference has no score',r.specificity_score,None)
r=screen_reference(S,{'reverse':'AAA'+reverse_complement(L)},intended_locus=TargetLocus('reverse',7,'-'))
add('Reverse strand coordinate 7 is verified',r.intended_site_excluded,True)
near='T'+S[1:]
r=screen_reference(S,{'near':near+'AGG'},1)
add('One substitution at spacer position 1',list(r.hits[0].mismatch_positions),[1])
rec=SeqRecord(Seq(L),id='NM_DEMO.1');rec.annotations['molecule_type']='mRNA'
rejected('Transcript without exon boundaries',lambda:segments_from_ncbi_record(rec),'exon boundaries')
rec.features=[SeqFeature(FeatureLocation(0,10),type='exon'),SeqFeature(FeatureLocation(10,23),type='exon')]
rejected('Two short exons cannot create a junction target',lambda:segments_from_ncbi_record(rec),'short exons')
rec=SeqRecord(Seq('A'*30+L+'A'*30),id='NC_DEMO.1');rec.annotations['molecule_type']='DNA'
rec.features=[SeqFeature(CompoundLocation([FeatureLocation(0,30),FeatureLocation(53,83)]),type='CDS')]
segments,_=segments_from_ncbi_record(rec)
add('Compound CDS excludes intervening intron',[len(s) for _,s in segments],[30,30])
add('Non G spacer has unmodified option',next(x.sequence for x in guide_format_variants(S) if x.label=='X20'),S)
state={};reset_panel_if_changed(state,'first',3,None,250);state['oc_screens']['guide']='old';reset_panel_if_changed(state,'second',3,None,250)
add('Changed panel invalidates prior result',state['oc_screens'],{})
add('Resolved plus strand site uses coordinates 1 to 20',[(g.start,g.end,g.strand) for g in scan_opencrispr_sites(L)],[(1,20,'+')])
output={'kind':'Synthetic software correctness examples, not wet-lab validation','app_version':'1.4.0','examples':rows}
path=Path(__file__).resolve().parent/'results'/'audit_examples.json'
path.write_text(json.dumps(output,indent=2)+'\n')
print(f'{len(rows)} / {len(rows)} examples passed. Results: {path.relative_to(path.parents[2])}')
