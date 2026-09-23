import pytest
from opencrispr_designer import scan_opencrispr_sites, design_from_segments, clean_dna

def test_no_valid_pam_returns_empty():
    assert scan_opencrispr_sites("ATATATATATATATATATATATATATAT", "G", "e") == []

@pytest.mark.parametrize("code", ["N","R","Y"])
def test_iupac_in_candidate_window_is_skipped(code):
    raw="A"*10 + code + "A"*9 + "AGG"
    assert not any(g.start==1 for g in scan_opencrispr_sites(raw,"G","e"))
    assert "N" in clean_dna(raw)

def test_single_gene_input_is_supported():
    spacer="ACGTACGTACGTACGTACGT"
    gs=design_from_segments("ONE",[("ex1","AAAA"+spacer+"AGG"+"TTTT")],100)
    assert any(g.spacer==spacer for g in gs)

def test_unrelated_sequence_has_no_false_exact_candidate():
    spacer="ACGTACGTACGTACGTACGT"
    gs=design_from_segments("OTHER",[("ex1","T"*80)],100)
    assert not any(g.spacer==spacer for g in gs)

def test_exon_junction_avoidance_explicitly():
    gs=design_from_segments("G",[("ex1","A"*20),("ex2","AGG"+"C"*30)],100)
    assert not any(g.spacer=="A"*20 and g.pam=="AGG" for g in gs)
