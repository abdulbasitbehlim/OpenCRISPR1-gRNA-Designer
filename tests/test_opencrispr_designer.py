from opencrispr_designer import (
    clean_dna,
    reverse_complement,
    scan_opencrispr_sites,
    design_from_segments,
    guide_format_variants,
    mismatch_positions,
    mit_offtarget_pair_score,
    mit_specificity,
    screen_local_reference,
)


def test_clean_and_reverse_complement():
    assert clean_dna(">x\nacgu RY") == "ACGTNN"
    assert reverse_complement("ACGTN") == "NACGT"


def test_forward_ngg_discovery():
    spacer = "ACGTACGTACGTACGTACGT"
    seq = "AAAA" + spacer + "AGG" + "TTTT"
    guides = scan_opencrispr_sites(seq, gene="G1", segment_id="ex1")
    assert any(g.spacer == spacer and g.pam == "AGG" and g.strand == "+" for g in guides)


def test_reverse_ngg_discovery():
    desired = "GACTGACTGACTGACTGACT"
    downstream = reverse_complement(desired)
    seq = "AAAA" + "CCA" + downstream + "TTTT"
    guides = scan_opencrispr_sites(seq)
    assert any(g.spacer == desired and g.pam == "TGG" and g.strand == "-" for g in guides)


def test_segments_are_scanned_independently_no_junction_candidate():
    seg1 = "A" * 20
    seg2 = "AGG" + "C" * 20
    guides = design_from_segments("G", [("a", seg1), ("b", seg2)], max_guides=100)
    assert not any(g.spacer == "A" * 20 and g.pam == "AGG" for g in guides)


def test_guide_format_variants_natural_g():
    spacer = "G" + "A" * 19
    formats = guide_format_variants(spacer)
    assert formats[0].label == "GX19"
    assert formats[0].sequence == spacer
    assert any(x.label == "gX20" and len(x.sequence) == 21 for x in formats)


def test_guide_format_variants_non_g():
    spacer = "A" * 20
    formats = guide_format_variants(spacer)
    assert formats[0].label == "gX19"
    assert formats[0].sequence == "G" + "A" * 19
    assert formats[1].label == "gX20"
    assert formats[1].sequence == "G" + spacer


def test_mit_exact_and_mismatch():
    a = "ACGTACGTACGTACGTACGT"
    b = a[:-1] + ("A" if a[-1] != "A" else "C")
    assert mismatch_positions(a, a) == ()
    assert mismatch_positions(a, b) == (20,)
    assert mit_offtarget_pair_score(a, a) == 1.0
    assert 0 <= mit_offtarget_pair_score(a, b) < 1.0
    assert mit_specificity([]) == 100.0


def test_local_reference_excludes_one_intended_exact_and_keeps_second_exact():
    spacer = "ACGTACGTACGTACGTACGT"
    locus = spacer + "AGG"
    panel = {"intended": locus, "duplicate": locus}
    spec, hits = screen_local_reference(spacer, panel, max_mismatches=0, exclude_one_exact=True)
    assert len(hits) == 1
    assert hits[0].mismatches == 0
    assert spec == 50.0


def test_local_reference_retains_near_match():
    guide = "ACGTACGTACGTACGTACGT"
    near = guide[:-1] + ("A" if guide[-1] != "A" else "C")
    panel = {"intended": guide + "AGG", "near": near + "TGG"}
    spec, hits = screen_local_reference(guide, panel, max_mismatches=1)
    assert any(h.mismatches == 1 for h in hits)
    assert 0 < spec < 100


def test_ambiguity_behavior_is_explicit():
    from opencrispr_designer import ambiguity_summary, scan_opencrispr_sites
    spacer = "ACGTACGTACGTACGTACGA"
    a = ambiguity_summary(">x\nACGTRYN")
    assert a["count"] == 3
    assert set(a["codes"]) == {"N", "R", "Y"}
    sites = scan_opencrispr_sites("AAAAA" + spacer + "NGG" + "AAAA", "G", "seg")
    assert any(g.spacer == spacer for g in sites)
    bad = spacer[:4] + "N" + spacer[5:]
    sites2 = scan_opencrispr_sites("AAAAA" + bad + "AGG" + "AAAA", "G", "seg")
    assert not any(g.start == 6 for g in sites2)


def test_specificity_method_is_exported_as_legacy_baseline_label():
    from opencrispr_designer import GuideRNA, guide_row
    g = GuideRNA("G", "seg", "ACGT"*5, "AGG", "+", 1, 20, 50.0, 80.0)
    g.specificity_score = 77.0
    g.specificity_method = "MIT/Hsu 2013 legacy baseline (supplied FASTA panel)"
    row = guide_row(g)
    assert row["Specificity method"].startswith("MIT/Hsu 2013 legacy baseline")
