from opencrispr_designer import GuideRNA, OffTargetHit
from validation import validate_opencrispr_guide


def _guide(spacer="ACGTACGTACGTACGTACGT", pam="AGG"):
    return GuideRNA("GENE", "seg", spacer, pam, "+", 1, 20, 50.0, 80.0)


def test_opencrispr_core_validation_passes():
    r = validate_opencrispr_guide(_guide())
    assert r.status == "PASS"
    assert r.fail_count == 0
    assert r.specificity_status == "NOT SCREENED"


def test_non_ngg_pam_fails():
    r = validate_opencrispr_guide(_guide(pam="AGA"))
    assert r.status == "FAIL"
    assert any(c.check == "OpenCRISPR-1 PAM compatibility" and c.status == "FAIL" for c in r.checks)


def test_additional_exact_local_hit_fails_specificity():
    hit = OffTargetHit("other", "ACGTACGTACGTACGTACGT", "TGG", "+", 10, 0, tuple(), 0, 1.0)
    r = validate_opencrispr_guide(_guide(), local_hits=[hit], specificity_score=50.0)
    assert r.status == "FAIL"
    assert r.specificity_status.startswith("FAIL")


def test_degenerate_n_pam_position_is_structurally_accepted():
    from opencrispr_designer import GuideRNA
    from validation import validate_opencrispr_guide
    g = GuideRNA("G", "seg", "ACGT"*5, "NGG", "+", 1, 20, 50.0, 80.0)
    r = validate_opencrispr_guide(g)
    pam_check = next(c for c in r.checks if c.check == "OpenCRISPR-1 PAM compatibility")
    assert pam_check.status == "PASS"
