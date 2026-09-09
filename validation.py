#!/usr/bin/env python3
"""Validation rules for OpenCRISPR-1 gRNA Designer.

The validator checks structural OpenCRISPR-1 compatibility, transparent sequence
quality flags, guide-expression formatting, and optional local-reference
specificity. MIT/Hsu is explicitly treated as a legacy supplied-panel baseline,
not as state-of-the-art off-target prediction. The validator does not claim to
predict experimental editing efficacy.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence
import re

from opencrispr_designer import GuideRNA, OffTargetHit, guide_format_variants


@dataclass(frozen=True)
class ValidationCheck:
    check: str
    status: str
    value: str
    explanation: str


@dataclass(frozen=True)
class ValidationReport:
    status: str
    checks: Sequence[ValidationCheck]
    pass_count: int
    review_count: int
    fail_count: int
    specificity_status: str

    @property
    def issues(self) -> str:
        msgs = [c.check for c in self.checks if c.status in {"REVIEW", "FAIL"}]
        return "; ".join(msgs) if msgs else "None"


def _add(checks: List[ValidationCheck], name: str, status: str, value: object, explanation: str) -> None:
    checks.append(ValidationCheck(name, status, str(value), explanation))


def validate_opencrispr_guide(
    guide: GuideRNA,
    min_sequence_score: float = 35.0,
    local_hits: Optional[Iterable[OffTargetHit]] = None,
    specificity_score: Optional[float] = None,
    min_specificity_review: float = 50.0,
) -> ValidationReport:
    checks: List[ValidationCheck] = []
    spacer = (guide.spacer or "").upper()
    pam = (guide.pam or "").upper()

    _add(
        checks,
        "Spacer length",
        "PASS" if len(spacer) == 20 else "FAIL",
        f"{len(spacer)} nt",
        "This OpenCRISPR-1 workflow designs 20-nt targeting spacers.",
    )
    canonical = bool(re.fullmatch(r"[ACGT]{20}", spacer))
    _add(
        checks,
        "Canonical DNA spacer",
        "PASS" if canonical else "FAIL",
        "ACGT only" if canonical else "Ambiguous/invalid bases detected",
        "The targeting spacer should contain only resolved DNA bases.",
    )
    pam_ok = len(pam) == 3 and pam[1:] == "GG" and pam[0] in "ACGTN"
    _add(
        checks,
        "OpenCRISPR-1 PAM compatibility",
        "PASS" if pam_ok else "FAIL",
        pam or "Missing",
        "OpenCRISPR-1 shows its strongest published compatibility at NGG PAM targets; this tool therefore treats NGG as the supported design mode. The first PAM base may be represented as N because it is the degenerate position, but the two G bases must be resolved.",
    )
    _add(
        checks,
        "OpenCRISPR compatibility flag",
        "PASS" if guide.opencrispr_compatible else "FAIL",
        str(bool(guide.opencrispr_compatible)),
        "The guide object must originate from or satisfy the supported OpenCRISPR-1 design path.",
    )

    gc = guide.gc_percent
    _add(
        checks,
        "GC content",
        "PASS" if 40.0 <= gc <= 60.0 else "REVIEW",
        f"{gc:.1f}%",
        "40–60% is used as a preferred review band, not as an absolute experimental cutoff.",
    )

    _add(
        checks,
        "Sequence-quality threshold",
        "PASS" if guide.sequence_score >= min_sequence_score else "REVIEW",
        f"{guide.sequence_score:.1f} / threshold {min_sequence_score:.1f}",
        "This is the app's transparent heuristic ranking threshold, not a learned OpenCRISPR efficacy probability.",
    )

    poly_t = "TTTT" in spacer
    _add(
        checks,
        "Poly-T motif",
        "REVIEW" if poly_t else "PASS",
        "TTTT present" if poly_t else "No TTTT",
        "A TTTT motif can interfere with common Pol III/U6 guide-expression systems and should be reviewed in the chosen construct.",
    )
    long_homopolymer = bool(re.search(r"A{5,}|C{5,}|G{5,}|T{5,}", spacer))
    _add(
        checks,
        "Long homopolymer",
        "REVIEW" if long_homopolymer else "PASS",
        "≥5-base run present" if long_homopolymer else "No ≥5-base run",
        "Long homopolymers are a synthesis/expression review flag, not an absolute failure rule.",
    )

    try:
        variants = guide_format_variants(spacer)
        preserves_full = any(v.label in {"GX19", "gX20"} for v in variants)
        if spacer.startswith("G"):
            fmt_value = "GX19 available (native 5′ G)"
            fmt_note = "The spacer already begins with G, so a native 20-nt GX19 expression format is available."
        else:
            fmt_value = "gX20 preserves full spacer; gX19 changes position 1"
            fmt_note = "For a non-G-starting spacer, gX20 preserves all 20 targeting bases; gX19 introduces a position-1 substitution and is therefore flagged for careful review."
        _add(checks, "5′-G expression format", "PASS" if preserves_full else "FAIL", fmt_value, fmt_note)
    except Exception as exc:
        _add(checks, "5′-G expression format", "FAIL", "Unavailable", str(exc))

    if local_hits is None and specificity_score is None:
        specificity_status = "NOT SCREENED"
        _add(
            checks,
            "Local specificity",
            "INFO",
            "Not screened",
            "No supplied-FASTA specificity screen has been run. This does not imply genome-wide specificity.",
        )
    else:
        hits = list(local_hits or [])
        extra_exact = sum(h.mismatches == 0 for h in hits)
        one_mm = sum(h.mismatches == 1 for h in hits)
        low_spec = specificity_score is not None and specificity_score < min_specificity_review
        if extra_exact:
            specificity_status = "FAIL (SUPPLIED PANEL)"
            _add(
                checks,
                "Local specificity",
                "FAIL",
                f"{extra_exact} additional exact hit(s); MIT={specificity_score}",
                "After excluding one presumed intended site, an additional exact PAM-compatible hit remains in the supplied panel and requires resolution before prioritization.",
            )
        elif one_mm or low_spec:
            specificity_status = "REVIEW (SUPPLIED PANEL)"
            reason = f"{one_mm} one-mismatch hit(s); MIT={specificity_score}"
            _add(
                checks,
                "Local specificity",
                "REVIEW",
                reason,
                "Close near matches or a low user-defined MIT/Hsu review score require manual review. MIT/Hsu (2013) is a legacy baseline; the supplied panel is not a whole-genome screen and newer empirical/ML methods may perform better.",
            )
        else:
            specificity_status = "PASS (SUPPLIED PANEL)"
            _add(
                checks,
                "Local specificity",
                "PASS",
                f"{len(hits)} retained hit(s); MIT={specificity_score}",
                "No extra exact or one-mismatch hit was retained and the legacy MIT/Hsu score met the review threshold in the supplied panel. This does not establish genome-wide specificity.",
            )

    fail_count = sum(c.status == "FAIL" for c in checks)
    review_count = sum(c.status == "REVIEW" for c in checks)
    pass_count = sum(c.status == "PASS" for c in checks)
    status = "FAIL" if fail_count else "REVIEW" if review_count else "PASS"
    return ValidationReport(status, tuple(checks), pass_count, review_count, fail_count, specificity_status)


def validation_summary_row(report: ValidationReport) -> dict:
    return {
        "Validation": report.status,
        "Validation PASS": report.pass_count,
        "Validation REVIEW": report.review_count,
        "Validation FAIL": report.fail_count,
        "Specificity validation": report.specificity_status,
        "Validation issues": report.issues,
    }
