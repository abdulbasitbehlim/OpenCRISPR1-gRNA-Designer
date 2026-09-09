#!/usr/bin/env python3
"""OpenCRISPR-1 guide spacer discovery and transparent ranking.

Scientific scope:
- OpenCRISPR-1 is treated as an NGG-preferring Cas9-like nuclease.
- Candidate spacers are 20 nt and are discovered on both strands.
- The built-in sequence score is an explainable prioritization heuristic, not an
  OpenCRISPR-specific learned activity model.
- Local FASTA specificity uses the published MIT/Hsu mismatch framework as a
  legacy baseline for transparent panel screening. It is not state-of-the-art
  off-target prediction and a supplied panel is not a whole-genome analysis.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import re

DNA = frozenset("ACGTN")
IUPAC_AMBIGUOUS = frozenset("NRYWSKMBDHVX")
MIT_WEIGHTS = (
    0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583,
)


@dataclass
class GuideRNA:
    gene: str
    segment_id: str
    spacer: str
    pam: str
    strand: str
    start: int
    end: int
    gc_percent: float
    sequence_score: float
    opencrispr_compatible: bool = True
    specificity_score: Optional[float] = None
    specificity_method: Optional[str] = None
    off_target_count: Optional[int] = None
    notes: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class OffTargetHit:
    contig: str
    spacer: str
    pam: str
    strand: str
    start: int
    mismatches: int
    mismatch_positions: Tuple[int, ...]
    seed_mismatches: int
    mit_pair_risk: float


@dataclass(frozen=True)
class GuideFormat:
    label: str
    sequence: str
    length: int
    description: str


def clean_dna(raw: str) -> str:
    if not raw:
        return ""
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip() and not ln.lstrip().startswith(">")]
    seq = "".join(lines).upper().replace("U", "T")
    seq = re.sub(r"[\s\d]", "", seq)
    seq = re.sub(r"[RYSWKMBDHVX]", "N", seq)
    bad = sorted(set(seq) - DNA)
    if bad:
        raise ValueError("Unsupported sequence character(s): " + ", ".join(bad))
    return seq


def ambiguity_summary(raw: str) -> Dict[str, object]:
    """Report IUPAC ambiguity and define how scanning handles it.

    Ambiguous symbols are normalized to N. Candidate spacers containing N are
    skipped rather than scored as mismatches. At the PAM, ambiguity is accepted
    only at the degenerate N position of an otherwise resolved NGG/CCN motif.
    """
    body = "".join(
        ln.strip() for ln in (raw or "").splitlines()
        if ln.strip() and not ln.lstrip().startswith(">")
    ).upper().replace("U", "T")
    body = re.sub(r"[\s\d]", "", body)
    codes = sorted({c for c in body if c in IUPAC_AMBIGUOUS})
    positions = tuple(i + 1 for i, c in enumerate(body) if c in IUPAC_AMBIGUOUS)
    return {"count": len(positions), "codes": codes, "positions": positions}


def reverse_complement(seq: str) -> str:
    return clean_dna(seq).translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def gc_percent(seq: str) -> float:
    s = clean_dna(seq)
    return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def sequence_quality_score(spacer: str) -> float:
    """Transparent 0-100 ranking, deliberately not called OpenCRISPR efficiency."""
    s = clean_dna(spacer)
    if len(s) != 20 or "N" in s:
        return 0.0
    gc = gc_percent(s)
    score = 58.0
    if 40 <= gc <= 60:
        score += 22
    elif 30 <= gc <= 70:
        score += 10
    else:
        score -= 15
    if "TTTT" in s:
        score -= 20
    if re.search(r"A{5,}|C{5,}|G{5,}|T{5,}", s):
        score -= 12
    if s[-1] == "G":
        score += 5
    if s[-2] in "AG":
        score += 3
    if s[0] == "G":
        score += 2
    return round(max(0.0, min(100.0, score)), 1)


def scan_opencrispr_sites(sequence: str, gene: str = "target", segment_id: str = "segment") -> List[GuideRNA]:
    """Discover OpenCRISPR-1-compatible 20-nt spacers with NGG PAMs on both strands."""
    seq = clean_dna(sequence)
    guides: List[GuideRNA] = []
    for pam_start in range(20, len(seq) - 2):
        pam = seq[pam_start:pam_start + 3]
        if len(pam) == 3 and pam[1:] == "GG":
            spacer = seq[pam_start - 20:pam_start]
            if "N" not in spacer:
                guides.append(GuideRNA(
                    gene=gene, segment_id=segment_id, spacer=spacer, pam=pam,
                    strand="+", start=pam_start - 20 + 1, end=pam_start,
                    gc_percent=round(gc_percent(spacer), 1),
                    sequence_score=sequence_quality_score(spacer),
                    notes=["20-nt spacer with NGG PAM; compatible with the published OpenCRISPR-1 target preference."],
                ))
    for pam_start in range(0, len(seq) - 22):
        raw_pam = seq[pam_start:pam_start + 3]
        if len(raw_pam) == 3 and raw_pam[:2] == "CC":
            downstream = seq[pam_start + 3:pam_start + 23]
            if len(downstream) == 20 and "N" not in downstream:
                spacer = reverse_complement(downstream)
                pam = reverse_complement(raw_pam)
                guides.append(GuideRNA(
                    gene=gene, segment_id=segment_id, spacer=spacer, pam=pam,
                    strand="-", start=pam_start + 4, end=pam_start + 23,
                    gc_percent=round(gc_percent(spacer), 1),
                    sequence_score=sequence_quality_score(spacer),
                    notes=["20-nt spacer with reverse-strand NGG PAM; reported in synthesis orientation."],
                ))
    guides.sort(key=lambda g: (g.sequence_score, -abs(g.gc_percent - 50.0)), reverse=True)
    return guides


def design_from_segments(gene: str, segments: Sequence[Tuple[str, str]], max_guides: int = 50, min_score: float = 0.0) -> List[GuideRNA]:
    guides: List[GuideRNA] = []
    for segment_id, seq in segments:
        guides.extend(scan_opencrispr_sites(seq, gene=gene, segment_id=segment_id))
    unique: Dict[Tuple[str, str, int, str], GuideRNA] = {}
    for g in guides:
        if g.sequence_score >= min_score:
            unique[(g.spacer, g.segment_id, g.start, g.strand)] = g
    out = list(unique.values())
    out.sort(key=lambda g: (g.sequence_score, -abs(g.gc_percent - 50.0), g.spacer), reverse=True)
    return out[:max_guides]


def guide_format_variants(spacer: str) -> List[GuideFormat]:
    """Return common 5'-G guide-expression formats explicitly evaluated in a 2026 OpenCRISPR study."""
    s = clean_dna(spacer)
    if len(s) != 20:
        raise ValueError("Guide spacer must be 20 nt.")
    variants: List[GuideFormat] = []
    if s.startswith("G"):
        variants.append(GuideFormat(
            "GX19", s, 20,
            "Natural 5'-G 20-nt spacer; no targeting base is changed.",
        ))
    else:
        variants.append(GuideFormat(
            "gX19", "G" + s[1:], 20,
            "5' base substituted with G for promoter compatibility; this creates one spacer-target mismatch at position 1 and must be evaluated carefully.",
        ))
    variants.append(GuideFormat(
        "gX20", "G" + s, 21,
        "An additional 5'-G is appended while retaining the full 20-nt targeting spacer.",
    ))
    return variants


def mismatch_positions(a: str, b: str) -> Tuple[int, ...]:
    a, b = clean_dna(a), clean_dna(b)
    if len(a) != 20 or len(b) != 20:
        raise ValueError("MIT comparison requires two 20-nt spacers.")
    return tuple(i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y)


def mit_offtarget_pair_score(guide: str, target: str) -> float:
    pos = mismatch_positions(guide, target)
    if not pos:
        return 1.0
    m = len(pos)
    d = 19 if m == 1 else (max(pos) - min(pos)) / (m - 1)
    term = 1.0
    for p in pos:
        term *= 1.0 - MIT_WEIGHTS[p - 1]
    return max(0.0, min(1.0, term * (1 / (m * m)) * (1 / ((((19 - d) / 19) * 4) + 1))))


def mit_specificity(pair_scores: Iterable[float]) -> float:
    return round(100.0 / (1.0 + sum(pair_scores)), 2)


def _scan_panel_sites(panel: Mapping[str, str]) -> List[GuideRNA]:
    sites: List[GuideRNA] = []
    for name, seq in panel.items():
        sites.extend(scan_opencrispr_sites(seq, gene=name, segment_id=name))
    return sites


def screen_local_reference(
    guide: str,
    panel: Mapping[str, str],
    max_mismatches: int = 3,
    exclude_one_exact: bool = True,
    max_hits: int = 250,
) -> Tuple[float, List[OffTargetHit]]:
    """PAM-aware local FASTA screen using MIT/Hsu as a legacy baseline.

    MIT/Hsu (2013) is retained because it is transparent and dependency-free,
    but later empirical metrics such as CFD and newer ML predictors can better
    model off-target activity in many settings. Treat this result as a supplied-
    panel baseline, not as the best available genome-wide specificity estimate.
    """
    g = clean_dna(guide)
    if len(g) != 20:
        raise ValueError("Specificity screening requires a 20-nt spacer.")
    hits: List[OffTargetHit] = []
    exact_excluded = False
    for site in _scan_panel_sites(panel):
        pos = mismatch_positions(g, site.spacer)
        if len(pos) > max_mismatches:
            continue
        if exclude_one_exact and not pos and not exact_excluded:
            exact_excluded = True
            continue
        seed_mm = sum(p >= 13 for p in pos)
        hits.append(OffTargetHit(
            contig=site.gene, spacer=site.spacer, pam=site.pam,
            strand=site.strand, start=site.start, mismatches=len(pos),
            mismatch_positions=pos, seed_mismatches=seed_mm,
            mit_pair_risk=mit_offtarget_pair_score(g, site.spacer),
        ))
    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches, -h.mit_pair_risk, h.contig, h.start))
    hits = hits[:max_hits]
    return mit_specificity(h.mit_pair_risk for h in hits), hits


def guide_row(g: GuideRNA) -> Dict[str, object]:
    return {
        "Spacer (20 nt)": g.spacer,
        "PAM": g.pam,
        "Strand": g.strand,
        "Segment": g.segment_id,
        "Start": g.start,
        "End": g.end,
        "GC%": g.gc_percent,
        "Sequence quality": g.sequence_score,
        "OpenCRISPR-1 compatible": g.opencrispr_compatible,
        "MIT specificity": g.specificity_score,
        "Specificity method": g.specificity_method or "Not screened",
        "Off-target hits": g.off_target_count,
    }
