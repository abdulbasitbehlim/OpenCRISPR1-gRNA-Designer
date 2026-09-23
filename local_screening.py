# ============================================================================
# LOCAL SCREENING
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Performs local-reference screening used to review possible guide specificity concerns.
#
# HOW TO READ THIS FILE:
# 1. Start with imports and constants to see the dependencies and fixed settings.
# 2. Read one function/class at a time rather than the whole file at once.
# 3. Follow the workflow from sequence input -> candidate discovery -> screening -> validation.
# 4. Scientific calculations, thresholds, validation rules and public APIs are
#    intentionally preserved while readability explanations are added.
#
# MAIN TOP-LEVEL PARTS:
# - class: TargetLocus
# - class: PanelSite
# - class: PanelScreen
# - function: normalize_panel
# - function: panel_fingerprint
# - function: iter_panel_sites
# - function: exact_target_loci
# - function: _resolved_guide
# - function: screen_reference
# ============================================================================

"""Bounded, coordinate-aware NGG panel screening, not genome-wide prediction.

Only the display is capped. Counts and the legacy MIT summary include every
resolved NGG hit within the mismatch radius. No site is presumed on-target.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import heapq
import json
from typing import Mapping

from opencrispr_designer import (
    OffTargetHit, clean_dna, mismatch_positions, mit_offtarget_pair_score,
    mit_specificity, reverse_complement,
)

MAX_PANEL_BP = 2_000_000
MAX_DISPLAY_HITS = 5_000


@dataclass(frozen=True)
class TargetLocus:
    contig: str
    start: int  # 1-based leftmost spacer base on supplied sequence
    strand: str


@dataclass(frozen=True)
class PanelSite:
    contig: str
    start: int
    strand: str
    spacer: str
    pam: str

    @property
    def locus(self) -> TargetLocus:
        return TargetLocus(self.contig, self.start, self.strand)


@dataclass(frozen=True)
class PanelScreen:
    guide: str
    panel_sha256: str
    panel_bp: int
    contig_count: int
    ambiguous_bases: int
    scanned_sites: int
    max_mismatches: int
    display_limit: int
    intended_locus: TargetLocus | None
    intended_site_excluded: bool
    total_hits: int
    mismatch_counts: tuple[int, ...]
    specificity_score: float | None
    hits: tuple[OffTargetHit, ...]
    warnings: tuple[str, ...]
    screened_at_utc: str
    method: str = "MIT/Hsu legacy baseline; supplied FASTA; NGG; substitutions only"

    @property
    def truncated(self) -> bool:
        return self.total_hits > len(self.hits)

    def to_dict(self) -> dict:
        return {**asdict(self), "display_truncated": self.truncated,
                "coordinate_system": "1-based inclusive; leftmost spacer base on supplied contig"}



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: normalize_panel
# ----------------------------------------------------------------------------
def normalize_panel(panel: Mapping[str, str]) -> dict[str, str]:
    if not panel:
        raise ValueError("Supply a non-empty reference FASTA panel.")
    if sum(len(s) for s in panel.values()) > MAX_PANEL_BP:
        raise ValueError(f"Local panels are limited to {MAX_PANEL_BP:,} bases; use an indexed genome workflow for larger references.")
    out = {}
    for name, raw in panel.items():
        if not name or not name.strip():
            raise ValueError("Reference records need non-empty identifiers.")
        seq = clean_dna(raw)
        if not seq:
            raise ValueError(f"Reference record {name!r} is empty.")
        out[name] = seq
    return out



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: panel_fingerprint
# ----------------------------------------------------------------------------
def panel_fingerprint(panel: Mapping[str, str]) -> str:
    return hashlib.sha256(json.dumps(sorted(panel.items()), separators=(",", ":")).encode()).hexdigest()



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: iter_panel_sites
# ----------------------------------------------------------------------------
def iter_panel_sites(panel: Mapping[str, str]):
    """Yield sites without calculating sequence-quality scores or sorting a genome."""
    for name, seq in panel.items():
        for i in range(len(seq) - 22):
            window = seq[i:i + 23]
            if window[21:23] == "GG" and "N" not in window[:20]:
                yield PanelSite(name, i + 1, "+", window[:20], window[20:])
            if window[:2] == "CC" and "N" not in window[3:]:
                yield PanelSite(name, i + 4, "-", reverse_complement(window[3:]), reverse_complement(window[:3]))



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: exact_target_loci
# ----------------------------------------------------------------------------
def exact_target_loci(guide: str, panel: Mapping[str, str]) -> list[TargetLocus]:
    g = _resolved_guide(guide)
    return [s.locus for s in iter_panel_sites(normalize_panel(panel)) if s.spacer == g]



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _resolved_guide
# ----------------------------------------------------------------------------
def _resolved_guide(guide: str) -> str:
    g = clean_dna(guide)
    if len(g) != 20 or "N" in g:
        raise ValueError("Screening requires a resolved 20-base ACGT targeting spacer.")
    return g



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: screen_reference
# ----------------------------------------------------------------------------
def screen_reference(
    guide: str, panel: Mapping[str, str], max_mismatches: int = 3,
    intended_locus: TargetLocus | None = None, max_hits: int = 250,
) -> PanelScreen:
    g = _resolved_guide(guide)
    if not isinstance(max_mismatches, int) or not 0 <= max_mismatches <= 4:
        raise ValueError("Maximum mismatches must be an integer from 0 to 4.")
    if not isinstance(max_hits, int) or not 1 <= max_hits <= MAX_DISPLAY_HITS:
        raise ValueError(f"Display limit must be from 1 to {MAX_DISPLAY_HITS}.")
    panel = normalize_panel(panel)
    if intended_locus is not None and (
        intended_locus.contig not in panel or intended_locus.start < 1
        or intended_locus.strand not in {"+", "-"}
    ):
        raise ValueError("Intended locus must identify a supplied contig, positive start and strand.")

    counts = [0] * (max_mismatches + 1)
    scanned = 0
    excluded = False
    risk_sum = 0.0
    retained = []
    for site in iter_panel_sites(panel):
        scanned += 1
        pos = mismatch_positions(g, site.spacer)
        if site.locus == intended_locus:
            if pos:
                raise ValueError("The specified intended locus does not exactly match the targeting spacer.")
            excluded = True
            continue
        if len(pos) > max_mismatches:
            continue
        seed_mm = sum(p >= 13 for p in pos)
        risk = mit_offtarget_pair_score(g, site.spacer)
        hit = OffTargetHit(site.contig, site.spacer, site.pam, site.strand, site.start,
                           len(pos), pos, seed_mm, risk)
        counts[len(pos)] += 1
        risk_sum += risk
        # The heap root is the worst retained hit. Numeric serial breaks ties.
        quality = (-len(pos), -seed_mm, risk, -scanned)
        item = (quality, hit)
        if len(retained) < max_hits:
            heapq.heappush(retained, item)
        elif quality > retained[0][0]:
            heapq.heapreplace(retained, item)
    if intended_locus is not None and not excluded:
        raise ValueError("The intended locus was not found as a resolved NGG site; check coordinates, strand and reference.")
    ambiguity = sum(seq.count("N") for seq in panel.values())
    warnings = []
    if not excluded:
        warnings.append("Intended locus unconfirmed: no exact match was automatically excluded.")
    if ambiguity:
        warnings.append("Reference contains ambiguous bases; unresolved spacer windows are skipped and the screen is incomplete.")
    if not scanned:
        warnings.append("No resolved NGG spacer windows were searchable; no specificity score is available.")
    if sum(counts) > max_hits:
        warnings.append("Hit display is truncated; score and counts include all matching sites within the selected radius.")
    warnings.append("Scope: supplied sequences, NGG PAMs and substitutions only; no bulges, alternate PAMs, chromatin or allele model.")
    hits = sorted((h for _, h in retained), key=lambda h: (h.mismatches, h.seed_mismatches, -h.mit_pair_risk, h.contig, h.start, h.strand))
    return PanelScreen(g, panel_fingerprint(panel), sum(map(len, panel.values())), len(panel),
                       ambiguity, scanned, max_mismatches, max_hits, intended_locus, excluded,
                       sum(counts), tuple(counts), mit_specificity([risk_sum]) if scanned else None,
                       tuple(hits), tuple(warnings), datetime.now(timezone.utc).isoformat())
