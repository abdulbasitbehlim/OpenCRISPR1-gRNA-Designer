# ============================================================================
# CONCORDANCE
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Compares benchmark outputs so scientific behaviour can be checked for consistency.
#
# HOW TO READ THIS FILE:
# 1. Start with imports and constants to see the dependencies and fixed settings.
# 2. Read one function/class at a time rather than the whole file at once.
# 3. Follow the workflow from sequence input -> candidate discovery -> screening -> validation.
# 4. Scientific calculations, thresholds, validation rules and public APIs are
#    intentionally preserved while readability explanations are added.
#
# MAIN TOP-LEVEL PARTS:
# - function: load_guides
# - function: compare
# ============================================================================

from __future__ import annotations
import csv
from pathlib import Path
from typing import Iterable, Set


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: load_guides
# ----------------------------------------------------------------------------
def load_guides(path: str | Path, column: str = "guide") -> Set[str]:
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows=csv.DictReader(f)
        return {r[column].strip().upper().replace('U','T') for r in rows if r.get(column)}


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: compare
# ----------------------------------------------------------------------------
def compare(ours: Iterable[str], external: Iterable[str]):
    a,b=set(ours),set(external)
    return {
        "ours_total":len(a), "external_total":len(b), "found_by_both":len(a&b),
        "ours_only":len(a-b), "external_only":len(b-a),
        "both_guides":sorted(a&b), "ours_only_guides":sorted(a-b), "external_only_guides":sorted(b-a)
    }
