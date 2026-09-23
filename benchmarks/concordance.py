# ============================================================================
# CONCORDANCE
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Compares benchmark outputs so scientific behaviour can be checked for consistency.
#
# HOW TO READ THIS FILE:
# 1. Read the imports/constants first to see which tools and settings are used.
# 2. Read one top-level function or class at a time.
# 3. Follow the workflow from input sequence -> candidate guides -> validation -> output.
# 4. Scientific formulas, thresholds, validation decisions and public function names
#    are intentionally preserved while readability comments are added.
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
