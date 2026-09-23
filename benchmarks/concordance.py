from __future__ import annotations
import csv
from pathlib import Path
from typing import Iterable, Set

def load_guides(path: str | Path, column: str = "guide") -> Set[str]:
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows=csv.DictReader(f)
        return {r[column].strip().upper().replace('U','T') for r in rows if r.get(column)}

def compare(ours: Iterable[str], external: Iterable[str]):
    a,b=set(ours),set(external)
    return {
        "ours_total":len(a), "external_total":len(b), "found_by_both":len(a&b),
        "ours_only":len(a-b), "external_only":len(b-a),
        "both_guides":sorted(a&b), "ours_only_guides":sorted(a-b), "external_only_guides":sorted(b-a)
    }
