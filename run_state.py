# ============================================================================
# RUN STATE
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Keeps workflow state and run metadata organised while the interactive application is being used.
#
# HOW TO READ THIS FILE:
# 1. Read the imports/constants first to see which tools and settings are used.
# 2. Read one top-level function or class at a time.
# 3. Follow the workflow from input sequence -> candidate guides -> validation -> output.
# 4. Scientific formulas, thresholds, validation decisions and public function names
#    are intentionally preserved while readability comments are added.
#
# MAIN TOP-LEVEL PARTS:
# - function: create_run_snapshot
# - function: panel_result_key
# - function: export_input_fasta
# - function: export_run
# - function: cas_offinder_input
# ============================================================================

"""Immutable-by-copy run metadata and external specificity handoff."""
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import re

APP_VERSION = "1.4.0"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: create_run_snapshot
# ----------------------------------------------------------------------------
def create_run_snapshot(records, settings):
    inputs = {gene: {**rec.provenance_dict(), "segments": [list(s) for s in rec.segments]}
              for gene, rec in sorted(records.items())}
    payload = {"version": APP_VERSION, "settings": dict(settings), "inputs": inputs}
    # Stable rerun identity excludes retrieval time but includes names, sequences and settings.
    identity = {"version": APP_VERSION, "settings": settings,
                "segments": {g: inputs[g]["segments"] for g in inputs}}
    payload["run_id"] = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    payload["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    return payload



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: panel_result_key
# ----------------------------------------------------------------------------
def panel_result_key(run_id, spacer, raw_fasta, mismatch_radius):
    digest = hashlib.sha256(raw_fasta.encode()).hexdigest()
    return f"plant_panel::{run_id}::{spacer}::{mismatch_radius}::{digest}"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: export_input_fasta
# ----------------------------------------------------------------------------
def export_input_fasta(records):
    return "\n".join(f">{gene}|segment_{i}\n{seq}" for gene, rec in records.items()
                     for i, (_, seq) in enumerate(rec.segments, 1)) + "\n"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: export_run
# ----------------------------------------------------------------------------
def export_run(snapshot, guides, reports, panel_reports=None, fallback=None):
    return {"app": "Plant MultiGene gRNA Designer", **snapshot,
            "guides": [{"guide": asdict(g), "validation": asdict(r)} for g, r in zip(guides, reports)],
            "panel_screens": {s: asdict(r) for s, r in (panel_reports or {}).items()},
            "fallback": None if fallback is None else {
                **fallback, "guides": [asdict(g) for g in fallback["guides"]]}}



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: cas_offinder_input
# ----------------------------------------------------------------------------
def cas_offinder_input(spacers, genome_path, max_mismatches=3, include_nag=False):
    """Cas-OFFinder 2 no-bulge input; export only, never executes an external scan."""
    if not genome_path.strip() or "\n" in genome_path or "\r" in genome_path:
        raise ValueError("Provide a nonempty single-line local genome FASTA path.")
    if not isinstance(max_mismatches, int) or not 0 <= max_mismatches <= 20:
        raise ValueError("Mismatch radius must be an integer from 0 to 20.")
    spacers = list(dict.fromkeys(spacers))
    if not spacers or any(not re.fullmatch("[ACGT]{20}", s) for s in spacers):
        raise ValueError("Export requires at least one resolved 20-nt spacer.")
    pattern = "N" * 20 + ("NRG" if include_nag else "NGG")
    return "\n".join([genome_path.strip(), pattern] +
                     [f"{s}NNN {max_mismatches}" for s in spacers]) + "\n"
