#!/usr/bin/env python3
"""Validation rules for Plant MultiGene gRNA Designer.

Validation deliberately separates hard design requirements from review-level
sequence-quality flags and optional local-reference screening. A PASS is not a
claim of experimental efficacy and a supplied FASTA panel is not equivalent to
whole-genome specificity analysis.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence
import re

from plant_multiguide import PanelHit, SharedGuide, PanelScreen, mismatch_positions, seed_mismatch_count, compatibility_proxy, gc_percent


@dataclass(frozen=True)
class ValidationCheck:
    check: str
    status: str  # PASS / REVIEW / FAIL / INFO
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


def validate_shared_guide(
    guide: SharedGuide,
    expected_genes: int,
    max_mismatches_per_gene: int = 2,
    max_seed_mismatches_per_gene: int = 1,
    min_compatibility: float = 55.0,
    panel_hits: Optional[Iterable[PanelHit]] = None,
    expected_gene_ids: Optional[Iterable[str]] = None,
    input_warnings: Optional[Iterable[str]] = None,
) -> ValidationReport:
    """Validate one shared guide using explicit, inspectable rules.

    Hard failures concern structural/design requirements. GC, poly-T and long
    homopolymers are review flags rather than claims that the guide cannot work.
    """
    checks: List[ValidationCheck] = []
    spacer = (guide.spacer or "").upper()

    _add(
        checks,
        "Spacer length",
        "PASS" if len(spacer) == 20 else "FAIL",
        f"{len(spacer)} nt",
        "SpCas9 targeting in this project uses a 20-nt spacer.",
    )
    canonical = bool(re.fullmatch(r"[ACGT]{20}", spacer))
    _add(
        checks,
        "Canonical DNA spacer",
        "PASS" if canonical else "FAIL",
        "ACGT only" if canonical else "Ambiguous/invalid bases detected",
        "A production candidate should not contain N or other ambiguous bases.",
    )

    covered = guide.genes_covered
    _add(
        checks,
        "All requested genes covered",
        "PASS" if covered == expected_genes else "FAIL",
        f"{covered}/{expected_genes}",
        "A single multi-gene guide is only valid for this design if every requested gene has a PAM-compatible target site.",
    )

    pam_ok = bool(guide.matches) and all(re.fullmatch(r"[ACGTN]GG", m.pam.upper()) for m in guide.matches)
    _add(
        checks,
        "NGG PAM at every target",
        "PASS" if pam_ok else "FAIL",
        "All NGG" if pam_ok else "One or more targets lack NGG",
        "Every per-gene target used by the shared design must be PAM-compatible.",
    )

    unique_gene_matches = len({m.gene for m in guide.matches})
    _add(
        checks,
        "One validated target per gene",
        "PASS" if unique_gene_matches == expected_genes == len(guide.matches) else "FAIL",
        f"{unique_gene_matches} distinct gene match(es)",
        "The validation uses one best PAM-compatible target match for each requested gene.",
    )
    if expected_gene_ids is not None:
        identity_ok = {m.gene for m in guide.matches} == set(expected_gene_ids)
        _add(checks, "Requested gene identities", "PASS" if identity_ok else "FAIL",
             identity_ok, "Target names must equal the requested set; matching counts alone are insufficient.")
    evidence_ok = canonical and bool(guide.matches)
    for m in guide.matches:
        if not re.fullmatch(r"[ACGT]{20}", m.target_spacer):
            evidence_ok = False
            continue
        if canonical:
            pos = mismatch_positions(spacer, m.target_spacer)
            evidence_ok &= (pos == m.mismatch_positions and len(pos) == m.mismatches
                            and seed_mismatch_count(pos) == m.seed_mismatches
                            and compatibility_proxy(spacer, m.target_spacer) == m.compatibility_proxy)
    _add(checks, "Recomputed target evidence", "PASS" if evidence_ok else "FAIL", evidence_ok,
         "Mismatch fields and compatibility are recomputed from the spacer and target sequences.")
    provisional = [w for w in (input_warnings or []) if any(term in w.lower() for term in
                   ("spliced cds", "no cds", "not automatically restricted", "user supplied"))]
    if provisional:
        _add(checks, "Input region context", "REVIEW", "Review source annotation",
             "Confirm genomic continuity and coding context for these inputs before prioritization.")

    if guide.design_type == "Exact shared":
        exact_ok = bool(guide.matches) and all(m.mismatches == 0 for m in guide.matches)
        _add(
            checks,
            "Exact-shared identity",
            "PASS" if exact_ok else "FAIL",
            f"{guide.exact_gene_count}/{expected_genes} exact",
            "An Exact shared design must use the identical 20-nt spacer in every requested gene.",
        )
    else:
        mm_ok = bool(guide.matches) and all(m.mismatches <= max_mismatches_per_gene for m in guide.matches)
        seed_ok = bool(guide.matches) and all(m.seed_mismatches <= max_seed_mismatches_per_gene for m in guide.matches)
        compat_ok = bool(guide.matches) and all(m.compatibility_proxy >= min_compatibility for m in guide.matches)
        _add(
            checks,
            "Per-gene mismatch limit",
            "PASS" if mm_ok else "FAIL",
            f"worst={guide.worst_mismatch_count}; limit={max_mismatches_per_gene}",
            "Every gene must remain within the user-selected total mismatch limit.",
        )
        worst_seed = max((m.seed_mismatches for m in guide.matches), default=0)
        _add(
            checks,
            "PAM-proximal seed mismatch limit",
            "PASS" if seed_ok else "FAIL",
            f"worst={worst_seed}; limit={max_seed_mismatches_per_gene}",
            "PAM-proximal mismatches are treated more conservatively by this tool.",
        )
        _add(
            checks,
            "Minimum compatibility proxy",
            "PASS" if compat_ok else "FAIL",
            f"minimum={guide.minimum_compatibility:.1f}; threshold={min_compatibility:.1f}",
            "This is the app's transparent mismatch-ranking proxy, not a calibrated cleavage probability.",
        )

    gc = gc_percent(spacer) if canonical else 0.0
    gc_status = "PASS" if 40.0 <= gc <= 60.0 else "REVIEW"
    gc_note = "Within the preferred 40–60% review band." if gc_status == "PASS" else "Outside the preferred 40–60% band; do not reject automatically, but review experimentally."
    _add(checks, "GC content", gc_status, f"{gc:.1f}%", gc_note)

    poly_t = "TTTT" in spacer
    _add(
        checks,
        "Poly-T motif",
        "REVIEW" if poly_t else "PASS",
        "TTTT present" if poly_t else "No TTTT",
        "A TTTT motif can be problematic for common Pol III/U6 expression systems; relevance depends on the expression design.",
    )
    long_homopolymer = bool(re.search(r"A{5,}|C{5,}|G{5,}|T{5,}", spacer))
    _add(
        checks,
        "Long homopolymer",
        "REVIEW" if long_homopolymer else "PASS",
        "≥5-base run present" if long_homopolymer else "No ≥5-base run",
        "Long homopolymers are flagged for synthesis/expression review rather than treated as an absolute failure.",
    )

    if panel_hits is None:
        specificity_status = "NOT SCREENED"
        _add(
            checks,
            "Local specificity panel",
            "INFO",
            "Not screened",
            "Optional supplied-FASTA screening has not been run. Whole-genome specificity still requires a genome-aware external analysis.",
        )
    else:
        hits = list(panel_hits)
        if isinstance(panel_hits, PanelScreen):
            _add(checks, "Panel scan scope", "INFO",
                 f"{panel_hits.total_hits} total; {len(hits)} displayed; radius {panel_hits.mismatch_radius}; "
                 f"{panel_hits.panel_bp} bp; truncated={panel_hits.truncated}",
                 f"NGG only, substitutions only. Panel SHA-256: {panel_hits.panel_sha256}")
        if hits:
            specificity_status = "REVIEW"
            closest = min(h.mismatches for h in hits)
            _add(
                checks,
                "Local specificity panel",
                "REVIEW",
                f"{len(hits)} hit(s); closest={closest} mismatch(es)",
                "Classify these hits as intended paralog targets versus unwanted sites. A supplied panel is not a whole-genome screen.",
            )
        else:
            specificity_status = "PASS (SUPPLIED PANEL)"
            _add(
                checks,
                "Local specificity panel",
                "PASS",
                "No near matches in supplied panel",
                "No PAM-compatible near matches were found within the selected mismatch radius in the supplied FASTA panel.",
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
