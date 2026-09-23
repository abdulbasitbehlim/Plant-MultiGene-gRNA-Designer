# ============================================================================
# TEST VALIDATION
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains automated tests that protect the Plant MultiGene gRNA
# Designer from accidental behaviour changes.
#
# HOW TO READ THIS FILE:
# 1. A test prepares sample input.
# 2. It calls the function being checked.
# 3. Assertions compare the result with the expected behaviour.
# 4. Test expectations are kept unchanged; only explanatory structure is added.
#
# MAIN TOP-LEVEL PARTS:
# - function: _match
# - function: test_exact_shared_validation_passes_core_rules
# - function: test_missing_gene_is_hard_failure
# - function: test_poly_t_is_review_not_hard_failure
# ============================================================================

from plant_multiguide import GeneMatch, SharedGuide
from validation import validate_shared_guide


def _match(gene, mm=0, seed=0, compat=100.0, pam="AGG", target="ACGT"*5):
    return GeneMatch(gene, target, pam, gene, "+", 1, 20, mm, tuple(range(1, mm+1)), seed, compat)


def test_exact_shared_validation_passes_core_rules():
    g = SharedGuide(spacer="ACGT"*5, design_type="Exact shared", matches=[_match("g1"), _match("g2")], gc_percent=50.0, sequence_score=80.0, average_compatibility=100.0, minimum_compatibility=100.0, exact_gene_count=2, worst_mismatch_count=0)
    r = validate_shared_guide(g, expected_genes=2)
    assert r.fail_count == 0
    assert r.status == "PASS"
    assert r.specificity_status == "NOT SCREENED"


def test_missing_gene_is_hard_failure():
    g = SharedGuide(spacer="ACGT"*5, design_type="Exact shared", matches=[_match("g1")], gc_percent=50.0, sequence_score=80.0, average_compatibility=100.0, minimum_compatibility=100.0, exact_gene_count=1, worst_mismatch_count=0)
    r = validate_shared_guide(g, expected_genes=2)
    assert r.status == "FAIL"
    assert any(c.check == "All requested genes covered" and c.status == "FAIL" for c in r.checks)


def test_poly_t_is_review_not_hard_failure():
    spacer = "TTTT" + "ACGT"*4
    g = SharedGuide(spacer=spacer, design_type="Exact shared", matches=[_match("g1", target=spacer), _match("g2", target=spacer)], gc_percent=40.0, sequence_score=55.0, average_compatibility=100.0, minimum_compatibility=100.0, exact_gene_count=2, worst_mismatch_count=0)
    r = validate_shared_guide(g, expected_genes=2)
    assert r.fail_count == 0
    assert r.status == "REVIEW"
