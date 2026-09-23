# ============================================================================
# TEST REVIEWER EDGE CASES
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains automated tests that protect the Plant MultiGene gRNA Designer workflow from accidental behaviour changes.
#
# HOW TO READ THIS FILE:
# 1. Tests first prepare an input or fixture.
# 2. The relevant program function is called.
# 3. Assertions check that the result still matches the expected behaviour.
# 4. Test logic and expected scientific results are intentionally unchanged.
#
# MAIN TOP-LEVEL PARTS:
# - function: test_no_valid_pam_in_any_gene
# - function: test_iupac_in_candidate_window_is_never_silently_used
# - function: test_single_gene_input_is_rejected_for_shared_design
# - function: test_zero_sequence_overlap_returns_no_shared_guide
# - function: test_exon_junction_is_not_created_by_concatenation
# ============================================================================

import pytest
from plant_multiguide import scan_spcas9, scan_gene_segments, design_shared_guides, clean_dna


# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_no_valid_pam_in_any_gene
# ----------------------------------------------------------------------------
def test_no_valid_pam_in_any_gene():
    genes={"G1":[("e1","A"*80)],"G2":[("e1","T"*80)]}
    sites,guides=design_shared_guides(genes, include_mismatch_aware=True)
    assert all(len(x)==0 for x in sites.values())
    assert guides==[]

@pytest.mark.parametrize("code", ["N","R","Y"])
def test_iupac_in_candidate_window_is_never_silently_used(code):
    raw="A"*10 + code + "A"*9 + "AGG"
    sites=scan_spcas9(raw,"G","e")
    assert not any(s.start==1 for s in sites)
    assert "N" in clean_dna(raw)


# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_single_gene_input_is_rejected_for_shared_design
# ----------------------------------------------------------------------------
def test_single_gene_input_is_rejected_for_shared_design():
    with pytest.raises(ValueError, match="at least two"):
        design_shared_guides({"G1":[("e1","A"*20+"AGG")]})


# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_zero_sequence_overlap_returns_no_shared_guide
# ----------------------------------------------------------------------------
def test_zero_sequence_overlap_returns_no_shared_guide():
    genes={"G1":[("e1","ACGTACGTACGTACGTACGAAGG")],"G2":[("e1","TTGCTTGCTTGCTTGCTTGCTGG")]}
    _,guides=design_shared_guides(genes, include_mismatch_aware=True, max_mismatches_per_gene=1, max_seed_mismatches_per_gene=0)
    assert guides==[]


# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_exon_junction_is_not_created_by_concatenation
# ----------------------------------------------------------------------------
def test_exon_junction_is_not_created_by_concatenation():
    sites=scan_gene_segments({"G":[("ex1","A"*20),("ex2","AGG"+"C"*30)]})["G"]
    assert not any(s.spacer=="A"*20 and s.pam=="AGG" for s in sites)
