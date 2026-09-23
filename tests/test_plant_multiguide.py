from plant_multiguide import (
    clean_dna,
    reverse_complement,
    scan_spcas9,
    design_shared_guides,
    exact_shared_guides,
    mismatch_positions,
    seed_mismatch_count,
    screen_reference_panel,
)


def flank(spacer, pam="AGG"):
    return "TTTTT" + spacer + pam + "AAAAA"


def test_clean_dna_and_reverse_complement():
    assert clean_dna(">x\nacgu rys") == "ACGTNNN"
    assert reverse_complement("ACGTN") == "NACGT"


def test_scan_plus_and_reverse_sites():
    spacer = "ACGTACGTACGTACGTACGA"
    plus = scan_spcas9(flank(spacer), "G1", "exon1")
    assert any(s.spacer == spacer and s.strand == "+" for s in plus)

    desired = "GCTAGCTAGCTAGCTAGCTA"
    downstream = reverse_complement(desired)
    seq = "AAAA" + "CCA" + downstream + "TTTT"
    rev = scan_spcas9(seq, "G2", "exon1")
    assert any(s.spacer == desired and s.strand == "-" and s.pam.endswith("GG") for s in rev)


def test_exact_shared_guide_across_three_genes():
    spacer = "ACGTACGTACGTACGTACGA"
    genes = {
        "G1": [("e1", flank(spacer, "AGG"))],
        "G2": [("e2", "CCCC" + spacer + "TGG" + "AAAA")],
        "G3": [("e3", "GGGG" + spacer + "CGG" + "TTTT")],
    }
    sites, guides = design_shared_guides(genes, include_mismatch_aware=False)
    assert all(sites[g] for g in genes)
    assert guides
    g = guides[0]
    assert g.spacer == spacer
    assert g.design_type == "Exact shared"
    assert g.genes_covered == 3
    assert g.exact_gene_count == 3


def test_mismatch_aware_consensus_targets_all_genes():
    s1 = "ACGTACGTACGTACGTACGA"
    s2 = "ACGTTCGTACGTACGTACGA"  # position 5 mismatch
    s3 = "ACGTACGTACGTACATACGA"  # distal/middle mismatch
    genes = {
        "G1": [("e1", flank(s1, "AGG"))],
        "G2": [("e1", flank(s2, "TGG"))],
        "G3": [("e1", flank(s3, "CGG"))],
    }
    _, guides = design_shared_guides(
        genes,
        include_mismatch_aware=True,
        max_mismatches_per_gene=2,
        max_seed_mismatches_per_gene=1,
        min_compatibility=50,
    )
    assert guides
    assert any(g.genes_covered == 3 for g in guides)
    assert any(g.design_type == "Mismatch-aware consensus" for g in guides)


def test_mismatch_helpers():
    a = "A" * 20
    b = "A" * 12 + "C" + "A" * 7
    pos = mismatch_positions(a, b)
    assert pos == (13,)
    assert seed_mismatch_count(pos) == 1


def test_no_forced_result_for_unrelated_sequences():
    genes = {
        "G1": [("e1", flank("ACGTACGTACGTACGTACGA"))],
        "G2": [("e1", flank("TTGCTTGCTTGCTTGCTTGC"))],
    }
    _, guides = design_shared_guides(
        genes, include_mismatch_aware=True, max_mismatches_per_gene=1, max_seed_mismatches_per_gene=0
    )
    assert guides == []


def test_reference_panel_screen():
    guide = "ACGTACGTACGTACGTACGA"
    panel = {"contig1": flank(guide), "contig2": "A" * 80}
    hits = screen_reference_panel(guide, panel, max_mismatches=0)
    assert hits
    assert hits[0].contig == "contig1"
    assert hits[0].mismatches == 0


def test_ambiguity_is_skipped_in_spacer_but_allowed_at_degenerate_pam_base():
    from plant_multiguide import ambiguity_summary, scan_spcas9
    spacer = "ACGTACGTACGTACGTACGA"
    amb = ambiguity_summary(">x\nACGTRYNN")
    assert amb["count"] == 4
    assert set(amb["codes"]) == {"N", "R", "Y"}

    sites = scan_spcas9("AAAAA" + spacer + "NGG" + "AAAA", "G1", "seg")
    assert any(s.spacer == spacer for s in sites)

    bad_spacer = spacer[:10] + "N" + spacer[11:]
    sites2 = scan_spcas9("AAAAA" + bad_spacer + "AGG" + "AAAA", "G1", "seg")
    assert not any(s.start == 6 for s in sites2)


def test_search_diagnostics_state_seed_exhaustive_strategy():
    from plant_multiguide import design_shared_guides, estimate_search_diagnostics
    s1 = "ACGTACGTACGTACGTACGA"
    s2 = "ACGTTCGTACGTACGTACGA"
    genes = {"G1": [("e1", flank(s1))], "G2": [("e1", flank(s2))]}
    sites, _ = design_shared_guides(genes, include_mismatch_aware=False)
    d = estimate_search_diagnostics(sites)
    assert d.genes == 2
    assert d.unique_seed_spacers >= 2
    assert "every distinct observed" in d.strategy
    assert "not exhaustive over all 4^20" in d.strategy


def test_mismatch_search_workload_guard_is_enforced():
    import pytest
    from plant_multiguide import design_shared_guides
    genes = {
        "G1": [("e1", flank("ACGTACGTACGTACGTACGA"))],
        "G2": [("e1", flank("ACGTTCGTACGTACGTACGA"))],
    }
    with pytest.raises(ValueError, match="workload is too large"):
        design_shared_guides(genes, include_mismatch_aware=True, max_pair_comparisons=1)


def test_hard_gene_count_ceiling():
    import pytest
    from plant_multiguide import design_shared_guides, MAX_GENES
    spacer = "ACGTACGTACGTACGTACGA"
    genes = {f"G{i}": [("e1", flank(spacer))] for i in range(MAX_GENES + 1)}
    with pytest.raises(ValueError, match="at most"):
        design_shared_guides(genes, include_mismatch_aware=False)
