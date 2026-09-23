import pytest
from plant_multiguide import scan_spcas9, scan_gene_segments, design_shared_guides, clean_dna

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

def test_single_gene_input_is_rejected_for_shared_design():
    with pytest.raises(ValueError, match="at least two"):
        design_shared_guides({"G1":[("e1","A"*20+"AGG")]})

def test_zero_sequence_overlap_returns_no_shared_guide():
    genes={"G1":[("e1","ACGTACGTACGTACGTACGAAGG")],"G2":[("e1","TTGCTTGCTTGCTTGCTTGCTGG")]}
    _,guides=design_shared_guides(genes, include_mismatch_aware=True, max_mismatches_per_gene=1, max_seed_mismatches_per_gene=0)
    assert guides==[]

def test_exon_junction_is_not_created_by_concatenation():
    sites=scan_gene_segments({"G":[("ex1","A"*20),("ex2","AGG"+"C"*30)]})["G"]
    assert not any(s.spacer=="A"*20 and s.pam=="AGG" for s in sites)
