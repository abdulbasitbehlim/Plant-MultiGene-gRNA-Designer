import json
from pathlib import Path

from plant_multiguide import design_shared_guides

PUBLISHED_GUIDE = "CTCTACTTTCTCCCTCATCT"


def _fixture_segments():
    p = Path(__file__).resolve().parents[1] / "benchmarks/data/multiknock_pup_fixture.fasta"
    seqs = {}
    name = None
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return {gene: [("benchmark_locus", seq)] for gene, seq in seqs.items()}


def test_hu_validation_reference_is_locked_and_not_mislabeled_as_library_guide():
    p = Path(__file__).resolve().parents[1] / "benchmarks/data/multiknock_pup_reference.json"
    d = json.loads(p.read_text())
    assert d["published_guide"] == PUBLISHED_GUIDE
    assert {g["locus"] for g in d["genes"]} == {"AT4G18197", "AT4G18195", "AT4G18205"}
    assert d["doi"] == "10.1038/s41477-023-01374-4"
    assert d["library_membership_verified"] is False
    assert d["library_rediscovery_claim_allowed"] is False
    assert "validation construct" in d["evidence_class"]


def test_designer_recovers_exact_published_validation_spacer_without_receiving_it():
    _, guides = design_shared_guides(
        _fixture_segments(),
        include_mismatch_aware=True,
        max_mismatches_per_gene=2,
        max_seed_mismatches_per_gene=1,
        min_compatibility=0.0,
        max_results=100,
    )
    matches = [g for g in guides if g.spacer == PUBLISHED_GUIDE]
    assert len(matches) == 1
    assert guides.index(matches[0]) == 0
    assert matches[0].genes_covered == 3
    assert matches[0].exact_gene_count == 2
    assert matches[0].worst_mismatch_count == 1
