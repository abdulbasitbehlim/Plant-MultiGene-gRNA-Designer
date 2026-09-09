"""Deterministic Hu et al. 2023 validation-construct sequence regression.

The designer receives only bundled PUP7/PUP8/PUP21 locus contexts. The published
20-nt validation spacer is used only after candidate generation to test exact sequence
recovery. This script does not claim the spacer is a verified member of the 5,635-guide
transportome library.
"""
from pathlib import Path
import csv
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plant_multiguide import design_shared_guides

PUBLISHED_GUIDE = "CTCTACTTTCTCCCTCATCT"
HERE = Path(__file__).resolve().parent

seqs = {}
name = None
for line in (HERE / "data/multiknock_pup_fixture.fasta").read_text().splitlines():
    if line.startswith(">"):
        name = line[1:]
        seqs[name] = ""
    else:
        seqs[name] += line.strip()

segments = {gene: [("benchmark_locus", seq)] for gene, seq in seqs.items()}
_, guides = design_shared_guides(
    segments,
    include_mismatch_aware=True,
    max_mismatches_per_gene=2,
    max_seed_mismatches_per_gene=1,
    min_compatibility=0.0,
    max_results=100,
)

rank = next((i + 1 for i, g in enumerate(guides) if g.spacer == PUBLISHED_GUIDE), None)
matched = next((g for g in guides if g.spacer == PUBLISHED_GUIDE), None)
exact_recovered = int(matched is not None)
row = {
    "benchmark": "Hu 2023 PUP7/PUP8/PUP21 validation-construct sequence regression",
    "evidence_class": "published follow-up validation construct; transportome-library membership unverified",
    "published_spacers_tested": 1,
    "exact_spacer_recovered": exact_recovered,
    "exact_spacer_recovery_percent": 100.0 * exact_recovered,
    "published_spacer": PUBLISHED_GUIDE,
    "generated_rank": rank if rank is not None else "not found",
    "total_candidates_generated": len(guides),
    "genes_with_exact_target_sequence": matched.exact_gene_count if matched else "NA",
    "genes_with_pam_compatible_target": matched.genes_covered if matched else "NA",
    "worst_mismatch_count": matched.worst_mismatch_count if matched else "NA",
    "library_membership_verified": False,
    "library_guide_rediscovery_claim": False,
    "rediscovery_definition": "exact spacer recovery requires identical 20-nt sequence; mismatch tolerance is reported separately",
    "scope": "deterministic locus-level validation-construct regression; not a verified 5,635-guide library benchmark",
}

(HERE / "results").mkdir(exist_ok=True)
with (HERE / "results/locked_benchmark.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=row.keys())
    writer.writeheader()
    writer.writerow(row)
(HERE / "results/locked_benchmark.json").write_text(json.dumps(row, indent=2))
print(json.dumps(row, indent=2))
