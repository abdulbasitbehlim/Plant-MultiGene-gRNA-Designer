"""Network-enabled exact-spacer recovery rerun for the Hu et al. validation construct.

Retrieves current versioned NCBI records, independently generates guide candidates, and
only then compares generated spacers with the published 20-nt sequence. This is not a
transportome-library membership test.
"""
from pathlib import Path
import json, sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sequence_sources import fetch_ncbi_gene
from plant_multiguide import design_shared_guides

GUIDE = "CTCTACTTTCTCCCTCATCT"
GENES = [("PUP7", "Arabidopsis thaliana"), ("PUP8", "Arabidopsis thaliana"), ("AT4G18205", "Arabidopsis thaliana")]
records = {g: fetch_ncbi_gene(g, org) for g, org in GENES}
segments = {g: r.segments for g, r in records.items()}
_, guides = design_shared_guides(
    segments,
    include_mismatch_aware=True,
    max_mismatches_per_gene=2,
    max_seed_mismatches_per_gene=1,
    min_compatibility=0.0,
    max_results=1000,
)
match = next((g for g in guides if g.spacer == GUIDE), None)
rank = next((i + 1 for i, g in enumerate(guides) if g.spacer == GUIDE), None)
result = {
    "published_spacer": GUIDE,
    "evidence_class": "published follow-up validation construct; library membership unverified",
    "exact_spacer_recovered": bool(match),
    "generated_rank": rank,
    "genes_with_exact_target_sequence": match.exact_gene_count if match else None,
    "genes_with_pam_compatible_target": match.genes_covered if match else None,
    "worst_mismatch_count": match.worst_mismatch_count if match else None,
    "library_guide_rediscovery_claim": False,
    "provenance": {g: r.provenance() for g, r in records.items()},
}
out = Path(__file__).resolve().parent / "results/live_full_gene_result.json"
out.write_text(json.dumps(result, indent=2, default=list))
print(json.dumps(result, indent=2, default=list))
