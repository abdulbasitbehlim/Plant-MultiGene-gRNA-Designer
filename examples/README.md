# Runnable synthetic examples

These short sequences demonstrate software behavior only; they are not plant loci or experimentally validated guides. Select Manual multi-FASTA and paste a file's contents.

| File | Settings | Expected lesson |
|---|---|---|
| exact_three_genes.fasta | Mismatch mode off | One spacer is present beside NGG in all three inputs |
| mismatch_pair.fasta | Default mismatch mode | A one-distal-mismatch candidate is retained |
| separate_exons.fasta | Group GeneID\|SegmentID enabled | Short pieces in G1 cannot create a joined genomic target |
| no_shared_guide.fasta | Mismatch mode off | No single exact shared guide; a two-guide exact fallback covers both inputs |
| pam_loss.fasta | Either mode | G2 lacks a compatible PAM; sequence similarity alone is insufficient |

Run `python benchmarks/run_audit_examples.py` for all 13 checked scenarios and their machine-readable results.
