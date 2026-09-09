# User Guide

## Gene lookup workflow

1. Run `streamlit run app.py`.
2. Enter at least two gene symbols/IDs.
3. Enter the plant organism.
4. Choose NCBI RefSeq or Ensembl REST.
5. Keep mismatch-aware consensus enabled if you want one guide even when the homologs are not perfectly identical.
6. Start with conservative settings (for example, no more than 1-2 mismatches per target and no more than one PAM-proximal seed mismatch).
7. Inspect the per-gene match table for every shortlisted guide.
8. Export the candidate list and perform genome-wide specificity review before experimental use.

## Manual FASTA workflow

Use one FASTA record per gene. Prefer reviewed genomic exon/CDS sequences rather than a spliced sequence if you want to eliminate the possibility of exon-junction candidates.

## Interpreting the result

- **Exact shared:** the same 20-nt guide is found next to an NGG PAM in every gene.
- **Mismatch-aware consensus:** the guide has a PAM-compatible near-match in every gene under the selected mismatch constraints.
- **Minimum compatibility:** the weakest per-gene heuristic compatibility among the target genes.
- **Exact genes:** how many of the requested genes match the guide with zero spacer mismatches.
- **Seed mismatches:** mismatches in guide positions 13-20, treated conservatively because these are PAM-proximal.

A high ranking is not experimental validation.


## Validate guides

After design, every guide has a **Validation** column with PASS, REVIEW or FAIL. Select a guide to see the full checklist. Hard failures indicate that the guide does not satisfy the current multi-gene design rules. Review flags identify GC, poly-T, homopolymer or supplied-panel findings that need human review.

To check a guide from another source, open **Validate a custom 20-nt guide against these genes**, paste the 20-nt spacer, and click **Validate custom guide**. The app reports the closest PAM-compatible target in every loaded gene and applies the same mismatch, seed and compatibility thresholds.

The optional FASTA panel screen updates the selected guide's specificity validation. Panel screening is local only and does not replace whole-genome analysis.

## v1.3.1 search interpretation

Mismatch-aware mode does **not** stop at the first acceptable consensus. Every distinct PAM-compatible spacer actually observed in the loaded genes is used as a seed. For each seed, the tool selects the nearest PAM-compatible target in every gene, builds a majority consensus, re-optimizes that consensus across every gene, and then ranks all surviving unique proposals before returning the requested top results.

This is exhaustive across **observed seed spacers**, not across every theoretical 20-nt sequence. Therefore the first result is the best candidate under this generated-proposal search, not a proof of the global optimum over all possible synthetic spacers.

### Scale limits

- Recommended interactive use: **2–12 genes**.
- Hard cap: **20 genes**.
- Mismatch-aware workload guard: **5,000,000 estimated pair comparisons**.

The app displays the PAM-site count, unique seed count and estimated upper comparison workload. For larger library-scale designs, use a dedicated scalable multi-target workflow.

## Reproducibility and provenance

Open **Sequence provenance, reproducibility and warnings** after a run. Record the following with any guide you keep:

- accession/source record version;
- assembly or genomic accession shown by the tool;
- annotation release/date;
- retrieval timestamp;
- sequence SHA-256 fingerprint.

The SHA-256 identifies the exact normalized sequence segments used by that run and is included in JSON exports. This matters because public records and annotations can change.

## Ambiguous IUPAC sequence

Manual/cultivar/draft FASTA may contain ambiguity codes. v1.3.1 behavior is deliberate:

- ambiguity codes are normalized to `N`;
- a 20-nt spacer containing `N` is skipped;
- ambiguity is not counted as a spacer mismatch;
- an unresolved base may occur only at the degenerate N position of an otherwise resolved NGG/CCN PAM;
- the ambiguity count/codes are shown in provenance.

## UI testing note

The delivered package includes a static UI contract test, but the build environment could not install or launch Streamlit. Before release/deployment, run the local checklist in `UI_TEST_CHECKLIST.md` and record the result.

## Benchmark terminology

For the Hu et al. PUP regression, **exact spacer recovered** means exact 20-nt sequence identity. Mismatch-tolerant targeting of one or more genes is reported separately. The PUP sequence in the locked fixture is treated as a published follow-up validation construct; v1.3.1 does not call it a verified member of the 5,635-guide transportome library. Read `benchmarks/BENCHMARK_EVIDENCE_RESOLUTION.md` before citing the benchmark.

## Reviewer/release workflow (v1.3.1)

For a manuscript-quality result, freeze the biological inputs first. Save the exported versioned accession, assembly/release metadata and sequence SHA-256. Run the locked Multi-Knock regression, then run `benchmarks/live_full_gene.py` with network access on the frozen full records. For external comparisons, export CRISPys/MultiTargeter guide lists from those same frozen sequences and calculate concordance with `benchmarks/concordance.py`. Do not report external concordance values that were not freshly executed.

Automated core QA at artifact build: **34 tests, 85% branch-aware coverage**. The numerical coverage excludes `app.py`; complete `UI_TEST_CHECKLIST.md` in a browser before release.
