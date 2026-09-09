# Plant MultiGene gRNA Designer v1.3.1 — architecture and scientific design

## 1. Scope

This Streamlit research tool asks a specific question: **can one SpCas9 spacer target every gene entered by the user?** It supports exact shared guides and a constrained mismatch-aware consensus mode for homologous plant genes. It never forces a single-guide answer when no candidate satisfies the rules.

A software PASS is a design-rule result, not proof of editing efficiency, phenotype, or genome-wide specificity.

## 2. Modules

```text
app.py
  ├─ plant_multiguide.py   # scanning, exact intersection, consensus search, ranking, panel scan
  ├─ sequence_sources.py   # NCBI / Ensembl / manual FASTA + reproducibility metadata
  └─ validation.py         # PASS / REVIEW / FAIL checks

tests/
  ├─ scientific-core tests
  ├─ sequence/provenance tests
  ├─ validation tests
  └─ static Streamlit UI contract test
```

The scientific core does not depend on Streamlit.

## 3. Sequence acquisition and reproducibility

`GeneSequenceRecord` stores the exact sequence segments used plus:

- source and accession;
- source record/version;
- assembly or genomic accession when available;
- annotation release/date when available;
- UTC retrieval timestamp;
- SHA-256 fingerprint of the exact segment set used by the run;
- ambiguity count/codes;
- source-specific warnings.

### NCBI

The selected RefSeq accession includes its version when NCBI supplies one (for example `NM_... .#`). The tool also attempts to capture the linked chromosome/genomic accession version from NCBI Gene and the GenBank/RefSeq record date. Transcript-level retrieval does **not** by itself prove a particular reference assembly, so unresolved assembly context remains an explicit warning/review item.

### Ensembl

The tool records the versioned transcript ID when available, the `assembly_name` returned by lookup, and the current Ensembl release reported by the REST service. Because public annotation can change, the sequence SHA-256 is retained in exports so a later rerun can verify whether the exact sequence input changed.

### Manual FASTA

Manual inputs are labelled `user-supplied / unspecified` for assembly unless the user encodes assembly/genotype information in the labels. The exact normalized sequences are fingerprinted.

## 4. Ambiguous IUPAC bases

Ambiguity is handled explicitly rather than silently as mismatches:

1. `R, Y, S, W, K, M, B, D, H, V, X` are normalized to `N`.
2. A candidate **spacer containing N is skipped**.
3. An ambiguous base is accepted only at the degenerate `N` position of an otherwise resolved `NGG`/`CCN` PAM; ambiguity in required G/G or C/C positions does not create a PAM.
4. Ambiguous bases are never converted into ordinary spacer mismatches for consensus scoring.
5. Sequence provenance reports ambiguity count and original ambiguity codes.

This behavior is covered by unit tests.

## 5. SpCas9 site discovery

Each sequence segment is scanned independently on both strands for 20-nt spacers adjacent to NGG PAMs. Separate exon/CDS segments are not concatenated, preventing artificial exon-junction targets.

## 6. Exact shared guides

For each gene, PAM-compatible spacers are indexed. The exact mode computes the set intersection across **all requested genes**. A reported exact guide therefore has the identical 20-nt spacer next to an NGG-compatible PAM in every gene.

## 7. Mismatch-aware consensus search — exact search semantics

The v1.3.1 algorithm is **seed-exhaustive over observed PAM-compatible spacers**, not “first-hit greedy” and not exhaustive over all theoretical 20-mers.

For every distinct spacer observed at a PAM-compatible site in any requested gene:

1. use that observed spacer as a seed;
2. inspect every PAM-compatible site in each gene and retain the nearest target under the weighted mismatch distance;
3. reject the seed if any gene already violates total/seed mismatch limits;
4. create a column-wise majority consensus from those selected targets, using the seed only to break nucleotide ties;
5. re-scan every gene against that consensus and retain the best PAM-compatible target;
6. reject the proposal if any gene violates mismatch, seed-mismatch, or minimum compatibility limits;
7. deduplicate identical consensus spacers;
8. **rank every surviving unique proposal before `max_results` truncation**.

Therefore, “best” means **best among all generated seed-derived proposals under the implemented rank key**. It is not a mathematical proof of the globally optimal spacer over the full `4^20` synthetic sequence space. This distinction is shown in the UI and exports.

The compatibility value is an explainable mismatch-position proxy, not CRISPys, CFD, MOFF, or a calibrated cleavage probability.

## 8. Computational complexity and hard ceilings

Let:

- `U` = number of distinct observed PAM-compatible seed spacers;
- `S` = total PAM-compatible sites across all entered genes.

One nearest-neighbour pass is approximately `O(U × S)`. The algorithm can perform two complete passes (seed selection + consensus re-optimization), so the tool reports a conservative upper workload of approximately `2 × U × S` pair comparisons.

Interactive safeguards:

- **recommended:** 2–12 genes;
- **hard gene cap:** 20 genes per run;
- **hard mismatch-search workload guard:** 5,000,000 estimated pair comparisons.

The guard is based on actual PAM-site counts, so a 15–20 gene run can still work when the scanned regions are small, while long/permissive inputs can be stopped earlier. Large genome-scale library design should use a dedicated scalable workflow such as the approaches represented by Multi-Knock/CRISPys rather than this interactive application.

## 9. Ranking

Surviving candidates are ranked by:

1. genes covered;
2. exact-match gene count;
3. minimum compatibility across genes;
4. average compatibility;
5. worst mismatch count;
6. transparent sequence-quality score;
7. closeness of GC to 50%.

Exact-shared designs are prioritized over mismatch-aware designs in the combined output.

## 10. Validation layer

`validation.py` separates hard requirements from review flags.

Hard failures include invalid spacer length/bases, incomplete gene coverage, missing NGG compatibility, or violation of user-selected mismatch constraints. Review flags include GC outside the preferred band, poly-T, long homopolymers, and findings in a supplied local panel.

The custom-guide validator evaluates an external 20-nt spacer against the already scanned PAM-compatible sites in every loaded gene.

## 11. Local reference-panel screen

The optional panel scanner is intended for a user-supplied paralog/off-target panel. It is **not** genome-wide analysis. Intended family targets must be distinguished from unwanted sites by the user.

## 12. Streamlit layer and exports

The UI exposes:

- exact versus mismatch-aware design controls;
- provenance/reproducibility metadata;
- ambiguity warnings;
- search strategy and workload diagnostics;
- ranked candidate and per-gene evidence tables;
- PASS/REVIEW/FAIL validation;
- custom-guide validation;
- CSV/FASTA/JSON exports.

JSON exports include full provenance and search diagnostics.

## 13. Testing status

v1.3.1 currently passes **34 automated tests** with **85% branch-aware core coverage** covering core scanning, exact and mismatch-aware design, workload guard, gene ceiling, ambiguity rules, reference-panel behavior, sequence-source parsing, provenance fingerprints, and validation.

`tests/test_ui_contract.py` also parses `app.py` with Python AST and checks that the critical reviewer-facing UI elements are present. **This is a static UI contract test, not an interactive Streamlit runtime test.** The artifact environment used to build this package does not contain Streamlit and cannot download it from PyPI, so an actual browser/runtime launch could not be performed here. Run `streamlit run app.py` locally and use `UI_TEST_CHECKLIST.md` before calling the interface release-tested.

## 14. Scientific boundaries

- Homologous genes can still differ in chromatin, expression, alleles, or repair outcomes.
- The current consensus proxy is not an experimentally trained efficacy model.
- Annotation provenance improves reproducibility but does not replace mapping the final guide to the exact reference/cultivar assembly.
- A local FASTA panel is not a genome-wide off-target screen.
- Large gene-family library construction is outside the intended interactive scale.

## 15. Key literature

- Prykhozhij et al. CRISPR MultiTargeter. *PLoS One*. 2015. doi:10.1371/journal.pone.0119372.
- Hyams et al. CRISPys. *J Mol Biol*. 2018. doi:10.1016/j.jmb.2018.03.019.
- Hu et al. Multi-Knock. *Nature Plants*. 2023. doi:10.1038/s41477-023-01374-4.
- Berman et al. Multi-targeted CRISPR libraries in tomato. *Nature Communications*. 2025. doi:10.1038/s41467-025-59280-6.

## 17. Published validation-construct regression and comparator protocol

The PUP benchmark was corrected in v1.3.1 to separate **screen-library evidence** from **follow-up validation evidence**. Hu et al.'s Figure 5 legend states that the CR8/21 screen sgRNA targets PUP8/PUP21 and was not designed for PUP7 because PUP7 lacks the corresponding PAM for that guide. The Methods separately describe `CTCTACTTTCTCCCTCATCT` as a 20-nt protospacer picked to target PUP7/PUP8/PUP21 at once in follow-up CRISPR cloning.

Accordingly, `benchmarks/run_locked_benchmark.py` is now named and reported as a **validation-construct sequence regression**. The design function receives only bundled PUP target contexts. The published spacer is compared with generated output only after candidate generation. In the locked fixture the exact 20-nt spacer is generated at rank 1; its target contexts contain two exact sequence matches and one one-mismatch PAM-compatible target.

The terms are strict:

- **exact spacer recovered** = the generated 20-nt sequence is identical to the published 20-nt sequence;
- **target compatibility** = per-gene PAM, mismatch and seed-mismatch properties, reported separately;
- **library guide rediscovered** = reserved for a future case in which both exact guide sequence and intended target assignment are directly verified from the published library data.

The Nature landing page identifies Supplementary Data 1 as the all-family sgRNA library, but the binary assignment table was not directly retrievable in the artifact environment. v1.3.1 therefore makes **no claim** that the PUP validation spacer belongs to the 5,635-guide transportome library. `benchmarks/BENCHMARK_EVIDENCE_RESOLUTION.md` records the evidence and manuscript-safe wording.

Fresh CRISPys and CRISPR MultiTargeter concordance counts are not fabricated. `benchmarks/EXTERNAL_COMPARISON_PROTOCOL.md` fixes the comparison rules and `benchmarks/concordance.py` computes intersection/ours-only/external-only counts after external outputs are exported with a `guide` column.

## 18. CI, packaging and release reproducibility

- Runtime and developer dependencies are exactly pinned.
- GitHub Actions tests Python 3.11, 3.12 and 3.13 on each push/pull request and enforces at least 80% core coverage.
- `coverage.txt` and `coverage.xml` preserve the measured artifact result.
- `pyproject.toml` installs the scientific core; Docker/Compose provides a one-command application environment.
- `CITATION.cff` and `.zenodo.json` are release metadata templates only; no DOI is claimed until an actual archived release exists.
- `CONTRIBUTING.md`, `SUPPORT.md`, `UI_TEST_CHECKLIST.md` and `RELEASE_CHECKLIST.md` define reviewer/community workflows.
