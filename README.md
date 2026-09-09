# Plant MultiGene gRNA Designer v1.3.1

A Python/Streamlit research tool for asking a specific CRISPR design question: **can one SpCas9 spacer plausibly target every member of a user-defined plant gene set?** The tool searches exact shared targets and mismatch-aware consensus targets, preserves sequence provenance, validates every candidate, and makes its search limits explicit.

## Live app and related projects

**Run the Plant MultiGene app:** https://plant-multigene-grna-designer.onrender.com/

**GitHub repository:** https://github.com/abdulbasitbehlim/Plant-MultiGene-gRNA-Designer

**Also available — OpenCRISPR-1 gRNA Designer:**
- Live app: https://opencrispr1-grna-designer.onrender.com/
- GitHub: https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer

> The Render free tier may sleep after inactivity, so the first load can take longer while the service wakes up.

## Statement of need

Single-gene CRISPR designers are optimized for choosing guides against one locus. Plant functional redundancy often requires a different workflow: a researcher may want one sgRNA that intentionally recognizes several homologous genes. CRISPR MultiTargeter, CRISPys and Multi-Knock established this multi-target design problem, but reproducible interpretation still depends on target provenance, exon boundaries, mismatch assumptions and search semantics. This project provides a small, auditable workflow that keeps those decisions visible instead of reporting an unexplained universal score. It is intended for research prioritization and teaching, not as a substitute for genome-wide off-target analysis or experimental validation.

## What is different

| Capability | This tool | CRISPR MultiTargeter | CRISPys / Multi-Knock |
|---|---|---|---|
| User-defined multi-gene shared-guide search | Yes | Yes | Yes |
| Exact shared-guide mode separated from mismatch-aware mode | Yes | Common-target workflow | Optimization-focused |
| Search semantics documented in code/docs | Seed-exhaustive over observed PAM-compatible spacers; no first-hit stopping | Alignment-oriented | Hierarchical/optimization strategy |
| Independent exon segments to avoid synthetic junction guides | Yes | Input-dependent | Input-dependent |
| Per-run accession/version + assembly/release + SHA-256 sequence fingerprint | Yes | Not the focus of the original publication | Not the focus of the original publication |
| Explicit workload ceiling | Yes | Different implementation | Designed for larger analyses |
| Local validation and custom-guide check | Yes | Tool-specific | Tool-specific |

This table describes design scope, not a claim that this implementation is globally superior.

## Scientific design modes

1. **Exact shared guide** — the identical 20-nt spacer occurs next to an NGG PAM in every requested gene.
2. **Mismatch-aware consensus guide** — a single spacer has a PAM-compatible near-match in every requested gene while satisfying user-selected total-mismatch, PAM-proximal seed-mismatch and compatibility-proxy constraints.

Mismatch-aware search is **seed-exhaustive over every distinct observed PAM-compatible spacer**. It does not stop after the first acceptable seed. Each seed is matched to its nearest target in every gene, a consensus is generated and re-optimized, all surviving unique proposals are ranked, and only then is `max_results` applied. It is **not** exhaustive over all theoretical `4^20` synthetic spacers, so “rank 1” means best among generated proposals rather than a proof of a mathematical global optimum.

## Published validation-construct regression

Hu et al. (*Nature Plants*, 2023; DOI `10.1038/s41477-023-01374-4`) report two distinct PUP guide contexts. The **CR8/21 screen line** is a PUP8/PUP21 guide and its Figure 5 legend explicitly states that it was not designed to target PUP7 because that guide lacks a corresponding PUP7 PAM. Separately, the cloning Methods describe the 20-nt protospacer

`CTCTACTTTCTCCCTCATCT`

as being picked to target **PUP7 / AT4G18197**, **PUP8 / AT4G18195** and **PUP21 / AT4G18205** at once for follow-up validation.

The v1.3.1 benchmark therefore treats this sequence as a **published three-gene validation construct**, not as a verified member of the 5,635-guide transportome library. The publisher identifies Supplementary Data 1 as the all-family sgRNA library, but that binary assignment table was not directly retrieved in this artifact environment, so library membership is deliberately left **unverified**.

The locked generation regression gives the designer only the bundled PUP target contexts. It does not pass the published spacer into the design function. Result:

| Published validation spacers | Exact spacer recovered | Rank | Exact target sequences | PAM-compatible targets | Worst mismatch |
|---:|---:|---:|---:|---:|---:|
| 1 | **1** | **1** | 2/3 | 3/3 | 1 |

**Definition:** “exact spacer recovered” means that the algorithm emitted the identical 20-nt sequence. Per-gene mismatches are a separate target-compatibility measurement and do not convert a near-hit into an exact recovery. The stronger phrase **“library guide rediscovered” is not used** unless both spacer sequence and intended target assignment are directly verified in the published library table. See `benchmarks/BENCHMARK_EVIDENCE_RESOLUTION.md`.

A stronger future library benchmark is the tomato NPF1.10/NPF1.11/NPF1.12 case reported by Berman et al. (2025), where the paper unambiguously states that one sgRNA targets all three genes with zero mismatches. Its exact guide sequence still needs direct Supplementary Data 1 verification before exact-guide concordance is claimed.

## Head-to-head comparison status

A machine-readable concordance workflow is included in `benchmarks/concordance.py` and `benchmarks/EXTERNAL_COMPARISON_PROTOCOL.md`.

| Comparator | Frozen input | Current artifact status | Concordance counts |
|---|---|---|---|
| CRISPys | Hu PUP validation trio | Published-method reference; fresh external rerun still required | **Not claimed** |
| CRISPR MultiTargeter | Hu PUP validation trio | Reproducible comparison protocol prepared; external service not executed in this build | **Not claimed** |

The package intentionally does not invent “both/ours-only/external-only” counts. Before a manuscript submission, run the comparators against the **same frozen sequences, PAM model, target regions and score thresholds**, export a `guide` column, then execute `benchmarks/concordance.py` and commit the resulting table.

## Inputs

### Gene lookup
Enter two or more gene symbols/IDs and one plant organism. NCBI RefSeq mode tries to retain coding exon pieces from a representative versioned RefSeq RNA. Ensembl REST mode scans individual transcript exons. Sequence retrieval records source/version, genomic assembly/accession or assembly name, annotation release/date when available, UTC retrieval time and a SHA-256 fingerprint of the exact segments used.

### Manual multi-FASTA
Provide one FASTA record per gene. This is useful for cultivar-specific sequences or frozen reviewed target sequences. Ambiguous IUPAC bases (`N`, `R`, `Y`, etc.) are normalized to `N`; any spacer window containing ambiguity is skipped rather than silently treated as a mismatch.

## Runnable example

In the Streamlit app choose **Gene lookup**, organism **Arabidopsis thaliana**, and enter:

```text
AT4G18197
AT4G18195
AT4G18205
```

Use NCBI RefSeq or Ensembl, enable mismatch-aware design, and inspect both the guide table and exported provenance. For a publication benchmark, freeze/export the exact retrieved sequences because public annotations can change.

## Validation

Each generated guide receives **PASS / REVIEW / FAIL** with a criterion-by-criterion explanation.

Hard checks include 20 resolved DNA bases, complete coverage of every requested gene, NGG evidence for every selected target, and the selected mismatch/seed/compatibility limits. Review checks include GC range, poly-T, long homopolymers and optional supplied-panel specificity review. A custom 20-nt guide can also be validated against all loaded genes.

A PASS means the candidate satisfies rules implemented here. It does **not** mean experimentally validated editing efficacy.

## Computational ceiling

- Recommended interactive range: **2–12 genes**.
- Hard maximum: **20 genes per run**.
- Mismatch-aware guard: **5,000,000 estimated seed-to-site comparisons**.

The guard prevents permissive searches over large gene lists from expanding without bound. For larger genome-scale library work, use a purpose-built pipeline such as CRISPys/Multi-Knock rather than bypassing the limit.

## Installation

### Clean Python environment

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

All direct runtime dependencies are exactly pinned in `requirements.txt`.

### One-command Docker run

```bash
docker compose up --build
```

Then open `http://localhost:8501`.

### Scientific-core package

```bash
pip install .
```

This installs the importable scientific modules. The full Streamlit application is run from the repository with the pinned requirements above.

## Tests and measured coverage

```bash
pip install -r requirements-dev.txt
pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-report=xml
```

v1.3.1 artifact result: **34/34 tests passed; 85% branch-aware unit coverage** across `plant_multiguide.py`, `sequence_sources.py` and `validation.py`. Coverage reports are committed as `coverage.txt` and `coverage.xml`.

Explicit reviewer edge cases include:

- no valid PAM in any gene;
- `N`, `R` and `Y` ambiguity inside a candidate window;
- single-gene input rejection for a shared-guide design;
- unrelated genes with zero acceptable shared target;
- exon-junction avoidance;
- mocked NCBI/Ensembl version, assembly/release and ambiguity provenance;
- exact recovery of the locked Hu et al. validation-construct spacer, with library membership explicitly left unverified.

`app.py` is excluded from the numerical unit-coverage percentage because UI lines are not meaningfully exercised by core unit tests. `tests/test_ui_contract.py` checks required UI wiring statically, and `UI_TEST_CHECKLIST.md` defines the browser-level release test.

## Continuous integration

`.github/workflows/ci.yml` runs the pinned test environment on Python 3.11, 3.12 and 3.13 for every push and pull request, generates coverage, and enforces a minimum total coverage threshold.

## Scoring boundaries

The multi-gene **compatibility proxy** is intentionally transparent and dependency-light. It is not CFD, MOFF, CRISTA or a calibrated cleavage probability. A final experimental shortlist should be re-evaluated using a genome-aware off-target workflow and, where appropriate, modern empirical specificity/activity models.

## Reproducibility

Exports retain, when available:

- versioned source accession/transcript;
- genomic chromosome/accession or assembly name;
- annotation release or RefSeq/GenBank record date;
- UTC retrieval timestamp;
- exact sequence-set SHA-256;
- ambiguity counts/codes and warnings.

The SHA-256 value allows a later rerun to detect when the biological input sequence has changed even if a public gene identifier is unchanged.

## Citation and archival release

`CITATION.cff` and `.zenodo.json` are included. **No DOI has been minted in this local package.** After final review, create a public tagged release, archive that exact release with Zenodo, then insert the minted DOI into `CITATION.cff`, README and any manuscript. See `RELEASE_CHECKLIST.md`.

## Key literature

- Prykhozhij SV, Rajan V, Gaston D, Berman JN. CRISPR MultiTargeter. *PLoS One*. 2015;10:e0119372. DOI: `10.1371/journal.pone.0119372`.
- Hyams G, et al. CRISPys: Optimal sgRNA Design for Editing Multiple Members of a Gene Family Using the CRISPR System. *Journal of Molecular Biology*. 2018;430:2184–2195. DOI: `10.1016/j.jmb.2018.03.019`.
- Hu Y, et al. Multi-Knock—a multi-targeted genome-scale CRISPR toolbox to overcome functional redundancy in plants. *Nature Plants*. 2023;9:572–587. DOI: `10.1038/s41477-023-01374-4`.
- Berman A, et al. Construction of multi-targeted CRISPR libraries in tomato to overcome functional redundancy at genome-scale level. *Nature Communications*. 2025. DOI: `10.1038/s41467-025-59280-6`.

## License and support

Original code in this repository is MIT licensed. See `LICENSE`, `CONTRIBUTING.md` and `SUPPORT.md`. Research software; no warranty of experimental performance.
