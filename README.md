# Plant MultiGene gRNA Designer v1.4.0

Plant MultiGene gRNA Designer is a Python and Streamlit research tool for a simple but important CRISPR question:

**Can one SpCas9 guide target several related plant genes at the same time?**

The program is designed for multi-gene targeting, especially when several homologous genes may have overlapping or redundant biological functions.

**Live app:** https://plant-multigene-grna-designer.onrender.com/  
**Repository:** https://github.com/abdulbasitbehlim/Plant-MultiGene-gRNA-Designer

> The free Render service may take a little time to wake up after inactivity.

---

## What this tool does

In simple terms, the program:

1. receives two or more plant gene sequences;
2. finds SpCas9-compatible 20 nt targets next to NGG PAMs;
3. looks for guides that can work across all requested genes;
4. checks exact matches first;
5. can also search mismatch-aware consensus guides;
6. validates every candidate;
7. keeps the exact sequence source and settings used in the run;
8. exports results for later review.

The tool is meant for **research prioritization and teaching**.

It does not replace whole-genome off-target analysis or experimental validation.

---

## Two design modes

### 1. Exact shared guide

The same 20 nt spacer must occur next to an NGG PAM in every requested gene.

This is the simplest case because one identical spacer is present in all genes.

### 2. Mismatch-aware consensus guide

The tool can also search for one spacer that has acceptable PAM-compatible near-matches in every gene.

The user can control things such as:

- total mismatches;
- PAM-proximal seed mismatches;
- compatibility limits.

The search checks every distinct observed PAM-compatible spacer as a possible starting point.

It does not search every theoretical sequence in the full 4^20 sequence space.

---

## Why multi-gene design is useful

In plants, related genes can sometimes perform overlapping functions.

Knocking out only one member of a gene family may therefore produce a weak or incomplete phenotype.

A shared guide can be useful when the scientific goal is to intentionally target several homologous genes together.

This software is designed around that specific use case.

---

## Input options

### Gene lookup

Enter two or more gene symbols or IDs from one plant organism.

The program can retrieve sequence data from supported public biological resources while keeping source and version information when available.

### Accession IDs

You can provide known NCBI nucleotide/RefSeq accessions or Ensembl stable IDs.

### Manual multi-FASTA

You can also paste or upload your own sequences.

Each FASTA record normally represents one gene.

For genes split into separate segments or exons, grouped identifiers can be written in the form:

\`\`\`text
GeneID|SegmentID
\`\`\`

Duplicate identifiers and empty records are rejected.

---

## Example

A simple Arabidopsis example is:

\`\`\`text
AT4G18197
AT4G18195
AT4G18205
\`\`\`

These correspond to the PUP7, PUP8 and PUP21 context used in the project benchmark.

When using database retrieval for publication-quality work, it is a good idea to export and freeze the exact sequences used, because public annotations can change over time.

---

## Validation labels

Each guide receives one of three labels:

- **PASS**
- **REVIEW**
- **FAIL**

These labels describe whether the candidate passed the rules implemented by the software.

They do **not** mean that the guide has been experimentally validated.

### Examples of hard checks

The software checks things such as:

- spacer length;
- resolved DNA bases;
- coverage of every requested gene;
- NGG PAM evidence;
- mismatch limits;
- seed-region mismatch limits.

### Review-level checks

The program also considers practical sequence features such as:

- GC content;
- poly-T sequences;
- long homopolymers;
- optional local specificity-panel results.

---

## Important safety and quality checks

The current version also protects against several common data problems.

- Duplicate FASTA IDs are rejected.
- Mixed-species gene lookup is rejected.
- Duplicate resolved biological records are rejected.
- Ambiguous sequence bases are tracked.
- Candidate windows containing unresolved bases are skipped.
- Run settings are stored so old results are not reused with incompatible settings.
- Large searches are bounded rather than allowed to grow without control.

---

## Computational limits

The interactive tool is intended for relatively small multi-gene problems.

Current limits include:

- recommended range: **2–12 genes**;
- hard maximum: **20 genes**;
- sequence/local panel size: up to **500,000 bases** each;
- mismatch-aware comparison guard: **5,000,000 estimated comparisons**.

For very large genome-scale library design, use a purpose-built large-scale workflow instead of bypassing these limits.

---

## How the ranking works

The tool uses a transparent **compatibility proxy** for multi-gene ranking.

It is intentionally simple and inspectable.

It is **not** the same as:

- CFD;
- MOFF;
- CRISTA;
- a calibrated cleavage probability;
- a guaranteed editing-efficiency model.

A final experimental shortlist should still be checked with appropriate genome-aware specificity and activity tools.

---

## Reproducibility

The software stores as much biological provenance as possible.

Exports may include:

- versioned accessions or transcript IDs;
- chromosome/accession or assembly information;
- annotation version/date;
- retrieval time;
- exact sequence-set SHA-256 fingerprint;
- ambiguity counts and warnings;
- run settings.

The SHA-256 fingerprint is useful because it can reveal when the actual biological input sequence has changed even if the public gene identifier has stayed the same.

---

## Published validation-context regression

The repository includes a regression example based on the PUP7/PUP8/PUP21 context reported by Hu et al. (2023).

The published spacer:

\`\`\`text
CTCTACTTTCTCCCTCATCT
\`\`\`

is used as a known reference context for software testing.

The benchmark checks whether the program can recover that spacer from the bundled target sequences without directly giving the spacer to the design function.

This is a **software regression benchmark**, not proof that the current software itself has been experimentally validated.

See:

- [benchmarks/BENCHMARK_EVIDENCE_RESOLUTION.md](benchmarks/BENCHMARK_EVIDENCE_RESOLUTION.md)
- [docs/RESEARCH_AUDIT.md](docs/RESEARCH_AUDIT.md)

---

## Installation

Clone the repository:

\`\`\`bash
git clone https://github.com/abdulbasitbehlim/Plant-MultiGene-gRNA-Designer.git
cd Plant-MultiGene-gRNA-Designer
\`\`\`

Create a virtual environment:

\`\`\`bash
python -m venv .venv
\`\`\`

Activate it.

### Windows PowerShell

\`\`\`powershell
.venv\Scripts\Activate.ps1
\`\`\`

### Linux/macOS

\`\`\`bash
source .venv/bin/activate
\`\`\`

Install the required packages:

\`\`\`bash
python -m pip install --upgrade pip
pip install -r requirements.txt
\`\`\`

Start the app:

\`\`\`bash
streamlit run app.py
\`\`\`

---

## Docker

You can also run the project with Docker:

\`\`\`bash
docker compose up --build
\`\`\`

Then open:

\`\`\`text
http://localhost:8501
\`\`\`

---

## Scientific-core package

The reusable Python modules can also be installed with:

\`\`\`bash
pip install .
\`\`\`

The full Streamlit interface is still run from the repository.

---

## Main project files

- \`app.py\` — Streamlit user interface.
- \`plant_multiguide.py\` — main multi-gene guide design logic.
- \`sequence_sources.py\` — biological sequence retrieval.
- \`accession_sources.py\` — accession handling.
- \`validation.py\` — guide and input validation.
- \`run_state.py\` — keeps run settings/results consistent.
- \`benchmarks/\` — reproducibility and benchmark scripts.
- \`tests/\` — automated tests.

---

## Running the tests

Install development requirements and run:

\`\`\`bash
pip install -r requirements-dev.txt
pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-report=xml
\`\`\`

The repository also contains continuous-integration checks for supported Python versions.

---

## Additional documentation

For more detailed explanations, see:

- [USER_GUIDE.md](USER_GUIDE.md)
- [docs/RESEARCH_AUDIT.md](docs/RESEARCH_AUDIT.md)
- [examples/README.md](examples/README.md)
- [UI_TEST_CHECKLIST.md](UI_TEST_CHECKLIST.md)
- [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md)

---

## Key references

- Prykhozhij SV, Rajan V, Gaston D, Berman JN. CRISPR MultiTargeter. *PLoS One*. 2015. DOI: 10.1371/journal.pone.0119372
- Hyams G, et al. CRISPys. *Journal of Molecular Biology*. 2018. DOI: 10.1016/j.jmb.2018.03.019
- Hu Y, et al. Multi-Knock. *Nature Plants*. 2023. DOI: 10.1038/s41477-023-01374-4
- Berman A, et al. Multi-targeted CRISPR libraries in tomato. *Nature Communications*. 2025. DOI: 10.1038/s41467-025-59280-6

---

## Citation and license

Citation metadata is available in:

- [CITATION.cff](CITATION.cff)
- [.zenodo.json](.zenodo.json)

A DOI should only be added after an actual archived release has been created and Zenodo has minted the DOI.

The original code in this repository is available under the **MIT License**.

This is research software and carries no guarantee of experimental performance.
