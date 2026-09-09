# GitHub Deployment Guide — Plant MultiGene gRNA Designer

## Recommended repository

**Repository name:** `Plant-MultiGene-gRNA-Designer`

**Description:** Research-oriented Python/Streamlit tool for finding and validating one SpCas9 guide that can target multiple homologous plant genes, with exact and mismatch-aware modes, provenance, validation, tests and reproducible benchmarks.

## Before making the repository public

1. Run `python -m pip install -r requirements-dev.txt` in a clean Python 3.11–3.13 environment.
2. Run `pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-fail-under=80`.
3. Run `streamlit run app.py` and complete `UI_TEST_CHECKLIST.md`.
4. Test at least one manual multi-FASTA run and one live NCBI/Ensembl lookup.
5. Verify CSV/FASTA/JSON downloads and the custom-guide validator.
6. Do not describe the PUP benchmark as a verified Multi-Knock library-guide rediscovery; use the wording in `benchmarks/BENCHMARK_EVIDENCE_RESOLUTION.md`.

## Upload to GitHub

Create an empty repository and upload **the contents of this folder at repository root**. Do not create an extra nested project folder inside the repository.

The included `.github/workflows/ci.yml` runs the test/coverage matrix on pushes and pull requests. No GitHub secret is required for the normal test suite.

## Streamlit Community Cloud

After the repository is public, create an app from the repository and set the entrypoint to:

`app.py`

If you want to identify your NCBI requests with your own contact email, add `NCBI_EMAIL` as an app secret/environment value where supported. Never commit private credentials.

## Docker

```bash
docker compose up --build
```

Then open `http://localhost:8501`.

## First release

For the first public research-preview release, a conservative tag is `v0.1.0`. Keep the software's internal scientific artifact version (v1.3.1) documented in the changelog/release notes if desired. Do not mint or claim a Zenodo DOI until that exact public release has actually been archived.

## Scientific boundary

A PASS is a computational design/validation status under the implemented assumptions; it is not proof of editing efficiency, whole-genome specificity or wet-lab success.
