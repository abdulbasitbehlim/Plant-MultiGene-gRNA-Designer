# Streamlit UI release checklist

This checklist exists because the artifact build environment cannot install Streamlit from PyPI. Complete it on a machine with the project requirements installed.

## Automated

```bash
python -m venv .venv
# activate the environment
pip install -r requirements.txt
pip install -r requirements-dev.txt
pytest -q
python -m streamlit run app.py
```

## Browser smoke test

- App loads without traceback.
- Dark/light mode controls remain readable.
- Gene lookup form submits and errors are human-readable when a lookup fails.
- Manual FASTA form accepts valid FASTA and reports ambiguity warnings when IUPAC ambiguity is present.
- Provenance section shows source/version, assembly/genomic record, annotation release/date, retrieval UTC and SHA-256.
- Ranked table renders without clipped columns.
- PASS / REVIEW / FAIL validation checklist opens and matches the selected guide.
- Custom/existing-guide validator accepts a valid 20-nt spacer and rejects malformed input.
- CSV, FASTA and JSON downloads open successfully.
- JSON contains provenance fields and the exact sequence SHA-256.
- No result is described as experimentally validated solely because software validation says PASS.

## Project-specific review

For Plant MultiGene, also verify:
- the search-strategy/workload expander appears;
- >20 entered genes are blocked;
- a workload-guard error is understandable;
- no-result behavior does not force a shared guide.

For OpenCRISPR-1, also verify:
- the 2026 conflicting-evidence warning is visible;
- the local specificity screen is labelled MIT/Hsu **legacy baseline**;
- no UI text calls MIT/Hsu state-of-the-art or calls compatibility an efficacy guarantee.

Record browser, OS, Streamlit version, test date, and any screenshots/issues in your project review notes.
