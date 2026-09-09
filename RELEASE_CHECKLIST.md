# Release / archive checklist

This repository is release-ready but **no GitHub release or Zenodo DOI has been created by this package build**.

1. Run `pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing` in a clean environment.
2. Run `docker compose up --build` and complete `UI_TEST_CHECKLIST.md` in a browser.
3. Confirm benchmark status and rerun any live/external comparison adapters.
4. Commit the exact benchmark result files used in the paper.
5. Create a version tag (for example `v1.3.1`) only after review.
6. Connect the public repository to Zenodo, enable GitHub release archiving, create the release, and copy the minted DOI into `CITATION.cff`, `.zenodo.json`, README, and JOSS `paper.md`.
7. Do not claim a DOI before Zenodo actually mints one.
