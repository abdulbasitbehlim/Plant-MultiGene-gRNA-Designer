# Contributing

Contributions are welcome after the project is made public. Please keep scientific behavior reproducible and reviewable.

1. Open an issue describing the bug, scientific question or proposed feature before a large change.
2. Create a focused branch and add or update tests for changed scientific behavior.
3. Run `pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing` and complete `UI_TEST_CHECKLIST.md` for UI changes.
4. Do not relabel heuristics as published scores. New scoring backends must cite their method, state required context and include reference-value tests.
5. Preserve accession/version, assembly/release and SHA-256 provenance in any new sequence source.
6. Do not commit patient-identifiable sequence data, secrets, large genome indexes or proprietary third-party assets.

For scientific changes, describe the biological assumption, source/reference and expected impact in the pull request.
