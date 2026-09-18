# Changelog

## 1.4.0

Correct target selection under mismatch constraints; retain feasible observed spacers; make rankings deterministic. Reject duplicate FASTA records, ambiguous NCBI gene results, multiple-CDS accessions, duplicate resolved records and mixed organisms. Add grouped exon FASTA, preserve short exons and compound CDS parts, and preserve plant transcript ID suffixes.

Save immutable run settings and actual sequence inputs. Scope panel evidence to run, guide, sequence text and radius. Reject empty panels and disclose hit truncation. Keep custom validation usable after no results. Add exact guide-set fallback, complete JSON and input FASTA exports, and Cas-OFFinder 2 input preparation. Include missing modules in the installable package.

Validation: 69 passing tests, 90.52% configured non-UI coverage, 13 synthetic demonstrations, and the existing locked PUP validation-construct regression. See docs/RESEARCH_AUDIT.md for scientific boundaries and research sources. No new editing-efficiency or whole-genome specificity claim.
