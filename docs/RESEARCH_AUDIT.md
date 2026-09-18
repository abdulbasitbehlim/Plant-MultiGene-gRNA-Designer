# Research audit and validation of version 1.4.0

Reviewed on 18 September 2026. Baseline: v1.3.1, commit `7cab4a80d8a7628d2c7ebcac39e27b9f13ba24f4`.

## Main conclusion

The tool is suitable for transparent candidate prioritization across a small, reviewed plant gene set. The update corrects input, target-selection and application-state failures. It does not establish experimental editing efficiency or genome-wide specificity. The scientific core is a heuristic shared-guide search, not a reimplementation of CRISPys.

## Research that informed the review

- **Prykhozhij et al. (2015), CRISPR MultiTargeter.** Searches common and unique targets in related sequences. Its annotated workflows check candidate sites against individual exons, and its paper directs users to downstream off-target tools. This supports keeping exon pieces separate and treating genome screening as a distinct step. DOI: https://doi.org/10.1371/journal.pone.0119372
- **Hyams et al. (2018), CRISPys.** Uses hierarchical clustering of potential targets and proposes guides for gene-family subgroups. Synthetic spacers and guide sets extend beyond a simple intersection of exact matches. The present app retains a smaller heuristic search and makes its incompleteness explicit. DOI: https://doi.org/10.1016/j.jmb.2018.03.019
- **Hu et al. (2023), Multi-Knock.** Demonstrates the biological value of targeting redundant Arabidopsis gene families. The repository's existing PUP7/PUP8/PUP21 fixture concerns a published follow-up validation construct; it is not a verified library-wide benchmark. DOI: https://doi.org/10.1038/s41477-023-01374-4
- **Berman et al. (2025), tomato multi-targeted libraries.** Combines CRISPys, CFD filtering, coding-region restrictions and genome-wide off-target analysis. Its experiments show lower cleavage efficiency as mismatch counts increase, including among guides passing CFD thresholds. Its NPF1.10/NPF1.11/NPF1.12 example targets three genes without mismatches. This motivates offering exact guide sets rather than making mismatch relaxation the only next step. DOI: https://doi.org/10.1038/s41467-025-59280-6

The paper summaries above are not a head-to-head accuracy comparison. No external CRISPys or CRISPR MultiTargeter run was performed during this update.

## Concrete findings and changes

| Limitation in the baseline | Change | Evidence |
|---|---|---|
| A nearest weighted-distance target could violate limits even when another site was feasible | Filter sites by every constraint before selecting; custom validation uses the same limits | Adversarial two-gene test |
| A feasible observed spacer could be replaced by the consensus proposal | Retain observed feasible spacers as well as accepted consensus proposals | Regression tests and PUP fixture |
| Equal-score ordering depended on set/input ordering | Deterministic gene, site and spacer tie breaks | Reordered inputs produce identical objects |
| Duplicate FASTA identifiers silently overwrote records | Reject duplicate and empty records | Parser rejection tests |
| Manual exons were either concatenated or counted as different genes | Optional GeneID\|SegmentID grouping with independent segment scans | Two exons form one gene without a junction target |
| Short annotated exons could trigger a fallback to the spliced CDS | Preserve short exon pieces; never replace known boundaries with a joined CDS | 20 nt plus 3 nt exon example |
| Compound CDS bounds could include intervening introns | Extract individual CDS parts and intersections with exons | Both strand annotations tested |
| Multi-CDS records selected the first CDS | Reject ambiguous multi-CDS accessions | Explicit error test |
| Ambiguous gene lookups selected the first NCBI Gene result | Reject multiple matches; request a unique accession | Mocked retrieval test |
| Repeated gene labels or aliases could inflate gene coverage | Reject duplicate resolved labels/records and mixed known organisms | Identity tests |
| A plant transcript suffix could be removed as if it were an ENS version | Preserve plant transcript IDs such as AT4G18197.1 | Requested URL regression |
| Changing sidebar controls could reinterpret an old run | Snapshot settings and sequence inputs; warn until design is rerun | Streamlit AppTest |
| Panel results were keyed only by spacer | Key by run, guide, panel text and mismatch radius | Streamlit rerun and key tests |
| Empty panel could produce a specificity pass | Reject empty panel input | Core and Streamlit tests |
| Displayed hit cap concealed the total | Report total hits, displayed count, truncation, radius and panel hash | Seven hits with a two-row cap |
| Custom guide validation was inaccessible after zero results | Keep custom validation and exports available | Streamlit zero-result workflow |
| No constructive alternative when one guide fails | Greedy exact guide-set fallback with uncovered genes | Complete and incomplete fallback tests |
| JSON lacked actual sequences and chosen thresholds | Complete run JSON plus independent-segment FASTA | Snapshot/export tests |
| External genome screening required hand-formatting | Cas-OFFinder 2 input export, NGG or NGG plus NAG | Format tests; external execution not claimed |
| Installed package omitted the accession module | Include accession_sources and run_state | Wheel contents inspected |

## Measured verification

The unchanged baseline passed 37 tests. The revised suite passes 69 tests, including three Streamlit workflow tests. Measured statement-and-branch coverage is 90.52% over the configured non-UI modules. app.py is excluded from that percentage and is exercised separately through Streamlit AppTest. Tests ran on Python 3.12.14; the CI configuration also requests 3.11 and 3.13, but those environments were not locally rerun.

`python benchmarks/run_audit_examples.py` passes 13 synthetic demonstrations and writes `benchmarks/results/audit_examples.json`. These are controlled software checks, not new biological experiments. `python benchmarks/run_locked_benchmark.py` regenerates the published PUP validation spacer at rank 1 from the existing locked context fixture: 3 PAM-compatible targets, 2 exact targets and one target with 1 mismatch. There are 2 generated candidates after retaining observed alternatives. This is one selected locus-context fixture, not a sensitivity or specificity estimate, an independent genome-wide benchmark, or evidence of library membership.

## Remaining limits

The compatibility proxy is not CFD and is not calibrated to a plant species. Its position penalties are a transparent heuristic. Search still starts from observed spacers that are feasible in all genes; an unobserved optimal spacer can therefore be missed. A zero-result search means no accepted candidate among those searched, not a mathematical proof that none exists.

The app cannot verify homology or whether a manually named input represents a distinct locus. It normally uses one representative transcript; it does not prove coverage of all isoforms, alleles or homeologs. Manual sequences and transcripts lacking boundaries need genomic review. Ensembl exons may include untranslated sequence. Coordinates are local to each scanned segment, not chromosome coordinates. Missing flanking sequence can omit legitimate targets near exon ends.

The local panel checks only NGG sites and substitutions in the supplied DNA. NAG sites, bulges, structural variants and unprovided genomic regions are outside its scope. The Cas-OFFinder file is a handoff for external execution, not a completed genome screen. Every intended and unintended locus must be classified individually; a whole chromosome must not be exempted as an intended target.

The interactive input and panel bounds are each 500,000 bases; at most 20 genes and 5,000,000 estimated mismatch comparisons are supported. The greedy exact guide set is not guaranteed to be the smallest set and does not optimize expression or off-targets. There is no plant-calibrated activity model, chromatin model, cleavage-outcome predictor, vector design, all-isoform mapping, or experimental validation in this update.

## Reproduce

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-fail-under=80
python benchmarks/run_audit_examples.py
python benchmarks/run_locked_benchmark.py
streamlit run app.py
```

Cas-OFFinder input syntax follows its official documentation: https://github.com/snugel/cas-offinder . The upstream README directs production users to Cas-OFFinder 2 rather than its experimental version 3. No external scanner is bundled or executed by this app.
