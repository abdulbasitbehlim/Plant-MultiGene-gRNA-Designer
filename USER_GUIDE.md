# Plant MultiGene gRNA Designer user guide

Version 1.4.0 helps you shortlist one SpCas9 guide for a group of plant genes. It can also suggest several exact-match guides if one shared guide is unavailable. A result is a computational candidate, not proof of editing.

## Run locally

Use Python 3.11 or newer. Open a terminal in the repository folder.

```bash
python -m venv .venv
```

On Windows PowerShell activate it with `.venv\Scripts\Activate.ps1`. On Linux or macOS use `source .venv/bin/activate`. Then run:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open http://localhost:8501. Manual FASTA works without a sequence-database connection after installation. Gene and accession retrieval require internet access. The optional Cas-OFFinder handoff requires a separate installation and an appropriate OpenCL device/runtime.

## Prepare inputs

Use distinct, reviewed genes from one organism. Prefer the assembly and cultivar you will study. The app does not establish homology from gene names. Recommended interactive size is 2 to 12 genes; 20 is the hard gene limit. Limit total target input and any reference panel to 500,000 bases each.

Choose Gene lookup for symbols and organism, Accession ID for known nucleotide or transcript identifiers, or Manual multi-FASTA for reviewed sequences. Ambiguous database lookup is rejected. Multiple accessions from one gene should not be counted as distinct genes. One selected transcript does not demonstrate coverage of every isoform.

The simplest FASTA format is one continuous region per gene:

```text
>GeneA
ACGTACGTACGTACGTACGAAGG
>GeneB
ACGTACGTACGTACGTACGATGG
```

For separate exon pieces, enable **Group separate exons using GeneID|SegmentID headers**:

```text
>GeneA|exon1
ACGTACGTACGTACGTACGAAGG
>GeneA|exon2
ATATATATATATATATATATATA
>GeneB|exon1
ACGTACGTACGTACGTACGATGG
```

These are synthetic teaching sequences. The app counts two genes in the second example and never joins the exon pieces. Keep every header unique. Unresolved bases are normalized to N; spacer windows containing them are skipped. An unresolved base at the first, degenerate position of NGG is allowed.

## Choose a search

Start with **Allow mismatch-aware consensus** switched off for an exact shared-guide search. If no guide is found, inspect PAM availability and sequence context. You can inspect the exact guide-set fallback before allowing mismatches.

With mismatch mode enabled, choose total mismatches, seed mismatches and minimum compatibility. The default values are 2, 1 and 55. They are software thresholds, not universally validated biological settings. Positions 13 to 20 of the spacer are the app's eight-base PAM-proximal seed. A stricter starting point is zero seed mismatches, with deliberate review before relaxing it.

Click **Design shared guides**. Both DNA strands are scanned for a 20-base spacer next to NGG. The algorithm retains feasible observed spacers and checks majority-consensus proposals. It does not search every possible synthetic 20-mer. Exact shared guides are ranked first; the remaining order favors gene coverage, number of exact targets, minimum compatibility, average compatibility, mismatch count and sequence quality, with deterministic ties.

## Read the result

Genes covered means distinct input gene labels with a selected compatible target. Exact genes means zero spacer mismatches. Minimum compatibility describes the least compatible selected target under the app's heuristic; it is not a percentage probability. Sequence quality is also a heuristic. GC content, poly-T and long homopolymers are review aids.

PASS means the implemented rules passed. REVIEW flags a candidate or input needing closer examination. FAIL means a hard design requirement was not met. NOT SCREENED refers specifically to specificity, independently of the core design status. Manual inputs and missing CDS/exon context can create review flags even for exact matches.

Inspect Per-gene target evidence. Coordinates are 1-based, inclusive positions within each input segment. They are not chromosome coordinates. Mismatch positions follow the guide orientation, including reverse-strand hits.

Changing sidebar settings does not alter a saved run. A warning appears and the previous settings remain attached to the results and exports until you submit a new design.

## If no single guide is found

A zero result applies to the current inputs, limits and searched proposals. It does not prove that no synthetic spacer could work. **Exact guide set fallback** proposes a set of exact-match guides and lists any genes left uncovered. This is a greedy set cover: it aims to cover the remaining genes quickly but does not guarantee the smallest or biologically best set. Inspect each member's sequence quality and off-targets.

**Validate a custom 20-nt guide against these genes** remains available after zero results. It first chooses a target satisfying the saved limits; if none exists, it shows the best proxy match to explain failure. Enter only the spacer, without the PAM or scaffold.

## Screen a supplied reference panel

Paste nonempty FASTA in Optional reference-panel near-match screen. Choose a radius and click Screen selected guide. The result includes total hits, displayed hits and whether the displayed list is capped. Exact copies of intended target sites are included; classify individual loci yourself. The panel is not an exclusion list, and a chromosome containing an intended gene is not entirely exempt.

Changing the panel text or radius invalidates the displayed screen for the current context. Empty panels cannot receive a pass. No hits means no matching NGG sites within the supplied DNA and substitution radius. It says nothing about DNA that was not supplied, alternative PAMs or bulges.

## Export and reproduce

| Export | Purpose |
|---|---|
| CSV | Ranked guides and validation summaries |
| Guide FASTA | Spacer sequences, without PAM or scaffold |
| Complete run JSON | Version, saved settings, exact input segments, provenance, target evidence, validation, panel metadata and fallback results |
| Input segments FASTA | Reuse with grouped GeneID\|SegmentID mode enabled |
| Fallback CSV | Exact guide-set members and covered genes |
| Displayed panel hits CSV | The displayed hit list; consult JSON for total count and truncation |
| Cas-OFFinder input | External genome-wide search preparation only |

For reproducibility, retain the JSON together with CSV/FASTA exports. Input FASTA includes all scanned segments; original segment names are retained in JSON. Reopen the app, paste the input FASTA, enable grouping, restore saved settings from JSON, and rerun. A JSON import button is not implemented.

## Prepare a genome specificity check

Open Prepare an external genome specificity search. Enter the local path to the exact plant reference FASTA file or directory, choose the mismatch radius, and decide whether to include NAG with NGG. Download the input file and run a separately installed Cas-OFFinder 2:

```bash
cas-offinder cas_offinder_input.txt C cas_offinder_hits.tsv
```

Replace the example genome path before running. This template uses substitutions only and does not request bulges. The app exports the returned shared guides, or the fallback members when no shared guides exist. Genome scanning, locus annotation, intended-target classification and interpretation of the output happen outside this app.

## Examples and scientific limits

Five paste-ready synthetic FASTA files are in examples/. Run `python benchmarks/run_audit_examples.py` for 13 synthetic checks. Run `python benchmarks/run_locked_benchmark.py` for the selected published PUP validation-construct fixture. Neither is a new plant experiment.

Read docs/RESEARCH_AUDIT.md for the original limitations, research links and complete verification results. The update does not add calibrated activity prediction, chromosome mapping, all-isoform validation, cultivar variation analysis, automatic genome screening or proof of experimental success.
