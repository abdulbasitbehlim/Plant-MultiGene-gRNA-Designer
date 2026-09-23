# ============================================================================
# APP
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Builds the Streamlit interface and connects sequence input, multi-gene guide design, validation and reporting.
#
# HOW TO READ THIS FILE:
# 1. Read the imports/constants first to see which tools and settings are used.
# 2. Read one top-level function or class at a time.
# 3. Follow the workflow from input sequence -> candidate guides -> validation -> output.
# 4. Scientific formulas, thresholds, validation decisions and public function names
#    are intentionally preserved while readability comments are added.
#
# MAIN TOP-LEVEL PARTS:
# - function: _table_value
# - function: key_value_table
# - function: clear_design_state
# ============================================================================

#!/usr/bin/env python3
"""Streamlit dashboard for plant multi-gene shared sgRNA design."""
from __future__ import annotations

import json
import re
from typing import Dict, List

import pandas as pd
import plotly.express as px
import streamlit as st

from plant_multiguide import (
    design_shared_guides,
    estimate_search_diagnostics,
    MAX_GENES, RECOMMENDED_MAX_GENES, MAX_PAIR_COMPARISONS,
    evaluate_candidate_guide,
    guide_summary_row,
    match_rows,
    screen_reference_panel_report,
    suggest_exact_guide_set,
)
from sequence_sources import fetch_gene, manual_records, parse_multifasta, validate_record_set
from accession_sources import fetch_accession
from validation import validate_shared_guide, validation_summary_row

from run_state import (APP_VERSION, create_run_snapshot, panel_result_key, export_input_fasta, export_run, cas_offinder_input)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _table_value
# ----------------------------------------------------------------------------
def _table_value(value):
    """Convert structured values into readable text for UI tables."""
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, dict):
        return " | ".join(
            f"{str(key).replace('_', ' ').title()}: {_table_value(item)}"
            for key, item in value.items()
        ) or "None"
    if isinstance(value, (list, tuple, set)):
        items = list(value)
        if isinstance(value, set):
            items = sorted(items, key=str)
        return ", ".join(_table_value(item) for item in items) or "None"
    return str(value)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: key_value_table
# ----------------------------------------------------------------------------
def key_value_table(data, field_label="Field", value_label="Value"):
    """Return a consistent two-column table instead of exposing raw JSON in the UI."""
    return pd.DataFrame(
        [
            {
                field_label: str(key).replace("_", " ").title(),
                value_label: _table_value(value),
            }
            for key, value in data.items()
        ]
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: clear_design_state
# ----------------------------------------------------------------------------
def clear_design_state():
    for key in list(st.session_state):
        if key.startswith("plant_") and key not in {"plant_custom_spacer"}:
            del st.session_state[key]

st.set_page_config(
    page_title="Plant MultiGene CRISPR Designer",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------- Theme ----------
dark_mode = st.sidebar.toggle("Dark mode", value=True)
P = (
    dict(app="#07130f", side="#0b1b14", panel="#10251b", panel2="#153022", text="#f1f8f4", muted="#a6b9ae", border="#2c4c3a", accent="#4ade80", accent2="#22d3ee", warn="#fbbf24", plot="plotly_dark")
    if dark_mode
    else dict(app="#f5faf7", side="#edf7f0", panel="#ffffff", panel2="#f0f7f2", text="#17251d", muted="#66756c", border="#d3e4d8", accent="#15803d", accent2="#0891b2", warn="#a16207", plot="plotly_white")
)
st.markdown(f"""
<style>
:root{{--app:{P['app']};--side:{P['side']};--panel:{P['panel']};--panel2:{P['panel2']};--text:{P['text']};--muted:{P['muted']};--border:{P['border']};--accent:{P['accent']};--accent2:{P['accent2']};--warn:{P['warn']};}}
html,body,[data-testid="stAppViewContainer"],.stApp{{background:var(--app)!important;color:var(--text)!important}}
[data-testid="stSidebar"]{{background:var(--side)!important;border-right:1px solid var(--border)!important}}
[data-testid="stHeader"]{{background:transparent!important}}
.block-container{{max-width:1500px;padding-top:1.2rem;padding-bottom:4rem}}
.stApp h1,.stApp h2,.stApp h3,.stApp h4,.stApp p,.stApp li,.stApp label{{color:var(--text)!important}}
[data-testid="stCaptionContainer"] p{{color:var(--muted)!important}}
[data-baseweb="input"] input,[data-baseweb="textarea"] textarea,.stApp input,.stApp textarea{{color:var(--text)!important;-webkit-text-fill-color:var(--text)!important}}
[data-baseweb="select"] span,[data-baseweb="select"] div,[data-baseweb="select"] svg{{color:var(--text)!important;fill:var(--text)!important}}
.stApp input::placeholder,.stApp textarea::placeholder{{color:var(--muted)!important;opacity:.85!important}}
.hero{{padding:2rem 2.2rem;border:1px solid var(--border);border-radius:22px;background:linear-gradient(130deg,var(--panel),var(--panel2));margin-bottom:1.2rem}}
.hero .eyebrow{{font-size:.75rem;letter-spacing:.15em;text-transform:uppercase;color:var(--accent)!important;font-weight:800}}
.hero h1{{font-size:clamp(2rem,4.5vw,4rem);line-height:1;margin:.6rem 0 .8rem;letter-spacing:-.04em}}
.hero p{{color:var(--muted)!important;max-width:940px;line-height:1.65}}
.pills{{display:flex;gap:.5rem;flex-wrap:wrap;margin-top:1rem}} .pills span{{border:1px solid var(--border);background:var(--panel2);border-radius:999px;padding:.35rem .65rem;font-size:.78rem}}
[data-testid="stForm"],[data-testid="stExpander"]{{background:var(--panel)!important;border:1px solid var(--border)!important;border-radius:16px!important}}
div[data-testid="stMetric"]{{background:var(--panel);border:1px solid var(--border);padding:.9rem 1rem;border-radius:14px}}
.stButton>button,.stDownloadButton>button{{border-radius:10px;min-height:2.7rem;font-weight:750;border:1px solid var(--border)}}
.stButton>button[kind="primary"]{{background:linear-gradient(90deg,var(--accent),var(--accent2));color:#04110a!important;border:0}}
[data-testid="stDataFrame"]{{border:1px solid var(--border);border-radius:12px;overflow:hidden}}
.note{{background:var(--panel2);border-left:3px solid var(--accent);padding:.85rem 1rem;border-radius:10px;margin:.5rem 0}}
.warn{{background:var(--panel2);border-left:3px solid var(--warn);padding:.85rem 1rem;border-radius:10px;margin:.5rem 0}}
</style>
""", unsafe_allow_html=True)

st.sidebar.markdown(f"### Plant MultiGene Designer\n**Shared SpCas9 guides · v{APP_VERSION}**")
st.sidebar.subheader("Design settings")
include_mismatch = st.sidebar.toggle("Allow mismatch-aware consensus", value=True)
max_mm = st.sidebar.slider("Max mismatches per target gene", 0, 3, 2, disabled=not include_mismatch)
max_seed_mm = st.sidebar.slider("Max PAM-proximal seed mismatches", 0, 2, 1, disabled=not include_mismatch)
min_compat = st.sidebar.slider("Minimum compatibility proxy", 0, 100, 55, 5, disabled=not include_mismatch)
max_results = st.sidebar.slider("Maximum returned guides", 10, 200, 50, 10)
with st.sidebar.expander("Validation rules"):
    st.markdown("- Hard checks: 20 nt, A/C/G/T only, all genes covered, NGG at every target, mismatch/seed thresholds.\n- Review checks: preferred GC 40–60%, poly-T, long homopolymers.\n- Local specificity: optional supplied FASTA panel; not a whole-genome claim.")
st.sidebar.info(f"Research-use prioritization only. Recommended interactive range: 2–{RECOMMENDED_MAX_GENES} genes; hard cap: {MAX_GENES}. Mismatch-aware search also stops above {MAX_PAIR_COMPARISONS:,} estimated pair comparisons. Whole-genome off-target review remains essential.")

st.markdown(f"""
<div class="hero">
  <div class="eyebrow">Gene-family targeting</div>
  <h1>Plant <span style="color:{P['accent']}">MultiGene</span> gRNA Designer</h1>
  <p>Enter two, three, or more homologous plant genes by gene lookup, database accession ID, or reviewed multi-FASTA. The app first searches for an exact 20-nt SpCas9 spacer with an NGG PAM in every gene. When enabled, it can additionally rank mismatch-aware consensus guides that retain a PAM-compatible target in every requested gene.</p>
  <div class="pills"><span>20 nt + NGG</span><span>Both strands</span><span>Accession ID</span><span>Exact shared guides</span><span>Consensus mode</span><span>Guide validation</span><span>CSV / FASTA / JSON</span></div>
</div>
""", unsafe_allow_html=True)

input_mode = st.radio("Input mode", ["Gene lookup", "Accession ID", "Manual multi-FASTA"], horizontal=True)
records = None
submitted = False

if input_mode == "Gene lookup":
    with st.form("lookup_form"):
        c1, c2 = st.columns([2, 1])
        with c1:
            genes_raw = st.text_area("Gene symbols / IDs", placeholder="Example:\nPUP7\nPUP8\nPUP21", height=140, help="Separate genes with commas, spaces, or new lines.")
        with c2:
            organism = st.text_input("Plant organism", value="Arabidopsis thaliana")
            source = st.selectbox("Sequence source", ["NCBI RefSeq", "Ensembl REST"])
        submitted = st.form_submit_button("Design shared guides", type="primary", use_container_width=True)
    if submitted:
        clear_design_state()
        genes: List[str] = [g for g in re.split(r"[\s,;]+", genes_raw.strip()) if g]
        genes = list(dict.fromkeys(genes))
        if len(genes) < 2:
            st.error("Enter at least two distinct genes.")
            st.stop()
        if len(genes) > MAX_GENES:
            st.error(f"This interactive tool supports at most {MAX_GENES} genes per run. For larger libraries use a dedicated scalable workflow.")
            st.stop()
        if len(genes) > RECOMMENDED_MAX_GENES:
            st.warning(f"You entered {len(genes)} genes. The recommended interactive range is 2–{RECOMMENDED_MAX_GENES}; the workload guard may stop permissive searches.")
        records = {}
        progress = st.progress(0, text="Retrieving gene sequences...")
        errors = []
        for i, gene in enumerate(genes, start=1):
            try:
                records[gene] = fetch_gene(gene, organism, source)
            except Exception as exc:
                errors.append(f"{gene}: {exc}")
            progress.progress(i / len(genes), text=f"Retrieved {i}/{len(genes)} gene records")
        progress.empty()
        if errors:
            st.error("One or more genes could not be resolved:\n\n" + "\n".join(f"- {e}" for e in errors))
            st.info("Use Accession ID or Manual multi-FASTA for genes that are absent or ambiguously annotated in the selected database.")
            st.stop()
elif input_mode == "Accession ID":
    with st.form("accession_form"):
        c1, c2 = st.columns([2, 1])
        with c1:
            accessions_raw = st.text_area(
                "Gene accession IDs",
                placeholder="Example:\nNM_...\nXM_...\n\nor Ensembl stable gene/transcript IDs",
                height=140,
                help="Enter at least two accession IDs. Separate them with commas, spaces, semicolons, or new lines.",
            )
        with c2:
            source = st.selectbox("Accession source", ["NCBI RefSeq / Nucleotide", "Ensembl REST"])
        submitted = st.form_submit_button("Fetch accessions and design shared guides", type="primary", use_container_width=True)
    if submitted:
        clear_design_state()
        accession_ids: List[str] = [x for x in re.split(r"[\s,;]+", accessions_raw.strip()) if x]
        accession_ids = list(dict.fromkeys(accession_ids))
        if len(accession_ids) < 2:
            st.error("Enter at least two distinct accession IDs.")
            st.stop()
        if len(accession_ids) > MAX_GENES:
            st.error(f"This interactive tool supports at most {MAX_GENES} accession IDs per run.")
            st.stop()
        if len(accession_ids) > RECOMMENDED_MAX_GENES:
            st.warning(f"You entered {len(accession_ids)} accessions. The recommended interactive range is 2–{RECOMMENDED_MAX_GENES} genes.")
        records = {}
        errors = []
        progress = st.progress(0, text="Retrieving accession records...")
        for i, accession in enumerate(accession_ids, start=1):
            try:
                rec = fetch_accession(accession, source)
                key = rec.gene or accession
                if key in records:
                    key = f"{key}|{accession}"
                records[key] = rec
            except Exception as exc:
                errors.append(f"{accession}: {exc}")
            progress.progress(i / len(accession_ids), text=f"Retrieved {i}/{len(accession_ids)} accession records")
        progress.empty()
        if errors:
            st.error("One or more accessions could not be resolved:\n\n" + "\n".join(f"- {e}" for e in errors))
            st.stop()
        if len(records) < 2:
            st.error("Fewer than two usable gene records were retrieved.")
            st.stop()
else:
    with st.form("fasta_form"):
        organism = st.text_input("Plant organism / cultivar label", value="Plant species")
        group_segments = st.checkbox("Group separate exons using GeneID|SegmentID headers", value=False)
        raw_fasta = st.text_area(
            "One FASTA record per gene",
            height=260,
            placeholder=">GeneA\nACGT...\n>GeneB\nACGT...\n>GeneC\nACGT...",
            help="Use genomic exon/CDS segments when possible. A single record per gene is treated as one continuous segment.",
        )
        submitted = st.form_submit_button("Design shared guides", type="primary", use_container_width=True)
    if submitted:
        clear_design_state()
        try:
            records = manual_records(raw_fasta, organism=organism, group_segments=group_segments)
        except Exception as exc:
            st.error(str(exc))
            st.stop()
        if len(records) < 2:
            st.error("Provide at least two FASTA records.")
            st.stop()
        if len(records) > MAX_GENES:
            st.error(f"This interactive tool supports at most {MAX_GENES} FASTA records/genes per run.")
            st.stop()
        if len(records) > RECOMMENDED_MAX_GENES:
            st.warning(f"{len(records)} genes exceeds the recommended interactive range of {RECOMMENDED_MAX_GENES}; permissive mismatch-aware searches may hit the workload guard.")

current_settings = {"include_mismatch_aware": include_mismatch, "max_mismatches_per_gene": max_mm,
                    "max_seed_mismatches_per_gene": max_seed_mm, "min_compatibility": min_compat,
                    "max_results": max_results}

if submitted and records:
    gene_segments = {gene: rec.segments for gene, rec in records.items()}
    try:
        validate_record_set(records)
        sites_by_gene, guides = design_shared_guides(
            gene_segments,
            include_mismatch_aware=include_mismatch,
            max_mismatches_per_gene=max_mm,
            max_seed_mismatches_per_gene=max_seed_mm,
            min_compatibility=min_compat,
            max_results=max_results,
        )
    except Exception as exc:
        st.error(f"Design failed: {exc}")
        st.stop()
    diagnostics = estimate_search_diagnostics(sites_by_gene)
    st.session_state["plant_snapshot"] = create_run_snapshot(records, current_settings)
    st.session_state["plant_records"] = records
    st.session_state["plant_sites"] = sites_by_gene
    st.session_state["plant_guides"] = guides
    st.session_state["plant_search_diagnostics"] = diagnostics

if "plant_guides" in st.session_state:
    records = st.session_state["plant_records"]
    sites_by_gene = st.session_state["plant_sites"]
    guides = st.session_state["plant_guides"]
    snapshot = st.session_state["plant_snapshot"]
    settings = snapshot["settings"]
    diagnostics = st.session_state["plant_search_diagnostics"]
    if settings != current_settings:
        st.warning("Settings have changed. Displayed results and exports still use the saved run settings below. Submit Design shared guides again to apply the new settings.")
    st.caption(f"Saved run {snapshot['run_id'][:12]} · {snapshot['created_at_utc']}")
    st.markdown("#### Saved design settings")
    st.dataframe(
        key_value_table(settings, "Setting", "Saved value"),
        hide_index=True,
        use_container_width=True,
    )
    validation_kwargs = dict(expected_genes=len(records), expected_gene_ids=list(records),
        max_mismatches_per_gene=settings["max_mismatches_per_gene"],
        max_seed_mismatches_per_gene=settings["max_seed_mismatches_per_gene"],
        min_compatibility=settings["min_compatibility"],
        input_warnings=[w for rec in records.values() for w in rec.warnings])

    st.subheader("Design summary")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Genes", len(records))
    c2.metric("PAM-compatible sites", sum(map(len, sites_by_gene.values())))
    c3.metric("Returned shared guides", len(guides))
    c4.metric("Exact shared", sum(g.design_type == "Exact shared" for g in guides))
    st.caption("Returned guides may be limited by Maximum returned guides. Coordinates are 1-based within each input segment, not chromosome coordinates.")

    with st.expander("Sequence provenance, reproducibility and warnings"):
        for gene, rec in records.items():
            st.markdown(f"**{gene}** — {rec.source} · `{rec.accession}` · {rec.total_bp:,} bp in {len(rec.segments)} segments")
            st.caption(f"sequence SHA-256: {rec.sequence_sha256}")
            st.dataframe(
                key_value_table(rec.provenance_dict(), "Provenance field", "Value"),
                hide_index=True,
                use_container_width=True,
            )
            for warning in rec.warnings:
                st.warning(warning)
    if settings["include_mismatch_aware"]:
        with st.expander("Mismatch-aware search strategy and computational workload"):
            st.write(diagnostics.strategy)
            st.write(f"Upper bound: {diagnostics.upper_pair_comparisons:,} candidate-to-site comparisons.")

    with st.expander("Exact guide set fallback", expanded=not bool(guides)):
        st.caption("A separate option using several guides. This greedy exact-match set cover is not guaranteed to use the fewest guides and does not predict editing success.")
        fallback_limit = st.slider("Maximum guides in fallback set", 1, 20, 5)
        fallback = suggest_exact_guide_set(sites_by_gene, fallback_limit)
        rows = [{**guide_summary_row(g), "Target genes": ", ".join(m.gene for m in g.matches)} for g in fallback["guides"]]
        if rows:
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
            st.download_button("Download fallback CSV", pd.DataFrame(rows).to_csv(index=False), "plant_exact_guide_set.csv", "text/csv")
        if fallback["uncovered_genes"]:
            st.warning("Uncovered genes: " + ", ".join(fallback["uncovered_genes"]))
        else:
            st.info("Every input gene has an exact match in this guide set. Specificity and sequence-quality review are still required for each member.")

    # Panel evidence is scoped to the current run, spacer, FASTA text and radius.
    panel_reports = {}
    selected = None
    if not guides:
        st.warning("No single guide satisfied the current all-gene constraints among the candidates searched. The synthetic spacer search is not exhaustive.")
    else:
        labels = [f"{i+1}. {g.spacer} · {g.design_type}" for i, g in enumerate(guides)]
        selected_label = st.selectbox("Inspect one guide", labels)
        selected = guides[labels.index(selected_label)]
        with st.expander("Optional reference-panel near-match screen"):
            st.caption("NGG sites and substitutions only; up to 500,000 input bases. Intended targets are also listed and must be classified by their individual loci. This is not a whole-genome clearance.")
            panel_raw = st.text_area("Reference panel FASTA", key="panel_fasta", height=160)
            panel_mm = st.slider("Panel mismatch radius", 0, 4, 3, key="panel_mm")
            key = panel_result_key(snapshot["run_id"], selected.spacer, panel_raw, panel_mm)
            if st.button("Screen selected guide", key="screen_panel"):
                try:
                    report = screen_reference_panel_report(selected.spacer, parse_multifasta(panel_raw), panel_mm)
                    st.session_state[key] = report
                except Exception as exc:
                    st.error(str(exc))
            for candidate in guides:
                candidate_key = panel_result_key(snapshot["run_id"], candidate.spacer, panel_raw, panel_mm)
                if candidate_key in st.session_state:
                    panel_reports[candidate.spacer] = st.session_state[candidate_key]
            report = panel_reports.get(selected.spacer)
            if report is not None:
                st.write(f"{report.total_hits} total hits; {len(report.hits)} displayed; {report.scanned_sites} NGG sites examined.")
                if report.truncated:
                    st.warning("The displayed hit list is truncated. Total hit count includes every hit in this bounded scan.")
                if report.hits:
                    st.dataframe(pd.DataFrame([h.__dict__ for h in report.hits]), hide_index=True, use_container_width=True)
                    st.download_button("Download displayed panel hits", pd.DataFrame([h.__dict__ for h in report.hits]).to_csv(index=False), "plant_panel_hits.csv", "text/csv")
                else:
                    st.info("No matches found within this panel and radius. Whole-genome specificity has not been established.")
            else:
                st.caption("No saved screen matches this run, guide, panel text and radius.")

    reports = [validate_shared_guide(g, panel_hits=panel_reports.get(g.spacer), **validation_kwargs) for g in guides]
    summary_df = pd.DataFrame([{**guide_summary_row(g), **validation_summary_row(r)} for g, r in zip(guides, reports)])
    if guides:
        st.markdown("#### Ranked guides with validation")
        st.dataframe(summary_df, use_container_width=True, hide_index=True)
        st.caption("Sequence quality and compatibility are uncalibrated ranking aids. PASS means software rules passed, not that editing was experimentally demonstrated.")
        selected_report = reports[guides.index(selected)]
        st.markdown("#### Validation")
        st.write(f"{selected_report.status} · Specificity: {selected_report.specificity_status}")
        st.dataframe(pd.DataFrame([c.__dict__ for c in selected_report.checks]), hide_index=True, use_container_width=True)
        st.markdown("#### Per-gene target evidence")
        st.dataframe(pd.DataFrame(match_rows(selected)), hide_index=True, use_container_width=True)
        fig = px.scatter(summary_df, x="GC%", y="Minimum compatibility", hover_name="Spacer (20 nt)", color="Design type")
        fig.update_layout(template=P["plot"], height=320)
        st.plotly_chart(fig, use_container_width=True)

    with st.expander("Validate a custom 20-nt guide against these genes"):
        st.caption("Selects a site satisfying the saved run limits when one exists. Otherwise shows the best proxy match to explain failure. Available even when the search returns no guides.")
        custom_spacer = st.text_input("Custom spacer (20 nt, no PAM)", key="plant_custom_spacer", max_chars=20).strip().upper()
        if st.button("Validate custom guide", key="plant_custom_validate"):
            try:
                custom = evaluate_candidate_guide(custom_spacer, sites_by_gene,
                    settings["max_mismatches_per_gene"], settings["max_seed_mismatches_per_gene"], settings["min_compatibility"])
                custom_report = validate_shared_guide(custom, **validation_kwargs)
                st.metric("Custom guide validation", custom_report.status)
                st.dataframe(pd.DataFrame([c.__dict__ for c in custom_report.checks]), hide_index=True)
                st.dataframe(pd.DataFrame(match_rows(custom)), hide_index=True)
            except Exception as exc:
                st.error(str(exc))

    st.subheader("Exports")
    if guides:
        st.download_button("Download CSV", summary_df.to_csv(index=False), "plant_shared_guides.csv", "text/csv")
        guide_fasta = "\n".join(f">shared_guide_{i+1}\n{g.spacer}" for i, g in enumerate(guides)) + "\n"
        st.download_button("Download guide FASTA", guide_fasta, "plant_shared_guides.fasta", "text/plain")
    payload = export_run(snapshot, guides, reports, panel_reports, fallback)
    payload["search"] = diagnostics.__dict__
    st.download_button("Download complete run JSON", json.dumps(payload, indent=2), "plant_design_run.json", "application/json")
    st.download_button("Download input segments FASTA", export_input_fasta(records), "plant_input_segments.fasta", "text/plain")
    st.caption("The JSON includes the actual input segments, saved settings, per-gene evidence, validation and available panel results. To rerun the exported input FASTA, enable grouped GeneID|SegmentID mode.")

    with st.expander("Prepare an external genome specificity search"):
        st.caption("Exports Cas-OFFinder 2 input without running it. Use your exact plant assembly and inspect intended versus unwanted loci. Bulges and cultivar variation require additional analysis.")
        genome_path = st.text_input("Local genome FASTA file or directory", value="/path/to/plant_genome_fasta")
        external_mm = st.slider("External mismatch radius", 0, 6, 3)
        include_nag = st.checkbox("Include NAG alongside NGG for specificity review", value=True)
        external_guides = guides or fallback["guides"]
        if external_guides:
            try:
                external_text = cas_offinder_input([g.spacer for g in external_guides], genome_path, external_mm, include_nag)
                st.download_button("Download Cas-OFFinder input", external_text, "cas_offinder_input.txt", "text/plain")
                st.code("cas-offinder cas_offinder_input.txt C cas_offinder_hits.tsv", language="bash")
            except ValueError as exc:
                st.error(str(exc))
    st.markdown("**Scientific boundary:** verify genomic continuity, coding context, cultivar/alleles and genome-wide specificity. One representative transcript does not establish coverage of every isoform. No editing efficiency, frameshift outcome or phenotype is guaranteed.")
