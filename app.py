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
    screen_reference_panel,
)
from sequence_sources import fetch_gene, manual_records, parse_multifasta
from accession_sources import fetch_accession
from validation import validate_shared_guide, validation_summary_row

APP_VERSION = "1.3.1"

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
        raw_fasta = st.text_area(
            "One FASTA record per gene",
            height=260,
            placeholder=">GeneA\nACGT...\n>GeneB\nACGT...\n>GeneC\nACGT...",
            help="Use genomic exon/CDS segments when possible. A single record per gene is treated as one continuous segment.",
        )
        submitted = st.form_submit_button("Design shared guides", type="primary", use_container_width=True)
    if submitted:
        try:
            records = manual_records(raw_fasta, organism=organism)
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

if submitted and records:
    gene_segments = {gene: rec.segments for gene, rec in records.items()}
    try:
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
    st.session_state["plant_records"] = records
    st.session_state["plant_sites"] = sites_by_gene
    st.session_state["plant_guides"] = guides
    st.session_state["plant_search_diagnostics"] = diagnostics

if "plant_guides" in st.session_state:
    records = st.session_state["plant_records"]
    sites_by_gene = st.session_state["plant_sites"]
    guides = st.session_state["plant_guides"]
    diagnostics = st.session_state.get("plant_search_diagnostics", estimate_search_diagnostics(sites_by_gene))

    st.subheader("Design summary")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Genes", len(records))
    m2.metric("PAM-compatible sites", sum(len(v) for v in sites_by_gene.values()))
    m3.metric("Shared-guide designs", len(guides))
    m4.metric("Exact shared", sum(g.design_type == "Exact shared" for g in guides))

    with st.expander("Sequence provenance, reproducibility and warnings"):
        for gene, rec in records.items():
            st.markdown(f"**{gene}** — {rec.source} · `{rec.accession}` · {rec.total_bp:,} bp across {len(rec.segments)} scanned segment(s)")
            st.caption(rec.description)
            st.code(
                f"record/version: {rec.source_record_version}\n"
                f"assembly/genomic record: {rec.assembly}\n"
                f"annotation release: {rec.annotation_release}\n"
                f"retrieved UTC: {rec.retrieved_at_utc}\n"
                f"sequence SHA-256: {rec.sequence_sha256}\n"
                f"ambiguous bases: {rec.ambiguity_count} ({', '.join(rec.ambiguity_codes) or 'none'})",
                language=None,
            )
            for warning in rec.warnings:
                st.warning(warning)

    if include_mismatch:
        with st.expander("Mismatch-aware search strategy and computational workload"):
            st.markdown(f"**Strategy:** {diagnostics.strategy}")
            st.markdown(
                f"**This run:** {diagnostics.genes} genes · {diagnostics.total_pam_sites:,} PAM-compatible sites · "
                f"{diagnostics.unique_seed_spacers:,} distinct observed seeds · up to "
                f"{diagnostics.upper_pair_comparisons:,} pair comparisons in the two nearest-neighbour passes."
            )
            st.caption(f"Recommended interactive range: 2–{RECOMMENDED_MAX_GENES} genes. Hard limits: {MAX_GENES} genes and {MAX_PAIR_COMPARISONS:,} estimated pair comparisons.")

    if not guides:
        st.warning("No single guide satisfied the current all-gene constraints. This is a scientifically valid result: the genes may not share a suitable PAM-compatible target under these settings.")
        st.markdown("Try reviewed genomic/exonic FASTA sequences, verify that the genes are homologous, or relax mismatch-aware thresholds carefully. For unrelated genes, use multiple sgRNAs rather than forcing one shared guide.")
        st.stop()

    validation_reports = [
        validate_shared_guide(
            g, expected_genes=len(records),
            max_mismatches_per_gene=max_mm,
            max_seed_mismatches_per_gene=max_seed_mm,
            min_compatibility=min_compat,
        )
        for g in guides
    ]
    summary_df = pd.DataFrame([
        {**guide_summary_row(g), **validation_summary_row(r)}
        for g, r in zip(guides, validation_reports)
    ])
    st.markdown("#### Ranked guides with validation")
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    c1, c2 = st.columns(2)
    with c1:
        fig = px.scatter(summary_df, x="GC%", y="Minimum compatibility", size="Sequence quality", hover_name="Spacer (20 nt)", symbol="Design type", title="Guide quality landscape")
        fig.update_layout(template=P["plot"], height=390, margin=dict(l=20, r=20, t=50, b=20))
        st.plotly_chart(fig, use_container_width=True)
    with c2:
        top = summary_df.head(20).copy()
        fig2 = px.bar(top, x="Spacer (20 nt)", y="Average compatibility", pattern_shape="Design type", title="Top shared-guide compatibility")
        fig2.update_layout(template=P["plot"], height=390, margin=dict(l=20, r=20, t=50, b=20), xaxis_tickangle=-45)
        st.plotly_chart(fig2, use_container_width=True)

    labels = [f"{i+1}. {g.spacer} · {g.design_type} · min {g.minimum_compatibility:.0f}" for i, g in enumerate(guides)]
    selected_label = st.selectbox("Inspect one guide", labels)
    selected = guides[labels.index(selected_label)]

    st.markdown(f"### `{selected.spacer}`")
    d1, d2, d3, d4 = st.columns(4)
    d1.metric("Design type", selected.design_type)
    d2.metric("Genes covered", selected.genes_covered)
    d3.metric("Worst mismatches", selected.worst_mismatch_count)
    d4.metric("GC", f"{selected.gc_percent:.1f}%")

    panel_state_key = f"plant_panel_hits::{selected.spacer}"
    selected_panel_hits = st.session_state.get(panel_state_key)
    selected_validation = validate_shared_guide(
        selected, expected_genes=len(records),
        max_mismatches_per_gene=max_mm,
        max_seed_mismatches_per_gene=max_seed_mm,
        min_compatibility=min_compat,
        panel_hits=selected_panel_hits,
    )
    st.markdown("#### Validation")
    v1, v2, v3, v4 = st.columns(4)
    v1.metric("Overall", selected_validation.status)
    v2.metric("PASS checks", selected_validation.pass_count)
    v3.metric("REVIEW checks", selected_validation.review_count)
    v4.metric("FAIL checks", selected_validation.fail_count)
    validation_df = pd.DataFrame([c.__dict__ for c in selected_validation.checks])
    st.dataframe(validation_df, use_container_width=True, hide_index=True)
    if selected_validation.status == "PASS":
        st.success("Core guide validation passed. Specificity status: " + selected_validation.specificity_status)
    elif selected_validation.status == "REVIEW":
        st.warning("Guide passes hard design requirements but has one or more review flags. Specificity status: " + selected_validation.specificity_status)
    else:
        st.error("Guide failed one or more hard validation requirements. Review the failed checks before prioritizing it.")

    st.markdown("#### Per-gene target evidence")
    st.dataframe(pd.DataFrame(match_rows(selected)), use_container_width=True, hide_index=True)
    for note in selected.notes:
        st.caption(note)

    with st.expander("Optional reference-panel near-match screen"):
        st.caption("This scans only the FASTA you provide. It is useful for a paralog/off-target panel, but it is not a whole-genome specificity analysis.")
        panel_raw = st.text_area("Reference panel FASTA", key="panel_fasta", height=160)
        panel_mm = st.slider("Panel mismatch radius", 0, 4, 3, key="panel_mm")
        if st.button("Screen selected guide", key="screen_panel"):
            try:
                panel = parse_multifasta(panel_raw)
                hits = screen_reference_panel(selected.spacer, panel, panel_mm)
                st.session_state[panel_state_key] = hits
                hit_df = pd.DataFrame([h.__dict__ for h in hits])
                if hit_df.empty:
                    st.success("No PAM-compatible hits were found within the selected mismatch radius in this supplied panel. Validation will record PASS for this supplied panel.")
                else:
                    st.warning("Near matches were found. Validation marks the panel as REVIEW because intended multi-gene targets must be distinguished from unwanted sites.")
                    st.dataframe(hit_df, use_container_width=True, hide_index=True)
            except Exception as exc:
                st.error(str(exc))

    with st.expander("Validate a custom 20-nt guide against these genes"):
        st.caption("This reuses the already scanned PAM-compatible sites for the loaded genes and reports the closest target in every gene. It is useful for checking a guide proposed by another tool or paper.")
        custom_spacer = st.text_input("Custom spacer (20 nt, no PAM)", key="plant_custom_spacer", max_chars=20).strip().upper()
        if st.button("Validate custom guide", key="plant_custom_validate"):
            try:
                custom_guide = evaluate_candidate_guide(custom_spacer, sites_by_gene)
                custom_report = validate_shared_guide(
                    custom_guide, expected_genes=len(records),
                    max_mismatches_per_gene=max_mm,
                    max_seed_mismatches_per_gene=max_seed_mm,
                    min_compatibility=min_compat,
                )
                st.metric("Custom guide validation", custom_report.status)
                st.dataframe(pd.DataFrame([c.__dict__ for c in custom_report.checks]), use_container_width=True, hide_index=True)
                st.markdown("**Closest PAM-compatible target in each gene**")
                st.dataframe(pd.DataFrame(match_rows(custom_guide)), use_container_width=True, hide_index=True)
            except Exception as exc:
                st.error(str(exc))

    st.subheader("Exports")
    export_validation_reports = [
        validate_shared_guide(
            candidate, expected_genes=len(records),
            max_mismatches_per_gene=max_mm,
            max_seed_mismatches_per_gene=max_seed_mm,
            min_compatibility=min_compat,
            panel_hits=st.session_state.get(f"plant_panel_hits::{candidate.spacer}"),
        )
        for candidate in guides
    ]
    export_summary_df = pd.DataFrame([
        {**guide_summary_row(candidate), **validation_summary_row(report)}
        for candidate, report in zip(guides, export_validation_reports)
    ])
    export_summary_df["Input source records"] = "; ".join(f"{g}:{r.source_record_version}" for g, r in records.items())
    export_summary_df["Input assemblies/genomic records"] = "; ".join(f"{g}:{r.assembly}" for g, r in records.items())
    export_summary_df["Input annotation releases"] = "; ".join(f"{g}:{r.annotation_release}" for g, r in records.items())
    export_summary_df["Input sequence SHA-256"] = "; ".join(f"{g}:{r.sequence_sha256}" for g, r in records.items())
    export_summary_df["Search strategy"] = diagnostics.strategy
    export_summary_df["Search upper pair comparisons"] = diagnostics.upper_pair_comparisons
    export_csv = export_summary_df.to_csv(index=False).encode("utf-8")
    export_fasta = "\n".join(
        f">shared_guide_{i+1}|{g.design_type.replace(' ', '_')}|genes={g.genes_covered}|min_compat={g.minimum_compatibility}\n{g.spacer}"
        for i, g in enumerate(guides)
    ).encode("utf-8")
    export_json = json.dumps({
        "app": "Plant MultiGene gRNA Designer",
        "version": APP_VERSION,
        "search": diagnostics.__dict__,
        "genes": {gene: rec.provenance_dict() for gene, rec in records.items()},
        "guides": [
            {**guide_summary_row(g), **validation_summary_row(r), "validation_checks": [c.__dict__ for c in r.checks], "matches": match_rows(g), "notes": g.notes}
            for g, r in zip(guides, export_validation_reports)
        ],
    }, indent=2).encode("utf-8")
    b1, b2, b3 = st.columns(3)
    b1.download_button("Download CSV", export_csv, "plant_shared_guides.csv", "text/csv", use_container_width=True)
    b2.download_button("Download FASTA", export_fasta, "plant_shared_guides.fasta", "text/plain", use_container_width=True)
    b3.download_button("Download JSON", export_json, "plant_shared_guides.json", "application/json", use_container_width=True)

    st.markdown("---")
    st.markdown("**Scientific boundary:** this software prioritizes candidates. Mismatch-aware ranking is best among the generated seed-derived proposals, not a proof of a global optimum over all possible 20-mers. Confirm the exact assembly/genotype, intended coding region, allele/cultivar variation, and genome-wide off-targets before experimental use.")