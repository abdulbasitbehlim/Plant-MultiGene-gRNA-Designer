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
.hero{{padding:26px 28px;border:1px solid var(--border);border-radius:24px;background:linear-gradient(135deg,var(--panel),var(--panel2));margin-bottom:18px}}
.hero h1{{font-size:2.25rem;margin:0 0 .35rem 0}}
.hero p{{font-size:1.03rem;margin:.2rem 0;color:var(--muted)!important;max-width:1000px}}
.pill{{display:inline-block;border:1px solid var(--border);border-radius:999px;padding:5px 10px;margin:4px 6px 0 0;color:var(--text);font-size:.82rem}}
[data-testid="stMetric"]{{background:var(--panel);border:1px solid var(--border);border-radius:16px;padding:14px}}
div[data-testid="stDataFrame"]{{border:1px solid var(--border);border-radius:14px;overflow:hidden}}
.stButton>button,.stDownloadButton>button{{border-radius:12px!important;font-weight:700!important}}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero">
<h1>Plant MultiGene gRNA Designer</h1>
<p>Find one SpCas9 spacer that can intentionally target every gene in a user-defined plant set. Exact shared guides and mismatch-aware consensus candidates are kept separate, every candidate is validated, and sequence provenance is exported.</p>
<span class="pill">20 nt + NGG</span><span class="pill">Both strands</span><span class="pill">Exact + mismatch-aware</span><span class="pill">PASS / REVIEW / FAIL</span>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.markdown("### Design settings")
    input_mode = st.radio("Input", ["Gene lookup", "Manual multi-FASTA"], horizontal=True)
    source = st.selectbox("Sequence source", ["NCBI RefSeq", "Ensembl REST"], disabled=input_mode != "Gene lookup")
    organism = st.text_input("Plant organism", value="Arabidopsis thaliana")
    allow_mismatch = st.toggle("Enable mismatch-aware consensus", value=True)
    max_mm = st.slider("Maximum spacer mismatches per gene", 0, 3, 1, disabled=not allow_mismatch)
    max_seed_mm = st.slider("Maximum PAM-proximal seed mismatches per gene", 0, 2, 1, disabled=not allow_mismatch)
    min_compat = st.slider("Minimum compatibility proxy", 0, 100, 70, 1, disabled=not allow_mismatch)
    max_results = st.slider("Maximum guides", 1, 50, 20)
    st.caption(f"Recommended: 2–{RECOMMENDED_MAX_GENES} genes · hard cap: {MAX_GENES} · mismatch guard: {MAX_PAIR_COMPARISONS:,} estimated pair comparisons")

if input_mode == "Gene lookup":
    gene_text = st.text_area("Gene symbols / IDs (one per line)", value="AT4G18197\nAT4G18195\nAT4G18205", height=120)
    submit = st.button("Fetch genes and design guides", type="primary", use_container_width=True)
else:
    fasta_text = st.text_area("Multi-FASTA (one record per gene)", height=220, placeholder=">geneA\nACGT...\n>geneB\nACGT...")
    submit = st.button("Design from FASTA", type="primary", use_container_width=True)

if submit:
    st.session_state.pop("plant_records", None)
    st.session_state.pop("plant_guides", None)
    try:
        if input_mode == "Gene lookup":
            genes = [x.strip() for x in gene_text.splitlines() if x.strip()]
            if len(genes) < 2:
                raise ValueError("Enter at least two genes for a shared-guide design.")
            if len(genes) > MAX_GENES:
                raise ValueError(f"This interactive tool supports at most {MAX_GENES} genes per run.")
            if len(genes) > RECOMMENDED_MAX_GENES:
                st.warning(f"{len(genes)} genes exceeds the recommended interactive range. The hard workload guard will stop very large mismatch searches.")
            records = {}
            progress = st.progress(0, text="Retrieving gene sequences…")
            for i, gene in enumerate(genes):
                records[gene] = fetch_gene(gene, organism, source=source)
                progress.progress((i + 1) / len(genes), text=f"Retrieved {gene}")
            progress.empty()
        else:
            parsed = parse_multifasta(fasta_text)
            if len(parsed) < 2:
                raise ValueError("Provide at least two FASTA records.")
            if len(parsed) > MAX_GENES:
                raise ValueError(f"This interactive tool supports at most {MAX_GENES} genes per run.")
            records = manual_records(parsed)

        guides = design_shared_guides(
            records,
            allow_mismatch=allow_mismatch,
            max_mismatches_per_gene=max_mm,
            max_seed_mismatches_per_gene=max_seed_mm,
            min_compatibility=min_compat,
            max_results=max_results,
        )
        st.session_state["plant_records"] = records
        st.session_state["plant_guides"] = guides
        st.session_state["plant_settings"] = dict(max_mm=max_mm, max_seed_mm=max_seed_mm, min_compat=min_compat)
    except Exception as exc:
        st.error(str(exc))

records = st.session_state.get("plant_records")
guides = st.session_state.get("plant_guides")
settings = st.session_state.get("plant_settings", {})

if records:
    with st.expander("Sequence provenance, reproducibility and warnings", expanded=False):
        rows = []
        for gene, rec in records.items():
            rows.append({
                "Gene": gene,
                "Source": rec.source,
                "Source record/version": rec.source_record_version,
                "Assembly / genomic record": rec.assembly,
                "Annotation release/date": rec.annotation_release,
                "Retrieved UTC": rec.retrieved_utc,
                "SHA-256": rec.sequence_sha256,
                "Ambiguity count": rec.ambiguity_count,
                "Ambiguity codes": ",".join(rec.ambiguity_codes),
                "Warnings": "; ".join(rec.warnings),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

if records is not None and guides is not None:
    diagnostics = estimate_search_diagnostics(records)
    with st.expander("Search strategy and workload", expanded=False):
        st.write(diagnostics.strategy)
        st.write(f"PAM-compatible sites: **{diagnostics.total_sites:,}** · unique observed seeds: **{diagnostics.unique_seed_spacers:,}** · estimated upper pair comparisons: **{diagnostics.upper_pair_comparisons:,}**")

    if not guides:
        st.warning("No valid shared guide satisfied the current settings. The tool does not force a result. Consider reviewing the input sequences or carefully relaxing mismatch settings.")
    else:
        max_mm = settings.get("max_mm", max_mm)
        max_seed_mm = settings.get("max_seed_mm", max_seed_mm)
        min_compat = settings.get("min_compat", min_compat)
        reports = [validate_shared_guide(g, max_mm, max_seed_mm, min_compat) for g in guides]
        summary_df = pd.DataFrame([
            {**guide_summary_row(g), **validation_summary_row(r)} for g, r in zip(guides, reports)
        ])

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Candidates", len(guides))
        c2.metric("Exact shared", sum(g.design_type == "Exact shared" for g in guides))
        c3.metric("PASS", sum(r.status == "PASS" for r in reports))
        c4.metric("Genes", len(records))

        st.subheader("Ranked guides")
        st.dataframe(summary_df, use_container_width=True, hide_index=True)

        selected_idx = st.selectbox("Inspect guide", range(len(guides)), format_func=lambda i: f"#{i+1} · {guides[i].spacer} · {guides[i].design_type}")
        selected = guides[selected_idx]
        report = reports[selected_idx]

        left, right = st.columns([1.25, 1])
        with left:
            st.markdown("#### Per-gene target evidence")
            st.dataframe(pd.DataFrame(match_rows(selected)), use_container_width=True, hide_index=True)
        with right:
            st.markdown("#### Validation")
            badge = {"PASS":"✅", "REVIEW":"⚠️", "FAIL":"❌"}[report.status]
            st.markdown(f"### {badge} {report.status}")
            st.dataframe(pd.DataFrame([c.__dict__ for c in report.checks]), use_container_width=True, hide_index=True)

        st.markdown("#### Validate a custom 20-nt guide against these genes")
        custom = st.text_input("Custom spacer", max_chars=40, placeholder="20 nt DNA spacer")
        if st.button("Validate custom guide"):
            candidate = re.sub(r"\s+", "", custom.upper().replace("U", "T"))
            try:
                custom_guide = evaluate_candidate_guide(candidate, records)
                custom_report = validate_shared_guide(custom_guide, max_mm, max_seed_mm, min_compat)
                st.write(f"**{custom_report.status}** · {custom_guide.genes_covered}/{len(records)} genes covered · exact genes {custom_guide.exact_gene_count}")
                st.dataframe(pd.DataFrame(match_rows(custom_guide)), use_container_width=True, hide_index=True)
            except Exception as exc:
                st.error(str(exc))

        st.markdown("#### Optional supplied-FASTA panel screen")
        st.caption("This is a local panel review only, not whole-genome specificity. Intended paralog targets may also appear as hits and must be interpreted by the user.")
        panel = st.text_area("Reference/paralog FASTA", height=130, key="panel_fasta")
        if st.button("Screen selected guide against panel") and panel.strip():
            hits = screen_reference_panel(selected.spacer, panel)
            st.session_state[f"plant_panel_hits::{selected.spacer}"] = hits
            if hits:
                st.dataframe(pd.DataFrame([h.__dict__ for h in hits]), use_container_width=True, hide_index=True)
            else:
                st.success("No PAM-compatible ≤3-mismatch panel hits were found.")

        export_validation_reports = [
            validate_shared_guide(
                candidate,
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
