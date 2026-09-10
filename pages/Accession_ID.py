#!/usr/bin/env python3
"""Accession-ID input workflow for Plant MultiGene gRNA Designer v1.3.1."""
from __future__ import annotations

import re
import json
import pandas as pd
import streamlit as st

from accession_sources import fetch_accession
from plant_multiguide import (
    MAX_GENES,
    RECOMMENDED_MAX_GENES,
    design_shared_guides,
    estimate_search_diagnostics,
    guide_summary_row,
    match_rows,
)
from validation import validate_shared_guide, validation_summary_row

APP_VERSION = "1.3.1"

st.set_page_config(page_title="Accession ID · Plant MultiGene", page_icon="🧬", layout="wide")
st.title("Accession ID input")
st.caption(f"Plant MultiGene gRNA Designer · v{APP_VERSION}")
st.info(
    "Use this third input method when you already know the database accession/stable ID. "
    "For NCBI, enter nucleotide/RefSeq accessions such as NM_... or XM_.... "
    "For Ensembl, enter stable gene or transcript IDs."
)

with st.form("plant_accession_form"):
    source = st.selectbox("Accession source", ["NCBI RefSeq / Nucleotide", "Ensembl REST"])
    raw_ids = st.text_area(
        "Gene accession IDs",
        height=150,
        placeholder="Example:\nNM_...\nXM_...\n\nor Ensembl stable IDs, one per line",
        help="Enter at least two accessions because this tool searches for one guide shared across multiple genes. Separate IDs with commas, spaces, semicolons, or new lines.",
    )
    include_mismatch = st.toggle("Allow mismatch-aware consensus", value=True)
    c1, c2, c3 = st.columns(3)
    with c1:
        max_mm = st.slider("Max mismatches per gene", 0, 3, 2, disabled=not include_mismatch)
    with c2:
        max_seed_mm = st.slider("Max seed mismatches", 0, 2, 1, disabled=not include_mismatch)
    with c3:
        min_compat = st.slider("Minimum compatibility proxy", 0, 100, 55, 5, disabled=not include_mismatch)
    max_results = st.slider("Maximum returned guides", 10, 200, 50, 10)
    submitted = st.form_submit_button("Fetch accessions and design shared guides", type="primary", use_container_width=True)

if submitted:
    accession_ids = [x for x in re.split(r"[\s,;]+", raw_ids.strip()) if x]
    accession_ids = list(dict.fromkeys(accession_ids))
    if len(accession_ids) < 2:
        st.error("Enter at least two distinct accession IDs.")
        st.stop()
    if len(accession_ids) > MAX_GENES:
        st.error(f"At most {MAX_GENES} accessions are supported per interactive run.")
        st.stop()
    if len(accession_ids) > RECOMMENDED_MAX_GENES:
        st.warning(f"The recommended interactive range is 2–{RECOMMENDED_MAX_GENES} genes.")

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
        progress.progress(i / len(accession_ids), text=f"Retrieved {i}/{len(accession_ids)} accessions")
    progress.empty()

    if errors:
        st.error("One or more accessions could not be resolved:\n\n" + "\n".join(f"- {e}" for e in errors))
        st.stop()
    if len(records) < 2:
        st.error("Fewer than two usable gene records were retrieved.")
        st.stop()

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
    st.subheader("Accession retrieval summary")
    for label, rec in records.items():
        st.markdown(f"**{label}** — {rec.source} · requested/resolved accession `{rec.accession}`")
        st.caption(f"{rec.description} · {rec.total_bp:,} bp across {len(rec.segments)} scanned segment(s)")
        with st.expander(f"Provenance · {label}"):
            st.json(rec.provenance_dict())
            for warning in rec.warnings:
                st.warning(warning)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Genes", len(records))
    m2.metric("PAM-compatible sites", sum(len(v) for v in sites_by_gene.values()))
    m3.metric("Shared guides", len(guides))
    m4.metric("Exact shared", sum(g.design_type == "Exact shared" for g in guides))

    if not guides:
        st.warning("No single guide satisfied the selected all-gene constraints. This is a valid result.")
        st.stop()

    reports = [
        validate_shared_guide(
            g,
            expected_genes=len(records),
            max_mismatches_per_gene=max_mm,
            max_seed_mismatches_per_gene=max_seed_mm,
            min_compatibility=min_compat,
        )
        for g in guides
    ]
    df = pd.DataFrame([
        {**guide_summary_row(g), **validation_summary_row(r)}
        for g, r in zip(guides, reports)
    ])
    st.subheader("Ranked shared guides")
    st.dataframe(df, use_container_width=True, hide_index=True)

    labels = [f"{i+1}. {g.spacer} · {g.design_type}" for i, g in enumerate(guides)]
    selected = guides[labels.index(st.selectbox("Inspect one guide", labels))]
    st.markdown("#### Per-gene target evidence")
    st.dataframe(pd.DataFrame(match_rows(selected)), use_container_width=True, hide_index=True)

    df["Input accessions"] = "; ".join(accession_ids)
    df["Input source records"] = "; ".join(f"{g}:{r.source_record_version}" for g, r in records.items())
    df["Input sequence SHA-256"] = "; ".join(f"{g}:{r.sequence_sha256}" for g, r in records.items())
    st.download_button(
        "Download CSV",
        df.to_csv(index=False).encode("utf-8"),
        "plant_shared_guides_from_accessions.csv",
        "text/csv",
        use_container_width=True,
    )

    export = {
        "app": "Plant MultiGene gRNA Designer",
        "version": APP_VERSION,
        "input_mode": "Accession ID",
        "source": source,
        "requested_accessions": accession_ids,
        "search": diagnostics.__dict__,
        "records": {k: v.provenance_dict() for k, v in records.items()},
    }
    st.download_button(
        "Download accession provenance JSON",
        json.dumps(export, indent=2).encode("utf-8"),
        "plant_accession_provenance.json",
        "application/json",
        use_container_width=True,
    )

st.markdown("---")
st.caption(
    "Scientific boundary: accession retrieval supplies sequence input; it does not validate the biological identity, cultivar/allele, intended coding region, genome-wide off-target profile, or experimental editing efficacy."
)
