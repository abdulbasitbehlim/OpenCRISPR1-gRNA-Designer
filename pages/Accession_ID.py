#!/usr/bin/env python3
"""Accession-ID input workflow for OpenCRISPR-1 gRNA Designer v1.3.1."""
from __future__ import annotations

import json
import pandas as pd
import streamlit as st

from accession_sources import fetch_accession
from opencrispr_designer import design_from_segments, guide_row
from validation import validate_opencrispr_guide, validation_summary_row

APP_VERSION = "1.3.1"

st.set_page_config(page_title="Accession ID · OpenCRISPR-1", page_icon="🧬", layout="wide")
st.title("Accession ID input")
st.caption(f"OpenCRISPR-1 gRNA Designer · v{APP_VERSION}")
st.info(
    "Use this third input method when you already know the database accession/stable ID. "
    "For NCBI, enter a nucleotide/RefSeq accession such as NM_... or XM_.... "
    "For Ensembl, enter a stable gene or transcript ID such as ENSG... or ENST...."
)

with st.form("opencrispr_accession_form"):
    source = st.selectbox("Accession source", ["NCBI RefSeq / Nucleotide", "Ensembl REST"])
    accession = st.text_input("Gene accession ID", placeholder="Example: NM_000546.6 or ENST00000269305")
    c1, c2 = st.columns(2)
    with c1:
        max_guides = st.slider("Maximum guides", 5, 100, 30, 5)
    with c2:
        min_score = st.slider("Minimum sequence-quality score", 0, 95, 35, 5)
    submitted = st.form_submit_button("Fetch accession and design guides", type="primary", use_container_width=True)

if submitted:
    try:
        record = fetch_accession(accession, source)
    except Exception as exc:
        st.error(f"Accession retrieval failed: {exc}")
        st.stop()

    guides = design_from_segments(record.gene, record.segments, max_guides=max_guides, min_score=min_score)

    st.subheader("Accession retrieval summary")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Gene / target", record.gene)
    m2.metric("Scanned segments", len(record.segments))
    m3.metric("Scanned bp", f"{record.total_bp:,}")
    m4.metric("Returned guides", len(guides))

    st.caption(f"Source: {record.source} · accession `{record.accession}` · {record.description}")
    with st.expander("Sequence provenance and reproducibility"):
        st.json(record.provenance_dict())
    for warning in record.warnings:
        st.warning(warning)

    if not guides:
        st.warning("No NGG-compatible 20-nt guide met the selected sequence-quality threshold.")
        st.stop()

    reports = [
        validate_opencrispr_guide(g, min_sequence_score=min_score)
        for g in guides
    ]
    df = pd.DataFrame([
        {**guide_row(g), **validation_summary_row(r)}
        for g, r in zip(guides, reports)
    ])

    st.subheader("Ranked OpenCRISPR-1-compatible guides")
    st.dataframe(df, use_container_width=True, hide_index=True)

    df["Input mode"] = "Accession ID"
    df["Requested accession"] = accession.strip()
    df["Resolved source record/version"] = record.source_record_version
    df["Input sequence SHA-256"] = record.sequence_sha256
    st.download_button(
        "Download CSV",
        df.to_csv(index=False).encode("utf-8"),
        "opencrispr_guides_from_accession.csv",
        "text/csv",
        use_container_width=True,
    )

    provenance = {
        "app": "OpenCRISPR-1 gRNA Designer",
        "version": APP_VERSION,
        "input_mode": "Accession ID",
        "source": source,
        "requested_accession": accession.strip(),
        "target": record.provenance_dict(),
    }
    st.download_button(
        "Download accession provenance JSON",
        json.dumps(provenance, indent=2).encode("utf-8"),
        "opencrispr_accession_provenance.json",
        "application/json",
        use_container_width=True,
    )

st.markdown("---")
st.caption(
    "Scientific boundary: accession retrieval supplies the sequence used for guide discovery. It does not guarantee locus identity, biological efficacy, or genome-wide specificity; review genomic context and validate experimentally."
)
