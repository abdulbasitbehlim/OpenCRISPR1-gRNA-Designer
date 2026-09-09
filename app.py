#!/usr/bin/env python3
"""Professional Streamlit dashboard for OpenCRISPR-1 guide design."""
from __future__ import annotations

import json
import pandas as pd
import plotly.express as px
import streamlit as st

from opencrispr_designer import (
    GuideRNA, design_from_segments, gc_percent, guide_row, guide_format_variants,
    screen_local_reference, sequence_quality_score,
)
from sequence_sources import fetch_gene, manual_record, parse_multifasta
from validation import validate_opencrispr_guide, validation_summary_row

APP_VERSION = "1.3.1"
st.set_page_config(page_title="OpenCRISPR-1 gRNA Designer", page_icon="🧬", layout="wide")

dark = st.sidebar.toggle("Dark mode", value=True)
P = (dict(app="#07111f", side="#0a1628", panel="#0f1d2e", panel2="#13243a", text="#f3f7fb", muted="#9fb0c3", border="#29445f", accent="#60a5fa", accent2="#a78bfa", plot="plotly_dark") if dark else dict(app="#f6f8fc", side="#eef3f9", panel="#ffffff", panel2="#f7f9fc", text="#172033", muted="#64748b", border="#d6deea", accent="#2563eb", accent2="#7c3aed", plot="plotly_white"))
st.markdown(f"""
<style>
:root{{--app:{P['app']};--side:{P['side']};--panel:{P['panel']};--panel2:{P['panel2']};--text:{P['text']};--muted:{P['muted']};--border:{P['border']};--accent:{P['accent']};--accent2:{P['accent2']};}}
html,body,[data-testid="stAppViewContainer"],.stApp{{background:var(--app)!important;color:var(--text)!important}}
[data-testid="stSidebar"]{{background:var(--side)!important;border-right:1px solid var(--border)!important}}
[data-testid="stHeader"]{{background:transparent!important}} .block-container{{max-width:1480px;padding-top:1.25rem;padding-bottom:4rem}}
.stApp h1,.stApp h2,.stApp h3,.stApp h4,.stApp p,.stApp li,.stApp label{{color:var(--text)!important}}
[data-testid="stCaptionContainer"] p{{color:var(--muted)!important}}
.hero{{padding:2.1rem 2.3rem;border:1px solid var(--border);border-radius:22px;background:linear-gradient(135deg,var(--panel),var(--panel2));margin-bottom:1.2rem}}
.hero .eyebrow{{font-size:.75rem;letter-spacing:.15em;text-transform:uppercase;color:var(--accent)!important;font-weight:800}} .hero h1{{font-size:clamp(2rem,4.5vw,4rem);margin:.6rem 0 .8rem;letter-spacing:-.04em;line-height:1}} .hero p{{color:var(--muted)!important;max-width:950px;line-height:1.65}}
.pills{{display:flex;gap:.5rem;flex-wrap:wrap;margin-top:1rem}} .pills span{{border:1px solid var(--border);background:var(--panel2);border-radius:999px;padding:.35rem .65rem;font-size:.78rem}}
[data-testid="stForm"],[data-testid="stExpander"]{{background:var(--panel)!important;border:1px solid var(--border)!important;border-radius:16px!important}}
div[data-testid="stMetric"]{{background:var(--panel);border:1px solid var(--border);padding:.9rem 1rem;border-radius:14px}}
.stButton>button,.stDownloadButton>button{{border-radius:10px;min-height:2.7rem;font-weight:750;border:1px solid var(--border)}} .stButton>button[kind="primary"]{{background:linear-gradient(90deg,var(--accent),var(--accent2));color:white!important;border:0}}
[data-testid="stDataFrame"]{{border:1px solid var(--border);border-radius:12px;overflow:hidden}}
</style>
""", unsafe_allow_html=True)

st.sidebar.markdown(f"### OpenCRISPR-1 Designer\n**NGG guide workbench · v{APP_VERSION}**")
max_guides = st.sidebar.slider("Maximum guides", 5, 100, 30, 5)
min_score = st.sidebar.slider("Minimum sequence-quality score", 0, 95, 35, 5)
with st.sidebar.expander("Validation settings"):
    min_spec_review = st.slider("MIT/Hsu legacy-baseline review threshold", 0, 100, 50, 5, help="Used only after a supplied-FASTA local specificity screen. MIT/Hsu (2013) is a legacy baseline, not a state-of-the-art off-target predictor or OpenCRISPR efficacy cutoff.")
    st.markdown("Hard checks: 20 nt, resolved DNA, NGG PAM, OpenCRISPR compatibility. Review checks: GC 40–60%, sequence-quality threshold, poly-T/homopolymers, and local near matches.")
st.sidebar.info("OpenCRISPR-1 compatibility is kept separate from the heuristic ranking. Current literature is mixed across independent evaluations, so the app does not present OpenCRISPR-1 performance as universally settled or invent an OpenCRISPR-specific efficacy model.")

st.markdown(f"""
<div class="hero"><div class="eyebrow">OpenCRISPR-1</div><h1>Guide RNA <span style="color:{P['accent']}">Designer</span></h1>
<p>Enter a gene and organism or provide reviewed FASTA. The tool retrieves/scans independent sequence segments, finds 20-nt spacers next to an NGG PAM on both strands, ranks sequence quality transparently, and reports 5'-G guide-expression formats that have been experimentally tested with OpenCRISPR-1.</p>
<div class="pills"><span>OpenCRISPR-1</span><span>20 nt + NGG</span><span>Both strands</span><span>GX19 / gX19 / gX20</span><span>Guide validation</span><span>Legacy MIT/Hsu panel screen</span><span>CSV / FASTA / JSON</span></div></div>
""", unsafe_allow_html=True)

st.warning("Evidence status (2026): OpenCRISPR-1 has strong positive founding and follow-up reports, but independent evaluations have reported conflicting on-target/off-target performance. Treat compatibility as a design rule, not a guarantee of generalizable experimental performance; benchmark on your own loci and conditions.")

with st.expander("Quick validation of an existing OpenCRISPR guide"):
    st.caption("Sequence-level validation only: enter a spacer and its PAM. This checks supported OpenCRISPR-1 design rules and expression-format flags, but it cannot confirm genomic locus or whole-genome specificity without target/reference sequence context.")
    qc1, qc2 = st.columns([3, 1])
    with qc1:
        quick_spacer = st.text_input("Existing spacer (20 nt)", key="oc_quick_spacer", max_chars=20).strip().upper()
    with qc2:
        quick_pam = st.text_input("PAM", value="AGG", key="oc_quick_pam", max_chars=3).strip().upper()
    if st.button("Validate existing guide", key="oc_quick_validate"):
        try:
            if len(quick_spacer) != 20:
                raise ValueError("Enter exactly 20 spacer bases.")
            qg = GuideRNA(
                gene="manual_validation", segment_id="sequence_only", spacer=quick_spacer, pam=quick_pam,
                strand="?", start=0, end=0, gc_percent=round(gc_percent(quick_spacer), 1),
                sequence_score=sequence_quality_score(quick_spacer), opencrispr_compatible=True,
            )
            qr = validate_opencrispr_guide(qg, min_sequence_score=min_score, min_specificity_review=min_spec_review)
            st.metric("Validation result", qr.status)
            st.dataframe(pd.DataFrame([c.__dict__ for c in qr.checks]), use_container_width=True, hide_index=True)
        except Exception as exc:
            st.error(str(exc))

mode = st.radio("Input mode", ["Gene lookup", "Manual sequence / FASTA"], horizontal=True)
submitted = False
record = None

if mode == "Gene lookup":
    with st.form("gene_form"):
        c1, c2, c3 = st.columns([1, 1.4, 1])
        with c1: gene = st.text_input("Gene symbol / ID", value="TP53")
        with c2: organism = st.text_input("Organism", value="Homo sapiens")
        with c3: source = st.selectbox("Sequence source", ["NCBI RefSeq", "Ensembl REST"])
        submitted = st.form_submit_button("Design OpenCRISPR guides", type="primary", use_container_width=True)
    if submitted:
        try:
            record = fetch_gene(gene.strip(), organism.strip(), source)
        except Exception as exc:
            st.error(f"Gene retrieval failed: {exc}")
            st.stop()
else:
    with st.form("seq_form"):
        c1, c2 = st.columns(2)
        with c1: gene = st.text_input("Target label", value="target_gene")
        with c2: organism = st.text_input("Organism / sample label", value="manual")
        raw = st.text_area("DNA / RNA / FASTA", height=240, placeholder=">exon1\nACGT...\n>exon2\nACGT...")
        submitted = st.form_submit_button("Design OpenCRISPR guides", type="primary", use_container_width=True)
    if submitted:
        try:
            record = manual_record(raw, gene=gene.strip() or "target_gene", organism=organism)
        except Exception as exc:
            st.error(str(exc)); st.stop()

if submitted and record:
    guides = design_from_segments(record.gene, record.segments, max_guides=max_guides, min_score=min_score)
    st.session_state["oc_record"] = record
    st.session_state["oc_guides"] = guides

if "oc_guides" in st.session_state:
    record = st.session_state["oc_record"]
    guides = st.session_state["oc_guides"]
    st.subheader("Design summary")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Gene", record.gene)
    m2.metric("Scanned segments", len(record.segments))
    m3.metric("Scanned bp", f"{record.total_bp:,}")
    m4.metric("Returned guides", len(guides))
    st.caption(f"Source: {record.source} · {record.accession} · {record.description}")
    with st.expander("Sequence provenance and reproducibility"):
        st.code(
            f"record/version: {record.source_record_version}\n"
            f"assembly/genomic record: {record.assembly}\n"
            f"annotation release: {record.annotation_release}\n"
            f"retrieved UTC: {record.retrieved_at_utc}\n"
            f"sequence SHA-256: {record.sequence_sha256}\n"
            f"ambiguous bases: {record.ambiguity_count} ({', '.join(record.ambiguity_codes) or 'none'})",
            language=None,
        )
    for w in record.warnings: st.warning(w)

    if not guides:
        st.warning("No NGG-compatible 20-nt guide met the selected quality threshold.")
        st.stop()

    validation_reports = []
    for candidate in guides:
        local_state = st.session_state.get(f"oc_local::{candidate.spacer}")
        if local_state is None:
            report = validate_opencrispr_guide(candidate, min_sequence_score=min_score, min_specificity_review=min_spec_review)
        else:
            spec, hits = local_state
            report = validate_opencrispr_guide(
                candidate, min_sequence_score=min_score, local_hits=hits,
                specificity_score=spec, min_specificity_review=min_spec_review,
            )
        validation_reports.append(report)
    df = pd.DataFrame([
        {**guide_row(candidate), **validation_summary_row(report)}
        for candidate, report in zip(guides, validation_reports)
    ])
    st.markdown("#### Ranked guides with validation")
    st.dataframe(df, use_container_width=True, hide_index=True)
    c1, c2 = st.columns(2)
    with c1:
        fig = px.scatter(df, x="GC%", y="Sequence quality", hover_name="Spacer (20 nt)", symbol="Strand", title="Candidate guide landscape")
        fig.update_layout(template=P["plot"], height=390, margin=dict(l=20,r=20,t=50,b=20)); st.plotly_chart(fig, use_container_width=True)
    with c2:
        top = df.head(20)
        fig2 = px.bar(top, x="Spacer (20 nt)", y="Sequence quality", pattern_shape="Strand", title="Top OpenCRISPR-compatible spacers")
        fig2.update_layout(template=P["plot"], height=390, margin=dict(l=20,r=20,t=50,b=20), xaxis_tickangle=-45); st.plotly_chart(fig2, use_container_width=True)

    labels = [f"{i+1}. {g.spacer} · {g.pam} · {g.strand} · score {g.sequence_score:.0f}" for i,g in enumerate(guides)]
    selected_label = st.selectbox("Inspect one guide", labels)
    g = guides[labels.index(selected_label)]
    st.markdown(f"### `{g.spacer}` + `{g.pam}`")
    d1,d2,d3,d4 = st.columns(4)
    d1.metric("OpenCRISPR-1", "Compatible")
    d2.metric("GC", f"{g.gc_percent:.1f}%")
    d3.metric("Sequence quality", f"{g.sequence_score:.1f}")
    d4.metric("Target strand", g.strand)

    selected_local = st.session_state.get(f"oc_local::{g.spacer}")
    if selected_local is None:
        selected_validation = validate_opencrispr_guide(g, min_sequence_score=min_score, min_specificity_review=min_spec_review)
    else:
        selected_spec, selected_hits = selected_local
        selected_validation = validate_opencrispr_guide(
            g, min_sequence_score=min_score, local_hits=selected_hits, specificity_score=selected_spec,
            min_specificity_review=min_spec_review,
        )
    st.markdown("#### Validation")
    v1,v2,v3,v4 = st.columns(4)
    v1.metric("Overall", selected_validation.status)
    v2.metric("PASS checks", selected_validation.pass_count)
    v3.metric("REVIEW checks", selected_validation.review_count)
    v4.metric("FAIL checks", selected_validation.fail_count)
    st.dataframe(pd.DataFrame([c.__dict__ for c in selected_validation.checks]), use_container_width=True, hide_index=True)
    if selected_validation.status == "PASS":
        st.success("Core OpenCRISPR-1 validation passed. Specificity status: " + selected_validation.specificity_status)
    elif selected_validation.status == "REVIEW":
        st.warning("No hard compatibility failure, but one or more checks need review. Specificity status: " + selected_validation.specificity_status)
    else:
        st.error("One or more hard validation checks failed. Resolve them before prioritizing this guide.")

    st.markdown("#### Guide-expression formats reported with OpenCRISPR-1")
    fmt_df = pd.DataFrame([x.__dict__ for x in guide_format_variants(g.spacer)])
    st.dataframe(fmt_df, use_container_width=True, hide_index=True)
    st.caption("These labels describe 5'-G guide formats evaluated in the 2026 OpenCRISPR study; they are not a claim that one format is universally optimal.")

    with st.expander("Optional local-reference specificity screen"):
        st.caption("Upload/paste an intended locus plus nearby or suspected off-target sequences. MIT/Hsu (2013) is retained only as a transparent legacy baseline; later empirical scores such as CFD and newer ML models can outperform it. This is not a whole-genome specificity result.")
        ref = st.text_area("Reference FASTA", height=170, key="oc_ref")
        radius = st.slider("Maximum mismatches", 0, 4, 3, key="oc_mm")
        if st.button("Run local specificity", key="oc_screen"):
            try:
                panel = parse_multifasta(ref)
                spec, hits = screen_local_reference(g.spacer, panel, radius)
                g.specificity_score = spec; g.specificity_method = "MIT/Hsu 2013 legacy baseline (supplied FASTA panel)"; g.off_target_count = len(hits)
                st.session_state[f"oc_local::{g.spacer}"] = (spec, hits)
                st.metric("MIT/Hsu legacy specificity (supplied panel)", spec)
                st.info("Interpret this as a legacy baseline only. For experimental prioritization, use a genome-aware workflow and, where applicable, more modern off-target scoring/validation such as CFD-informed or validated ML/empirical methods.")
                screened_validation = validate_opencrispr_guide(
                    g, min_sequence_score=min_score, local_hits=hits, specificity_score=spec,
                    min_specificity_review=min_spec_review,
                )
                if screened_validation.specificity_status.startswith("FAIL"):
                    st.error("Local validation: " + screened_validation.specificity_status)
                elif screened_validation.specificity_status.startswith("REVIEW"):
                    st.warning("Local validation: " + screened_validation.specificity_status)
                else:
                    st.success("Local validation: " + screened_validation.specificity_status)
                hit_df = pd.DataFrame([{**h.__dict__, "mismatch_positions": ",".join(map(str,h.mismatch_positions)) or "Exact", "mit_pair_risk": round(h.mit_pair_risk*100,2)} for h in hits])
                if hit_df.empty: st.success("No retained near matches after excluding one exact intended target.")
                else: st.dataframe(hit_df, use_container_width=True, hide_index=True)
            except Exception as exc:
                st.error(str(exc))

    st.subheader("Exports")
    export_reports = []
    for candidate in guides:
        export_state = st.session_state.get(f"oc_local::{candidate.spacer}")
        if export_state is None:
            export_report = validate_opencrispr_guide(candidate, min_sequence_score=min_score, min_specificity_review=min_spec_review)
        else:
            export_spec, export_hits = export_state
            export_report = validate_opencrispr_guide(
                candidate, min_sequence_score=min_score, local_hits=export_hits,
                specificity_score=export_spec, min_specificity_review=min_spec_review,
            )
        export_reports.append(export_report)
    out_df = pd.DataFrame([
        {**guide_row(candidate), **validation_summary_row(report)}
        for candidate, report in zip(guides, export_reports)
    ])
    out_df["Input source record/version"] = record.source_record_version
    out_df["Input assembly/genomic record"] = record.assembly
    out_df["Input annotation release"] = record.annotation_release
    out_df["Input retrieved UTC"] = record.retrieved_at_utc
    out_df["Input sequence SHA-256"] = record.sequence_sha256
    out_df["Evidence status"] = "2025–2026 independent evaluations are mixed; compatibility is not an efficacy guarantee."
    out_df["Specificity method note"] = "MIT/Hsu 2013 is a legacy supplied-panel baseline, not state-of-the-art or genome-wide specificity."
    csv = out_df.to_csv(index=False).encode()
    fasta = "\n".join(f">opencrispr_guide_{i+1}|pam={x.pam}|strand={x.strand}|score={x.sequence_score}\n{x.spacer}" for i,x in enumerate(guides)).encode()
    js = json.dumps({
        "app":"OpenCRISPR-1 gRNA Designer",
        "version":APP_VERSION,
        "target": record.provenance_dict(),
        "evidence_status": "OpenCRISPR-1 performance remains context-dependent; 2026 literature includes conflicting independent evaluations. Compatibility is not an efficacy guarantee.",
        "specificity_method_note": "MIT/Hsu 2013 is a legacy supplied-panel baseline, not state-of-the-art or genome-wide specificity.",
        "guides":[
            {**guide_row(candidate), **validation_summary_row(report), "validation_checks":[c.__dict__ for c in report.checks]}
            for candidate, report in zip(guides, export_reports)
        ],
    }, indent=2).encode()
    b1,b2,b3 = st.columns(3)
    b1.download_button("Download CSV", csv, "opencrispr_guides.csv", "text/csv", use_container_width=True)
    b2.download_button("Download FASTA", fasta, "opencrispr_guides.fasta", "text/plain", use_container_width=True)
    b3.download_button("Download JSON", js, "opencrispr_guides.json", "application/json", use_container_width=True)

    st.markdown("---")
    st.markdown("**Scientific boundary:** candidate ranking does not replace genomic/exon validation, variant review, or genome-wide off-target analysis. OpenCRISPR-1 guide compatibility is based on published NGG/guide-architecture evidence, but 2026 studies report conflicting performance across experimental settings. The local MIT/Hsu score is a legacy supplied-panel baseline, not state-of-the-art specificity prediction.")
