#!/usr/bin/env python3
"""Professional Streamlit dashboard for OpenCRISPR-1 guide design."""
from __future__ import annotations

import json
import pandas as pd
import plotly.express as px
import streamlit as st

from opencrispr_designer import (
    GuideRNA, design_from_segments, gc_percent, guide_row, guide_format_variants,
    sequence_quality_score,
)
from sequence_sources import fetch_gene, manual_record, parse_multifasta
from accession_sources import fetch_accession
from validation import validate_opencrispr_guide, validation_summary_row

from local_screening import TargetLocus, screen_reference
from workflow import APP_VERSION, clear_design_state, guide_key, reset_panel_if_changed, export_bundle

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
[data-baseweb="input"] input,[data-baseweb="textarea"] textarea,.stApp input,.stApp textarea{{color:var(--text)!important;-webkit-text-fill-color:var(--text)!important}}
[data-baseweb="select"] span,[data-baseweb="select"] div,[data-baseweb="select"] svg{{color:var(--text)!important;fill:var(--text)!important}}
.stApp input::placeholder,.stApp textarea::placeholder{{color:var(--muted)!important;opacity:.85!important}}
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
<p>Enter a gene and organism, provide a database accession ID, or paste reviewed FASTA. The tool retrieves/scans independent sequence segments, finds 20-nt spacers next to an NGG PAM on both strands, ranks sequence quality transparently, and reports 5'-G guide-expression formats that have been experimentally tested with OpenCRISPR-1.</p>
<div class="pills"><span>OpenCRISPR-1</span><span>20 nt + NGG</span><span>Both strands</span><span>Accession ID</span><span>GX19 / gX19 / gX20</span><span>Guide validation</span><span>Legacy MIT/Hsu panel screen</span><span>CSV / FASTA / JSON</span></div></div>
""", unsafe_allow_html=True)

st.warning("Evidence status (2026): OpenCRISPR-1 has strong positive founding and follow-up reports, but independent evaluations have reported conflicting on-target/off-target performance. Treat compatibility as a design rule, not a guarantee of generalizable experimental performance; benchmark on your own loci and conditions.")

with st.expander("Quick validation of an existing OpenCRISPR guide"):
    st.caption("Sequence-level validation only: enter a spacer and its PAM. This checks supported OpenCRISPR-1 design rules and expression-format flags, but it cannot confirm genomic locus or whole-genome specificity without target/reference sequence context.")
    qc1, qc2 = st.columns([3, 1])
    with qc1:
        quick_spacer = st.text_input("Existing spacer (20 nt)", key="oc_quick_spacer", max_chars=20).strip().upper()
    with qc2:
        quick_pam = st.text_input("PAM", value="", key="oc_quick_pam", max_chars=3).strip().upper()
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

def load_offline_example():
    clear_design_state(st.session_state)
    st.session_state["oc_input_mode"] = "Manual sequence / FASTA"
    st.session_state["oc_target_raw"] = ">demo_forward\nACGTACGTACGTACGTACGTAGG\n>demo_reverse\nCCATCAGTCAGTCAGTCAGTCAG"

st.button("Load offline example", on_click=load_offline_example)
mode = st.radio("Input mode", ["Gene lookup", "Accession ID", "Manual sequence / FASTA"], horizontal=True, key="oc_input_mode")
submitted = False
record = None

if mode == "Gene lookup":
    with st.form("gene_form"):
        c1, c2, c3 = st.columns([1, 1.4, 1])
        with c1: gene = st.text_input("Gene symbol / ID", value="TP53")
        with c2: organism = st.text_input("Organism", value="Homo sapiens")
        with c3: source = st.selectbox("Sequence source", ["Ensembl REST", "NCBI RefSeq"])
        submitted = st.form_submit_button("Design OpenCRISPR guides", type="primary", use_container_width=True)
    if submitted:
        clear_design_state(st.session_state)
        try:
            record = fetch_gene(gene.strip(), organism.strip(), source)
        except Exception as exc:
            st.error(f"Gene retrieval failed: {exc}")
            st.stop()
elif mode == "Accession ID":
    with st.form("accession_form"):
        c1, c2 = st.columns([2, 1])
        with c1:
            accession = st.text_input(
                "Gene accession ID",
                placeholder="Example: ENST00000269305 (Ensembl source)",
                help="Choose the matching source. NCBI RNA records require exon annotations; transcripts without boundaries are rejected. Ensembl transcript IDs provide independent exons.",
            )
        with c2:
            source = st.selectbox("Accession source", ["NCBI RefSeq / Nucleotide", "Ensembl REST"])
        submitted = st.form_submit_button("Fetch accession and design OpenCRISPR guides", type="primary", use_container_width=True)
    if submitted:
        clear_design_state(st.session_state)
        try:
            record = fetch_accession(accession.strip(), source)
        except Exception as exc:
            st.error(f"Accession retrieval failed: {exc}")
            st.stop()
else:
    with st.form("seq_form"):
        c1, c2 = st.columns(2)
        with c1: gene = st.text_input("Target label", value="target_gene")
        with c2: organism = st.text_input("Organism / sample label", value="manual")
        raw = st.text_area("DNA / RNA / FASTA", key="oc_target_raw", height=240, placeholder=">exon1\nACGT...\n>exon2\nACGT...")
        submitted = st.form_submit_button("Design OpenCRISPR guides", type="primary", use_container_width=True)
    if submitted:
        clear_design_state(st.session_state)
        try:
            record = manual_record(raw, gene=gene.strip() or "target_gene", organism=organism)
        except Exception as exc:
            st.error(str(exc)); st.stop()

if submitted and record:
    try:
        guides = design_from_segments(record.gene, record.segments, max_guides=max_guides, min_score=min_score)
    except ValueError as exc:
        st.error(str(exc))
        st.stop()
    st.session_state["oc_record"] = record
    st.session_state["oc_guides"] = guides
    st.session_state["oc_settings"] = {"max_guides": max_guides, "min_score": min_score, "min_spec_review": min_spec_review}
    st.session_state["oc_screens"] = {}

if "oc_guides" in st.session_state:
    record = st.session_state["oc_record"]
    guides = st.session_state["oc_guides"]
    settings = st.session_state["oc_settings"]
    if settings != {"max_guides": max_guides, "min_score": min_score, "min_spec_review": min_spec_review}:
        st.warning("Design settings changed. Submit the input again to update the results and exports.")
        st.stop()
    st.subheader("Design summary")
    cols = st.columns(4)
    for col, label, value in zip(cols, ["Gene", "Scanned segments", "Scanned bp", "Returned guides"],
                                  [record.gene, len(record.segments), f"{record.total_bp:,}", len(guides)]):
        col.metric(label, value)
    st.caption(f"Source: {record.source} · {record.accession} · {record.description}")
    st.caption("Coordinates are 1-based and inclusive within each supplied segment. They are not chromosome coordinates.")
    with st.expander("Sequence provenance and reproducibility"):
        st.caption(f"sequence SHA-256: {record.sequence_sha256}")
        st.dataframe(
            key_value_table(record.provenance_dict(), "Provenance field", "Value"),
            hide_index=True,
            use_container_width=True,
        )
    for warning in record.warnings:
        st.warning(warning)
    if not guides:
        st.warning("No NGG-compatible 20-nt guide met the selected quality threshold.")
        st.stop()

    selected = st.selectbox("Inspect one guide", range(len(guides)),
                            format_func=lambda i: f"{i+1}. {guides[i].spacer} · {guides[i].segment_id}:{guides[i].start} ({guides[i].strand})")
    g = guides[selected]
    with st.expander("Optional local-reference specificity screen"):
        st.caption("Supply up to 2 million bases. The screen searches the unmodified 20-base targeting spacer on both strands at NGG sites, allowing 0–4 substitutions. It does not screen a complete sgRNA or model the effects of an added or substituted 5′ G.")
        ref = st.text_area("Reference FASTA", height=170, key="oc_ref")
        radius = st.slider("Maximum mismatches", 0, 4, 3, key="oc_mm")
        display_limit = st.number_input("Maximum displayed hits", min_value=1, max_value=5000, value=250)
        exclude = st.checkbox("Exclude a verified intended locus", value=False,
                              help="Only this exact contig, start and strand will be excluded, after the spacer and PAM are checked.")
        intended = None
        if exclude:
            c1, c2, c3 = st.columns(3)
            contig = c1.text_input("Intended reference FASTA ID")
            start = c2.number_input("Intended spacer start", min_value=1, value=1)
            strand = c3.selectbox("Intended spacer strand", ["+", "-"])
            intended = TargetLocus(contig.strip(), int(start), strand)
        reset_panel_if_changed(st.session_state, ref, radius, intended, int(display_limit))
        if st.button("Run local specificity", key="oc_screen"):
            key = guide_key(record, g)
            st.session_state["oc_screens"].pop(key, None)
            try:
                result = screen_reference(g.spacer, parse_multifasta(ref), radius, intended, int(display_limit))
                st.session_state["oc_screens"][key] = result
            except ValueError as exc:
                st.error(str(exc))
        result = st.session_state["oc_screens"].get(guide_key(record, g))
        if result:
            st.metric("MIT/Hsu legacy specificity (supplied panel)", result.specificity_score if result.specificity_score is not None else "Not available")
            st.caption(f"{result.total_hits} total hits · {len(result.hits)} displayed · {result.scanned_sites} searchable NGG sites · reference SHA-256 {result.panel_sha256}")
            for warning in result.warnings:
                st.warning(warning)
            if result.hits:
                st.dataframe(pd.DataFrame([h.__dict__ for h in result.hits]), use_container_width=True, hide_index=True)
            else:
                st.info("No matching hits were retained within this panel and mismatch radius. Review the search scope and intended-locus status.")

    screens = st.session_state["oc_screens"]
    rows, bundle = export_bundle(record, guides, screens, settings)
    df = pd.DataFrame(rows)
    st.markdown("#### Ranked guides with validation")
    st.dataframe(df, use_container_width=True, hide_index=True)
    fig = px.scatter(df, x="GC%", y="Sequence quality", hover_name="Spacer (20 nt)", symbol="Strand")
    fig.update_layout(template=P["plot"], height=340, margin=dict(l=20, r=20, t=20, b=20))
    st.plotly_chart(fig, use_container_width=True)

    st.markdown(f"### `{g.spacer}` + `{g.pam}`")
    selected_validation = validate_opencrispr_guide(g, min_sequence_score=min_score,
                                                   min_specificity_review=min_spec_review,
                                                   panel_screen=screens.get(guide_key(record, g)))
    st.metric("Sequence and panel rule result", selected_validation.status)
    st.caption("Specificity status: " + selected_validation.specificity_status + ". PASS means the implemented checks passed; it does not mean proven editing or genome-wide specificity.")
    st.dataframe(pd.DataFrame([c.__dict__ for c in selected_validation.checks]), use_container_width=True, hide_index=True)
    st.markdown("#### Guide-expression formats reported with OpenCRISPR-1")
    st.dataframe(pd.DataFrame([x.__dict__ for x in guide_format_variants(g.spacer)]), use_container_width=True, hide_index=True)
    st.caption("GX19 has a natural first G; gX19 replaces base 1; gX20 adds a G before the intact spacer. X20 is the unmodified spacer. These are spacer formats, not full sgRNA constructs; provide an appropriate scaffold separately. Format order is not an efficacy ranking.")

    st.subheader("Exports")
    out_df = df.copy()
    out_df["App version"] = APP_VERSION
    out_df["Input source record/version"] = record.source_record_version
    out_df["Input assembly/genomic record"] = record.assembly
    out_df["Input sequence SHA-256"] = record.sequence_sha256
    out_df["Design settings"] = json.dumps(settings, sort_keys=True)
    fasta = "\n".join(f">guide_{i+1}|segment={x.segment_id}|start={x.start}|end={x.end}|strand={x.strand}|pam={x.pam}\n{x.spacer}" for i, x in enumerate(guides))
    b1, b2, b3 = st.columns(3)
    b1.download_button("Download CSV", out_df.to_csv(index=False).encode(), "opencrispr_guides.csv", "text/csv", use_container_width=True)
    b2.download_button("Download FASTA", fasta.encode(), "opencrispr_guides.fasta", "text/plain", use_container_width=True)
    b3.download_button("Download JSON", json.dumps(bundle, indent=2).encode(), "opencrispr_guides.json", "application/json", use_container_width=True)
    st.caption("JSON includes validation checks, expression formats, reference fingerprint, screening settings, total counts and retained hit details. Keep the original reference FASTA alongside the export.")

with st.expander("Scientific sources and remaining limits"):
    st.markdown("""
- [Ruffolo et al. 2025 — founding OpenCRISPR study](https://www.nature.com/articles/s41586-025-09298-z)
- [Hwang et al. 2026 — independent evaluation and guide formats](https://link.springer.com/article/10.1186/s13073-026-01682-2)
- [Tian et al. 2025 — comparative evaluation](https://doi.org/10.1126/sciadv.adu7334)

The heuristic is not a learned OpenCRISPR activity model. The local screen does not cover a whole genome, RNA/DNA bulges, non-NGG PAMs, chromatin or sample variants. Genomic mapping, scaffold selection, delivery, and experimental confirmation remain separate tasks. Base editing and prime editing are not designed by this app.
""")
