# ============================================================================
# ACCESSION SOURCES
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Retrieves and normalizes accession-based sequence sources used by the OpenCRISPR-1 workflow.
#
# HOW TO READ THIS FILE:
# 1. Start with imports and constants to see the dependencies and fixed settings.
# 2. Read one function/class at a time rather than the whole file at once.
# 3. Follow the workflow from sequence input -> candidate discovery -> screening -> validation.
# 4. Scientific calculations, thresholds, validation rules and public APIs are
#    intentionally preserved while readability explanations are added.
#
# MAIN TOP-LEVEL PARTS:
# - function: _gene_label_from_genbank
# - function: _segments_from_ncbi_record
# - function: fetch_ncbi_accession
# - function: _select_ensembl_transcript
# - function: fetch_ensembl_accession
# - function: fetch_accession
# ============================================================================

#!/usr/bin/env python3
"""Direct accession-ID retrieval for the OpenCRISPR-1 gRNA Designer.

Adds accession-based sequence retrieval without changing the existing gene lookup
or manual FASTA workflows. NCBI nucleotide/RefSeq accessions are fetched from
nuccore and Ensembl stable gene/transcript IDs are resolved through Ensembl REST.
"""
from __future__ import annotations

from io import StringIO
from typing import List, Tuple
import re

from Bio import SeqIO

from opencrispr_designer import clean_dna
from sequence_sources import (
    ENSEMBL,
    NCBI_EMAIL,
    NCBI_EUTILS,
    NCBI_TOOL,
    GeneSequenceRecord,
    _ambiguity_fields,
    _append_ambiguity_warning,
    _ensembl_release,
    _get,
    segments_from_ncbi_record,
)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _gene_label_from_genbank
# ----------------------------------------------------------------------------
def _gene_label_from_genbank(rec, fallback: str) -> str:
    for feature_type in ("gene", "CDS"):
        for feature in rec.features:
            if feature.type == feature_type:
                values = feature.qualifiers.get("gene", [])
                if values:
                    return str(values[0])
    return fallback



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _segments_from_ncbi_record
# ----------------------------------------------------------------------------
def _segments_from_ncbi_record(rec) -> tuple[List[Tuple[str, str]], List[str]]:
    return segments_from_ncbi_record(rec)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ncbi_accession
# ----------------------------------------------------------------------------
def fetch_ncbi_accession(accession: str) -> GeneSequenceRecord:
    requested = accession.strip()
    if not requested:
        raise ValueError("Enter an NCBI nucleotide/RefSeq accession ID.")

    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
    response = _get(
        f"{NCBI_EUTILS}/efetch.fcgi",
        params={**common, "db": "nuccore", "id": requested, "rettype": "gb", "retmode": "text"},
    )
    records = list(SeqIO.parse(StringIO(response.text), "genbank"))
    if not records:
        raise ValueError(f"NCBI could not retrieve a readable nucleotide record for '{requested}'.")

    if len(records) != 1:
        raise ValueError("Provide exactly one nucleotide accession.")
    rec = records[0]
    if "." in requested and rec.id.upper() != requested.upper():
        raise ValueError(f"Requested accession {requested} but received {rec.id}; version mismatch.")
    segments, warnings = _segments_from_ncbi_record(rec)
    amb_count, amb_codes = _ambiguity_fields(str(rec.seq))
    _append_ambiguity_warning(warnings, amb_count, amb_codes)
    warnings.insert(0, f"Sequence retrieved directly from NCBI using accession '{requested}'.")

    return GeneSequenceRecord(
        gene=_gene_label_from_genbank(rec, requested),
        organism=str(rec.annotations.get("organism", "unknown")),
        source="NCBI accession",
        accession=rec.id,
        description=rec.description,
        segments=segments,
        warnings=warnings,
        assembly="not resolved from direct accession lookup",
        annotation_release=f"RefSeq/GenBank record date {rec.annotations.get('date', 'unknown')}",
        source_record_version=rec.id,
        ambiguity_count=amb_count,
        ambiguity_codes=amb_codes,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _select_ensembl_transcript
# ----------------------------------------------------------------------------
def _select_ensembl_transcript(data: dict) -> tuple[dict, str]:
    object_type = str(data.get("object_type", "")).lower()
    if object_type == "transcript" or (data.get("Exon") and not data.get("Transcript")):
        return data, str(data.get("Parent") or data.get("display_name") or data.get("id") or "unknown")

    transcripts = data.get("Transcript", [])
    if not transcripts:
        raise ValueError("Ensembl returned no transcript/exon annotation for this stable ID.")
    canonical = str(data.get("canonical_transcript", "")).split(".")[0]
    tx = next(
        (t for t in transcripts if t.get("is_canonical") or str(t.get("id", "")).split(".")[0] == canonical),
        max(transcripts, key=lambda t: abs(int(t.get("end", 0)) - int(t.get("start", 0)))),
    )
    return tx, str(data.get("display_name") or data.get("id") or "unknown")



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ensembl_accession
# ----------------------------------------------------------------------------
def fetch_ensembl_accession(accession: str) -> GeneSequenceRecord:
    requested = accession.strip()
    if not requested:
        raise ValueError("Enter an Ensembl stable gene or transcript ID.")

    stable_id = re.sub(r"\.\d+$", "", requested)
    headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": NCBI_TOOL}
    data = _get(f"{ENSEMBL}/lookup/id/{stable_id}", params={"expand": 1}, headers=headers).json()
    requested_version = requested.split(".")[-1] if "." in requested else None
    if requested_version is not None and str(data.get("version")) != requested_version:
        raise ValueError(f"Requested Ensembl version {requested_version}, received {data.get('version', 'unknown')}; use a matching archived record or reviewed FASTA.")
    tx, gene_label = _select_ensembl_transcript(data)
    exons = tx.get("Exon", [])
    if not exons:
        raise ValueError("The resolved Ensembl transcript has no expanded exon records.")

    segments: List[Tuple[str, str]] = []
    ambiguity_count = 0
    ambiguity_codes_set = set()
    for idx, exon in enumerate(exons, start=1):
        exon_id = exon.get("id")
        if not exon_id:
            continue
        raw_seq = _get(
            f"{ENSEMBL}/sequence/id/{exon_id}",
            headers={"Accept": "text/plain", "User-Agent": NCBI_TOOL},
        ).text
        count, codes = _ambiguity_fields(raw_seq)
        ambiguity_count += count
        ambiguity_codes_set.update(codes)
        seq = clean_dna(raw_seq)
        if len(seq) >= 23:
            segments.append((f"exon_{idx}:{exon_id}", seq))

    if not segments:
        raise ValueError("No Ensembl exon sequence long enough for 20 nt + NGG discovery was retrieved.")

    warnings = [
        f"Sequence retrieved directly from Ensembl using stable ID '{requested}'.",
        "Accession mode scans individual transcript exons to avoid synthetic exon-junction guides; "
        "review coding-region context before experimental use.",
    ]
    amb_codes = tuple(sorted(ambiguity_codes_set))
    _append_ambiguity_warning(warnings, ambiguity_count, amb_codes)

    tx_id = str(tx.get("id", stable_id))
    tx_version = tx.get("version")
    versioned_tx = f"{tx_id}.{tx_version}" if tx_version not in (None, "") else tx_id
    organism = str(data.get("species") or tx.get("species") or "unknown").replace("_", " ")

    return GeneSequenceRecord(
        gene=gene_label,
        organism=organism,
        source="Ensembl accession",
        accession=requested,
        description=str(data.get("description") or f"Ensembl record resolved from {requested}"),
        segments=segments,
        warnings=warnings,
        assembly=str(data.get("assembly_name") or tx.get("assembly_name") or "unknown"),
        annotation_release=f"Ensembl release {_ensembl_release(headers)}",
        source_record_version=versioned_tx,
        ambiguity_count=ambiguity_count,
        ambiguity_codes=amb_codes,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_accession
# ----------------------------------------------------------------------------
def fetch_accession(accession: str, source: str) -> GeneSequenceRecord:
    source_l = source.lower()
    if source_l.startswith("ncbi"):
        return fetch_ncbi_accession(accession)
    if source_l.startswith("ensembl"):
        return fetch_ensembl_accession(accession)
    raise ValueError("Accession source must be NCBI or Ensembl.")
