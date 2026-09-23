# ============================================================================
# SEQUENCE SOURCES
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Loads, cleans and records sequence data before guide discovery begins.
#
# HOW TO READ THIS FILE:
# 1. Start with imports and constants to see the dependencies and fixed settings.
# 2. Read one function/class at a time rather than the whole file at once.
# 3. Follow the workflow from sequence input -> candidate discovery -> screening -> validation.
# 4. Scientific calculations, thresholds, validation rules and public APIs are
#    intentionally preserved while readability explanations are added.
#
# MAIN TOP-LEVEL PARTS:
# - class: GeneSequenceRecord
# - function: _ambiguity_fields
# - function: _append_ambiguity_warning
# - function: _ensembl_release
# - function: _ncbi_genomic_record
# - function: normalize_species_name
# - function: parse_multifasta
# - function: segments_from_ncbi_record
# - function: manual_record
# - function: _get
# - function: fetch_ncbi_gene
# - function: fetch_ensembl_gene
# - function: fetch_gene
# ============================================================================

#!/usr/bin/env python3
"""Gene/sequence retrieval for OpenCRISPR-1 guide design."""
from __future__ import annotations

from dataclasses import dataclass, field
from io import StringIO
from typing import Dict, List, Tuple
from datetime import datetime, timezone
import hashlib
import os
import re
import time
import requests
from Bio import SeqIO

from opencrispr_designer import clean_dna, ambiguity_summary

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "opencrispr.designer@example.com")
NCBI_TOOL = "OpenCRISPR_gRNA_Designer"
ENSEMBL = "https://rest.ensembl.org"

SPECIES_ALIASES = {
    "human": "homo_sapiens", "homo sapiens": "homo_sapiens",
    "mouse": "mus_musculus", "mus musculus": "mus_musculus",
    "rat": "rattus_norvegicus", "rattus norvegicus": "rattus_norvegicus",
    "zebrafish": "danio_rerio", "danio rerio": "danio_rerio",
    "arabidopsis": "arabidopsis_thaliana", "arabidopsis thaliana": "arabidopsis_thaliana",
    "rice": "oryza_sativa", "oryza sativa": "oryza_sativa",
    "maize": "zea_mays", "corn": "zea_mays", "zea mays": "zea_mays",
    "tomato": "solanum_lycopersicum", "solanum lycopersicum": "solanum_lycopersicum",
}


@dataclass
class GeneSequenceRecord:
    gene: str
    organism: str
    source: str
    accession: str
    description: str
    segments: List[Tuple[str, str]]
    warnings: List[str] = field(default_factory=list)
    retrieved_at_utc: str = field(default_factory=lambda: datetime.now(timezone.utc).replace(microsecond=0).isoformat())
    assembly: str = "unknown"
    annotation_release: str = "unknown"
    source_record_version: str = "unknown"
    sequence_sha256: str = ""
    ambiguity_count: int = 0
    ambiguity_codes: Tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not self.sequence_sha256:
            payload = "\n".join(f">{name}\n{seq}" for name, seq in self.segments).encode("utf-8")
            self.sequence_sha256 = hashlib.sha256(payload).hexdigest()

    @property
    def total_bp(self) -> int:
        return sum(len(s) for _, s in self.segments)

    def provenance_dict(self) -> Dict[str, object]:
        return {
            "gene": self.gene, "organism": self.organism, "source": self.source,
            "accession": self.accession, "source_record_version": self.source_record_version,
            "assembly_or_genomic_record": self.assembly, "annotation_release": self.annotation_release,
            "retrieved_at_utc": self.retrieved_at_utc, "sequence_sha256": self.sequence_sha256,
            "ambiguity_count": self.ambiguity_count, "ambiguity_codes": list(self.ambiguity_codes),
            "warnings": list(self.warnings),
        }



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _ambiguity_fields
# ----------------------------------------------------------------------------
def _ambiguity_fields(raw_seq: str) -> Tuple[int, Tuple[str, ...]]:
    a = ambiguity_summary(raw_seq)
    return int(a["count"]), tuple(a["codes"])



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _append_ambiguity_warning
# ----------------------------------------------------------------------------
def _append_ambiguity_warning(warnings: List[str], count: int, codes: Tuple[str, ...]) -> None:
    if count:
        warnings.append(
            f"Input contains {count} ambiguous IUPAC base(s) ({', '.join(codes)}). "
            "They are normalized to N; guide spacers containing ambiguity are skipped, "
            "and ambiguity is never silently treated as a mismatch. An unresolved base is "
            "accepted only at the degenerate N position of an otherwise resolved NGG/CCN PAM."
        )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _ensembl_release
# ----------------------------------------------------------------------------
def _ensembl_release(headers: Dict[str, str]) -> str:
    try:
        data = _get(f"{ENSEMBL}/info/data", headers=headers).json()
        releases = data.get("releases", [])
        return str(releases[0]) if releases else "unknown"
    except Exception:
        return "unknown"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _ncbi_genomic_record
# ----------------------------------------------------------------------------
def _ncbi_genomic_record(gene_id: str, common: Dict[str, str]) -> str:
    try:
        data = _get(
            f"{NCBI_EUTILS}/esummary.fcgi",
            params={**common, "db": "gene", "id": gene_id, "retmode": "json"},
        ).json()
        info = data.get("result", {}).get(str(gene_id), {})
        genomic = info.get("genomicinfo", []) or []
        return str(genomic[0].get("chraccver", "unknown")) if genomic else "unknown"
    except Exception:
        return "unknown"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: normalize_species_name
# ----------------------------------------------------------------------------
def normalize_species_name(organism: str) -> str:
    key = organism.strip().lower().replace("_", " ")
    return SPECIES_ALIASES.get(key, key.replace(" ", "_"))



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: parse_multifasta
# ----------------------------------------------------------------------------
def parse_multifasta(raw: str) -> Dict[str, str]:
    if not raw.strip():
        return {}
    if not any(ln.lstrip().startswith(">") for ln in raw.splitlines()):
        return {"sequence_1": clean_dna(raw)}
    if not raw.lstrip().startswith(">"):
        raise ValueError("FASTA sequence data must follow a header; remove text before the first >.")
    out = {}
    name = None
    chunks = []
    def save_record():
        if name is None:
            return
        seq = clean_dna("\n".join(chunks))
        if not seq:
            raise ValueError(f"FASTA record {name!r} is empty.")
        out[name] = seq
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            save_record()
            header = line[1:].strip()
            if not header:
                raise ValueError("Every FASTA record needs a non-empty identifier.")
            name = header.split()[0]
            if name in out:
                raise ValueError(f"Duplicate FASTA identifier {name!r}; give every record a unique name.")
            chunks = []
        else:
            chunks.append(line)
    save_record()
    if not out:
        raise ValueError("No readable FASTA records were found.")
    return out



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: segments_from_ncbi_record
# ----------------------------------------------------------------------------
def segments_from_ncbi_record(rec) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Extract independent annotated intervals, never reconstruct spliced DNA.

    Coordinates remain relative to the supplied record, in its stored orientation.
    A transcript without exon boundaries cannot establish genomic adjacency.
    """
    seq = clean_dna(str(rec.seq))
    if len(seq) > 500_000:
        raise ValueError("Use a gene-specific accession or a target region of at most 500,000 bases.")
    cds = [f for f in rec.features if f.type == "CDS"]
    exons = [f for f in rec.features if f.type == "exon"]
    if len(cds) > 1:
        raise ValueError("Multiple CDS annotations were found; provide a single-gene accession or explicitly selected genomic FASTA region.")
    is_rna = ("RNA" in str(rec.annotations.get("molecule_type", "")).upper()
              or rec.id.upper().startswith(("NM_", "XM_", "NR_", "XR_")))
    if is_rna and not exons:
        raise ValueError("Transcript exon boundaries are unavailable; use annotated exons or a reviewed genomic FASTA region. Spliced transcript fallback is disabled.")
    features = exons or cds
    warnings = ["Coordinates refer to independent source-record intervals, not chromosome positions. Verify assembly and coding context."]
    if not cds:
        warnings.append("No CDS feature was present; coding-region status is unverified.")
    if not features:
        if not seq:
            raise ValueError("The accession contains no usable sequence.")
        return [("accession_sequence", seq)], warnings
    intervals = set()
    for feature in features:
        if feature.location is None:
            raise ValueError("An annotation is missing sequence coordinates.")
        for part in feature.location.parts:
            if part.ref is not None or part.ref_db is not None:
                raise ValueError("Remote feature locations require an explicitly resolved target sequence.")
            a, b = int(part.start), int(part.end)
            # Intersect each exon with each CDS part, not the CDS bounding span.
            if exons and cds:
                if cds[0].location is None:
                    raise ValueError("The CDS annotation is missing sequence coordinates.")
                for c in cds[0].location.parts:
                    if c.ref is not None or c.ref_db is not None:
                        raise ValueError("Remote CDS locations require an explicitly resolved target sequence.")
                    x, y = max(a, int(c.start)), min(b, int(c.end))
                    if y > x:
                        intervals.add((x, y))
            else:
                intervals.add((a, b))
    prefix = "coding_exon" if cds and exons else "exon" if exons else "CDS_part"
    segments = [(f"{prefix}_{i}:{a+1}-{b}", seq[a:b])
                for i, (a, b) in enumerate(sorted(intervals), 1) if b - a >= 23]
    if not segments:
        raise ValueError("No independent annotated interval is at least 23 bases; short exons are not joined for DNA targeting.")
    return segments, warnings



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: manual_record
# ----------------------------------------------------------------------------
def manual_record(raw: str, gene: str = "manual_target", organism: str = "manual") -> GeneSequenceRecord:
    seqs = parse_multifasta(raw)
    if not seqs:
        raise ValueError("No sequence provided.")
    segments = [(name, seq) for name, seq in seqs.items()]
    a = ambiguity_summary(raw)
    count, codes = int(a["count"]), tuple(a["codes"])
    warnings = ["Manual sequence annotation is user supplied; verify intended genomic/exonic context, genotype, and genome assembly."]
    _append_ambiguity_warning(warnings, count, codes)
    return GeneSequenceRecord(
        gene=gene, organism=organism, source="Manual FASTA", accession=",".join(seqs.keys()),
        description="User-supplied sequence", segments=segments, warnings=warnings,
        assembly="user-supplied / unspecified", annotation_release="not applicable",
        source_record_version="user-supplied", ambiguity_count=count, ambiguity_codes=codes,
    )


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _get
# ----------------------------------------------------------------------------
def _get(url: str, *, params=None, headers=None, retries: int = 3, timeout: int = 30):
    err = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            r.raise_for_status()
            return r
        except Exception as exc:
            err = exc
            if attempt + 1 < retries:
                time.sleep(0.8 + attempt)
    raise RuntimeError(f"Network request failed: {err}")



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ncbi_gene
# ----------------------------------------------------------------------------
def fetch_ncbi_gene(gene: str, organism: str) -> GeneSequenceRecord:
    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
    q = f"{gene}[Gene Name] AND {organism}[Organism] AND alive[prop]"
    s = _get(f"{NCBI_EUTILS}/esearch.fcgi", params={**common, "db": "gene", "term": q, "retmax": 5, "retmode": "json"})
    ids = s.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError(f"NCBI Gene could not resolve '{gene}' in '{organism}'.")
    if len(ids) > 1:
        raise ValueError("Gene lookup is ambiguous; use a specific accession ID.")
    l = _get(f"{NCBI_EUTILS}/elink.fcgi", params={**common, "dbfrom": "gene", "db": "nuccore", "id": ids[0], "linkname": "gene_nuccore_refseqrna", "retmode": "json"})
    dbs = l.json().get("linksets", [{}])[0].get("linksetdbs", [])
    nids = dbs[0].get("links", []) if dbs else []
    if not nids:
        raise ValueError("No RefSeq RNA record linked to the selected NCBI gene.")
    f = _get(f"{NCBI_EUTILS}/efetch.fcgi", params={**common, "db": "nuccore", "id": ",".join(nids[:25]), "rettype": "gb", "retmode": "text"})
    records = list(SeqIO.parse(StringIO(f.text), "genbank"))
    if not records:
        raise ValueError("NCBI returned no readable nucleotide record.")

    def rank(rec):
        has_cds = any(feat.type == "CDS" for feat in rec.features)
        prefix = 0 if rec.id.upper().startswith("NM_") else 1 if rec.id.upper().startswith("XM_") else 2
        return (0 if has_cds else 1, prefix, -len(rec.seq))

    rec = sorted(records, key=rank)[0]
    segments, warnings = segments_from_ncbi_record(rec)
    warnings.append("Gene lookup selects one transcript using a heuristic; use Accession ID to control the isoform.")
    amb_count, amb_codes = _ambiguity_fields(str(rec.seq))
    _append_ambiguity_warning(warnings, amb_count, amb_codes)
    genomic_record = _ncbi_genomic_record(ids[0], common)
    record_date = str(rec.annotations.get("date", "unknown"))
    return GeneSequenceRecord(
        gene=gene, organism=organism, source="NCBI RefSeq", accession=rec.id,
        description=rec.description, segments=segments, warnings=warnings, assembly=genomic_record,
        annotation_release=f"RefSeq/GenBank record date {record_date}", source_record_version=rec.id,
        ambiguity_count=amb_count, ambiguity_codes=amb_codes,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ensembl_gene
# ----------------------------------------------------------------------------
def fetch_ensembl_gene(gene: str, organism: str) -> GeneSequenceRecord:
    species = normalize_species_name(organism)
    headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": NCBI_TOOL}
    if re.match(r"^ENS[A-Z0-9]*G\d+", gene, re.I):
        url = f"{ENSEMBL}/lookup/id/{gene}"
    else:
        url = f"{ENSEMBL}/lookup/symbol/{species}/{gene}"
    data = _get(url, params={"expand": 1}, headers=headers).json()
    txs = data.get("Transcript", [])
    if not txs:
        raise ValueError("Ensembl returned no expanded transcript annotation.")
    canonical = str(data.get("canonical_transcript", "")).split(".")[0]
    tx = next((t for t in txs if t.get("is_canonical") or str(t.get("id", "")).split(".")[0] == canonical), max(txs, key=lambda t: abs(int(t.get("end", 0)) - int(t.get("start", 0)))))
    exons = tx.get("Exon", [])
    if not exons:
        raise ValueError("Selected Ensembl transcript has no exon records.")
    segments = []
    ambiguity_count = 0
    ambiguity_codes_set = set()
    for i, exon in enumerate(exons, start=1):
        exon_id = exon.get("id")
        if not exon_id:
            continue
        raw_seq = _get(f"{ENSEMBL}/sequence/id/{exon_id}", headers={"Accept": "text/plain", "User-Agent": NCBI_TOOL}).text
        acount, acodes = _ambiguity_fields(raw_seq)
        ambiguity_count += acount
        ambiguity_codes_set.update(acodes)
        seq = clean_dna(raw_seq)
        if len(seq) >= 23:
            segments.append((f"exon_{i}:{exon_id}", seq))
    if not segments:
        raise ValueError("No exon sequence long enough for 20 nt + NGG discovery was retrieved.")
    warnings = ["Ensembl mode scans individual exons to avoid exon-junction guides; review whether the selected site lies in the intended coding region."]
    amb_codes = tuple(sorted(ambiguity_codes_set))
    _append_ambiguity_warning(warnings, ambiguity_count, amb_codes)
    release = _ensembl_release(headers)
    tx_id = str(tx.get("id", data.get("id", gene)))
    tx_version = tx.get("version")
    versioned_tx = f"{tx_id}.{tx_version}" if tx_version not in (None, "") else tx_id
    return GeneSequenceRecord(
        gene=gene, organism=organism, source="Ensembl REST", accession=versioned_tx,
        description=data.get("description", "") or gene, segments=segments, warnings=warnings,
        assembly=str(data.get("assembly_name", "unknown")), annotation_release=f"Ensembl release {release}",
        source_record_version=versioned_tx, ambiguity_count=ambiguity_count, ambiguity_codes=amb_codes,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_gene
# ----------------------------------------------------------------------------
def fetch_gene(gene: str, organism: str, source: str = "NCBI RefSeq") -> GeneSequenceRecord:
    if source.lower().startswith("ncbi"):
        return fetch_ncbi_gene(gene, organism)
    if source.lower().startswith("ensembl"):
        return fetch_ensembl_gene(gene, organism)
    raise ValueError("Unknown sequence source.")
