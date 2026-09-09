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


def _ambiguity_fields(raw_seq: str) -> Tuple[int, Tuple[str, ...]]:
    a = ambiguity_summary(raw_seq)
    return int(a["count"]), tuple(a["codes"])


def _append_ambiguity_warning(warnings: List[str], count: int, codes: Tuple[str, ...]) -> None:
    if count:
        warnings.append(
            f"Input contains {count} ambiguous IUPAC base(s) ({', '.join(codes)}). "
            "They are normalized to N; guide spacers containing ambiguity are skipped, "
            "and ambiguity is never silently treated as a mismatch. An unresolved base is "
            "accepted only at the degenerate N position of an otherwise resolved NGG/CCN PAM."
        )


def _ensembl_release(headers: Dict[str, str]) -> str:
    try:
        data = _get(f"{ENSEMBL}/info/data", headers=headers).json()
        releases = data.get("releases", [])
        return str(releases[0]) if releases else "unknown"
    except Exception:
        return "unknown"


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


def normalize_species_name(organism: str) -> str:
    key = organism.strip().lower().replace("_", " ")
    return SPECIES_ALIASES.get(key, key.replace(" ", "_"))


def parse_multifasta(raw: str) -> Dict[str, str]:
    if not raw.strip():
        return {}
    if not any(ln.lstrip().startswith(">") for ln in raw.splitlines()):
        return {"sequence_1": clean_dna(raw)}
    out = {}
    for i, rec in enumerate(SeqIO.parse(StringIO(raw), "fasta"), start=1):
        out[rec.id or f"sequence_{i}"] = clean_dna(str(rec.seq))
    if not out:
        raise ValueError("No readable FASTA records were found.")
    return out


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


def fetch_ncbi_gene(gene: str, organism: str) -> GeneSequenceRecord:
    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
    q = f"{gene}[Gene Name] AND {organism}[Organism] AND alive[prop]"
    s = _get(f"{NCBI_EUTILS}/esearch.fcgi", params={**common, "db": "gene", "term": q, "retmax": 5, "retmode": "json"})
    ids = s.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError(f"NCBI Gene could not resolve '{gene}' in '{organism}'.")
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
    seq = clean_dna(str(rec.seq))
    cds_features = [x for x in rec.features if x.type == "CDS"]
    exon_features = [x for x in rec.features if x.type == "exon"]
    warnings: List[str] = []
    segments: List[Tuple[str, str]] = []
    if cds_features:
        cds = cds_features[0]
        c0, c1 = int(cds.location.start), int(cds.location.end)
        for idx, exon in enumerate(exon_features, start=1):
            e0, e1 = int(exon.location.start), int(exon.location.end)
            a, b = max(c0, e0), min(c1, e1)
            if b > a and b - a >= 23:
                segments.append((f"coding_exon_{idx}:{a+1}-{b}", seq[a:b]))
        if not segments:
            segments = [(f"spliced_CDS:{c0+1}-{c1}", seq[c0:c1])]
            warnings.append("Exon boundaries were unavailable in the chosen RefSeq RNA; review genomic mapping to exclude exon-junction candidates.")
    else:
        segments = [("transcript", seq)]
        warnings.append("No CDS feature was present; candidates are scanned from the transcript and require coding/genomic context review.")
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


def fetch_gene(gene: str, organism: str, source: str = "NCBI RefSeq") -> GeneSequenceRecord:
    if source.lower().startswith("ncbi"):
        return fetch_ncbi_gene(gene, organism)
    if source.lower().startswith("ensembl"):
        return fetch_ensembl_gene(gene, organism)
    raise ValueError("Unknown sequence source.")
