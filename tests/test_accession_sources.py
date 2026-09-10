from io import StringIO

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import accession_sources as ac


class FakeResponse:
    def __init__(self, *, json_data=None, text=""):
        self._json = json_data or {}
        self.text = text

    def json(self):
        return self._json


def _gb_text():
    seq = Seq("A" * 10 + "GACTGACTGACTGACTGACTTGG" + "C" * 50)
    rec = SeqRecord(seq, id="NM_TEST.2", name="TEST", description="mock accession")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["organism"] = "Homo sapiens"
    rec.annotations["date"] = "10-SEP-2026"
    rec.features = [
        SeqFeature(FeatureLocation(5, 70), type="CDS", qualifiers={"gene": ["TEST1"]}),
        SeqFeature(FeatureLocation(0, 40), type="exon"),
        SeqFeature(FeatureLocation(40, len(seq)), type="exon"),
    ]
    buf = StringIO()
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()


def test_fetch_ncbi_accession(monkeypatch):
    monkeypatch.setattr(ac, "_get", lambda *a, **k: FakeResponse(text=_gb_text()))
    rec = ac.fetch_ncbi_accession("NM_TEST.2")
    assert rec.gene == "TEST1"
    assert rec.organism == "Homo sapiens"
    assert rec.source == "NCBI accession"
    assert rec.accession == "NM_TEST.2"
    assert rec.segments
    assert len(rec.sequence_sha256) == 64


def test_fetch_ensembl_transcript_accession(monkeypatch):
    lookup = {
        "id": "ENST00000001",
        "object_type": "Transcript",
        "Parent": "ENSG00000001",
        "species": "homo_sapiens",
        "assembly_name": "GRCh38",
        "version": 4,
        "Exon": [{"id": "EX1"}],
    }

    def fake_get(url, *, params=None, headers=None, retries=3, timeout=30):
        if "/lookup/id/" in url:
            return FakeResponse(json_data=lookup)
        if url.endswith("/sequence/id/EX1"):
            return FakeResponse(text="A" * 10 + "GACTGACTGACTGACTGACTTGG" + "C" * 10)
        if url.endswith("/info/data"):
            return FakeResponse(json_data={"releases": [115]})
        raise AssertionError(url)

    monkeypatch.setattr(ac, "_get", fake_get)
    rec = ac.fetch_ensembl_accession("ENST00000001.4")
    assert rec.source == "Ensembl accession"
    assert rec.gene == "ENSG00000001"
    assert rec.assembly == "GRCh38"
    assert rec.annotation_release == "Ensembl release 115"
    assert rec.source_record_version.endswith(".4")


def test_fetch_accession_routes_and_rejects(monkeypatch):
    sentinel = object()
    monkeypatch.setattr(ac, "fetch_ncbi_accession", lambda accession: sentinel)
    monkeypatch.setattr(ac, "fetch_ensembl_accession", lambda accession: sentinel)
    assert ac.fetch_accession("NM_X", "NCBI") is sentinel
    assert ac.fetch_accession("ENSX", "Ensembl") is sentinel
    with pytest.raises(ValueError, match="NCBI or Ensembl"):
        ac.fetch_accession("X", "Other")
