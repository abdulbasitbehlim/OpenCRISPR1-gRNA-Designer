from io import StringIO

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import sequence_sources as ss


class FakeResponse:
    def __init__(self, *, json_data=None, text=""):
        self._json = json_data or {}
        self.text = text

    def json(self):
        return self._json

    def raise_for_status(self):
        return None


def _genbank_text(*, with_exons=True, with_cds=True, record_id="NM_TEST.1"):
    seq = Seq("A" * 10 + "GAGTCCGAGCAGAAGAAGAATGG" + "C" * 47)
    rec = SeqRecord(seq, id=record_id, name="TEST", description="mock RefSeq transcript")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["date"] = "09-SEP-2026"
    features = []
    if with_cds:
        features.append(SeqFeature(FeatureLocation(5, 65), type="CDS"))
    if with_exons:
        features.extend([
            SeqFeature(FeatureLocation(0, 40), type="exon"),
            SeqFeature(FeatureLocation(40, len(seq)), type="exon"),
        ])
    rec.features = features
    buf = StringIO()
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()


def test_fetch_ncbi_gene_mocked_preserves_version_assembly_and_exons(monkeypatch):
    gb = _genbank_text()

    def fake_get(url, *, params=None, headers=None, retries=3, timeout=30):
        if "esearch.fcgi" in url:
            return FakeResponse(json_data={"esearchresult": {"idlist": ["7157"]}})
        if "elink.fcgi" in url:
            return FakeResponse(json_data={"linksets": [{"linksetdbs": [{"links": ["999"]}]}]})
        if "efetch.fcgi" in url:
            return FakeResponse(text=gb)
        if "esummary.fcgi" in url:
            return FakeResponse(json_data={"result": {"7157": {"genomicinfo": [{"chraccver": "NC_000017.11"}]}}})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(ss, "_get", fake_get)
    rec = ss.fetch_ncbi_gene("TP53", "Homo sapiens")
    assert rec.accession == "NM_TEST.1"
    assert rec.source_record_version == "NM_TEST.1"
    assert rec.assembly == "NC_000017.11"
    assert "09-SEP-2026" in rec.annotation_release
    assert rec.segments and all("coding_exon" in n for n, _ in rec.segments)
    assert len(rec.sequence_sha256) == 64


def test_fetch_ncbi_gene_mocked_no_cds_uses_transcript(monkeypatch):
    gb = _genbank_text(with_exons=False, with_cds=False)

    def fake_get(url, *, params=None, headers=None, retries=3, timeout=30):
        if "esearch.fcgi" in url:
            return FakeResponse(json_data={"esearchresult": {"idlist": ["1"]}})
        if "elink.fcgi" in url:
            return FakeResponse(json_data={"linksets": [{"linksetdbs": [{"links": ["2"]}]}]})
        if "efetch.fcgi" in url:
            return FakeResponse(text=gb)
        if "esummary.fcgi" in url:
            return FakeResponse(json_data={"result": {}})
        raise AssertionError(url)

    monkeypatch.setattr(ss, "_get", fake_get)
    rec = ss.fetch_ncbi_gene("X", "Homo sapiens")
    assert rec.segments[0][0] == "transcript"
    assert any("no cds" in w.lower() for w in rec.warnings)


def test_fetch_ncbi_gene_error_branches(monkeypatch):
    monkeypatch.setattr(ss, "_get", lambda *a, **k: FakeResponse(json_data={"esearchresult": {"idlist": []}}))
    with pytest.raises(ValueError, match="could not resolve"):
        ss.fetch_ncbi_gene("NOTREAL", "Homo sapiens")


def test_fetch_ensembl_gene_mocked_release_and_ambiguity(monkeypatch):
    lookup = {
        "id": "ENSG00000141510",
        "assembly_name": "GRCh38",
        "description": "mock TP53",
        "canonical_transcript": "ENST00000269305",
        "Transcript": [
            {
                "id": "ENST00000269305",
                "version": 9,
                "is_canonical": 1,
                "start": 1,
                "end": 100,
                "Exon": [{"id": "ENSE0001"}, {"id": "ENSE0002"}, {}],
            }
        ],
    }

    def fake_get(url, *, params=None, headers=None, retries=3, timeout=30):
        if "/lookup/" in url:
            return FakeResponse(json_data=lookup)
        if url.endswith("/sequence/id/ENSE0001"):
            return FakeResponse(text="A" * 5 + "GAGTCCGAGCAGAAGAAGAATGG" + "A" * 10)
        if url.endswith("/sequence/id/ENSE0002"):
            return FakeResponse(text="ACGTRYN" + "C" * 25)
        if url.endswith("/info/data"):
            return FakeResponse(json_data={"releases": [115]})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(ss, "_get", fake_get)
    rec = ss.fetch_ensembl_gene("ENSG00000141510", "Homo sapiens")
    assert rec.assembly == "GRCh38"
    assert rec.annotation_release == "Ensembl release 115"
    assert rec.source_record_version.endswith(".9")
    assert rec.ambiguity_count == 3
    assert set(rec.ambiguity_codes) == {"N", "R", "Y"}
    assert len(rec.segments) == 2


def test_fetch_ensembl_gene_error_branches(monkeypatch):
    monkeypatch.setattr(ss, "_get", lambda *a, **k: FakeResponse(json_data={"Transcript": []}))
    with pytest.raises(ValueError, match="no expanded transcript"):
        ss.fetch_ensembl_gene("X", "human")


def test_fetch_gene_routes_and_rejects_unknown(monkeypatch):
    sentinel = ss.GeneSequenceRecord("G", "O", "mock", "A.1", "d", [("s", "A" * 23)])
    monkeypatch.setattr(ss, "fetch_ncbi_gene", lambda gene, organism: sentinel)
    monkeypatch.setattr(ss, "fetch_ensembl_gene", lambda gene, organism: sentinel)
    assert ss.fetch_gene("G", "O", "NCBI RefSeq") is sentinel
    assert ss.fetch_gene("G", "O", "Ensembl") is sentinel
    with pytest.raises(ValueError, match="Unknown sequence source"):
        ss.fetch_gene("G", "O", "Other")


def test_get_retry_success_and_failure(monkeypatch):
    calls = {"n": 0}

    def flaky(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] < 2:
            raise RuntimeError("temporary")
        return FakeResponse(text="ok")

    monkeypatch.setattr(ss.requests, "get", flaky)
    monkeypatch.setattr(ss.time, "sleep", lambda *_: None)
    assert ss._get("https://example.test", retries=2).text == "ok"
    assert calls["n"] == 2

    monkeypatch.setattr(ss.requests, "get", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("down")))
    with pytest.raises(RuntimeError, match="Network request failed"):
        ss._get("https://example.test", retries=1)
