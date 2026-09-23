from sequence_sources import parse_multifasta, manual_record, normalize_species_name


def test_parse_multifasta_preserves_records():
    raw = ">ex1\nACGTACGT\n>ex2\nTTTTGGGG\n"
    out = parse_multifasta(raw)
    assert list(out) == ["ex1", "ex2"]
    assert out["ex1"] == "ACGTACGT"


def test_manual_record_segments_are_independent():
    rec = manual_record(">one\n" + "A" * 30 + "\n>two\n" + "C" * 30, gene="X")
    assert rec.gene == "X"
    assert len(rec.segments) == 2
    assert rec.total_bp == 60


def test_species_alias():
    assert normalize_species_name("Homo sapiens") == "homo_sapiens"
    assert normalize_species_name("rice") == "oryza_sativa"


def test_manual_provenance_and_ambiguity_fingerprint():
    from sequence_sources import manual_record
    r = manual_record(">seg1\nACGTRYNACGTACGTACGTACGTAGG\n", gene="X", organism="sample")
    assert r.ambiguity_count == 3
    assert set(r.ambiguity_codes) == {"N", "R", "Y"}
    assert r.assembly == "user-supplied / unspecified"
    assert len(r.sequence_sha256) == 64
    p = r.provenance_dict()
    assert p["retrieved_at_utc"]
    assert p["sequence_sha256"] == r.sequence_sha256
