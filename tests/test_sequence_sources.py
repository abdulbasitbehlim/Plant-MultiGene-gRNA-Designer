from sequence_sources import parse_multifasta, manual_records, normalize_species_name


def test_parse_multifasta_and_manual_records():
    raw = ">GeneA\nACGTACGT\n>GeneB description\nTTTTCCCC\n"
    seqs = parse_multifasta(raw)
    assert seqs["GeneA"] == "ACGTACGT"
    assert seqs["GeneB"] == "TTTTCCCC"
    recs = manual_records(raw, organism="Arabidopsis thaliana")
    assert set(recs) == {"GeneA", "GeneB"}
    assert recs["GeneA"].segments[0][1] == "ACGTACGT"


def test_species_aliases():
    assert normalize_species_name("rice") == "oryza_sativa"
    assert normalize_species_name("Arabidopsis thaliana") == "arabidopsis_thaliana"


def test_manual_provenance_and_ambiguity_fingerprint():
    from sequence_sources import manual_records
    raw = ">GeneA\nACGTRYNACGT\n>GeneB\nTTTTCCCCAAAA\n"
    recs = manual_records(raw, organism="test cultivar")
    a = recs["GeneA"]
    assert a.ambiguity_count == 3
    assert set(a.ambiguity_codes) == {"N", "R", "Y"}
    assert a.assembly == "user-supplied / unspecified"
    assert a.source_record_version == "user-supplied"
    assert len(a.sequence_sha256) == 64
    assert a.provenance_dict()["sequence_sha256"] == a.sequence_sha256
