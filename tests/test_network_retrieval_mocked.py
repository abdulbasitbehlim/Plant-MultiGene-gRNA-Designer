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
    seq = Seq("A" * 10 + "CTCTACTTTCTCCCTCATCTTGG" + "C" * 47)
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

    def fake_get(url, *, params=None, headers=None, timeout=30, retries=3):
        if "esearch.fcgi" in url:
            return FakeResponse(json_data={"esearchresult": {"idlist": ["123"]}})
        if "elink.fcgi" in url:
            return FakeResponse(json_data={"linksets": [{"linksetdbs": [{"links": ["999"]}]}]})
        if "efetch.fcgi" in url:
            return FakeResponse(text=gb)
        if "esummary.fcgi" in url:
            return FakeResponse(json_data={"result": {"123": {"genomicinfo": [{"chraccver": "NC_003075.7"}]}}})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(ss, "_requests_get", fake_get)
    rec = ss.fetch_ncbi_gene("PUP7", "Arabidopsis thaliana")
    assert rec.source == "NCBI RefSeq"
    assert rec.accession == "NM_TEST.1"
    assert rec.source_record_version == "NM_TEST.1"
    assert rec.assembly == "NC_003075.7"
    assert "09-SEP-2026" in rec.annotation_release
    assert rec.segments and all("coding_exon" in name for name, _ in rec.segments)
    assert len(rec.sequence_sha256) == 64


def test_fetch_ncbi_gene_mocked_falls_back_without_exons(monkeypatch):
    gb = _genbank_text(with_exons=False)

    def fake_get(url, *, params=None, headers=None, timeout=30, retries=3):
        if "esearch.fcgi" in url:
            return FakeResponse(json_data={"esearchresult": {"idlist": ["123"]}})
        if "elink.fcgi" in url:
            return FakeResponse(json_data={"linksets": [{"linksetdbs": [{"links": ["999"]}]}]})
        if "efetch.fcgi" in url:
            return FakeResponse(text=gb)
        if "esummary.fcgi" in url:
            return FakeResponse(json_data={"result": {}})
        raise AssertionError(url)

    monkeypatch.setattr(ss, "_requests_get", fake_get)
    rec = ss.fetch_ncbi_gene("PUP7", "Arabidopsis thaliana")
    assert rec.segments[0][0].startswith("spliced_CDS")
    assert any("junction" in w.lower() for w in rec.warnings)
    assert rec.assembly == "unknown"


def test_fetch_ncbi_gene_error_branches(monkeypatch):
    monkeypatch.setattr(ss, "_requests_get", lambda *a, **k: FakeResponse(json_data={"esearchresult": {"idlist": []}}))
    with pytest.raises(ValueError, match="could not resolve"):
        ss.fetch_ncbi_gene("NOTREAL", "Arabidopsis thaliana")


def test_fetch_ensembl_gene_mocked_release_and_ambiguity(monkeypatch):
    lookup = {
        "id": "AT4G18197",
        "assembly_name": "TAIR10",
        "description": "mock PUP7",
        "canonical_transcript": "AT4G18197.1",
        "Transcript": [
            {
                "id": "AT4G18197.1",
                "version": 2,
                "is_canonical": 1,
                "start": 1,
                "end": 100,
                "Exon": [{"id": "EXON1"}, {"id": "EXON2"}, {}],
            }
        ],
    }

    def fake_get(url, *, params=None, headers=None, timeout=30, retries=3):
        if "/lookup/" in url:
            return FakeResponse(json_data=lookup)
        if url.endswith("/sequence/id/EXON1"):
            return FakeResponse(text="A" * 5 + "CTCTACTTTCTCCCTCATCTTGG" + "A" * 10)
        if url.endswith("/sequence/id/EXON2"):
            return FakeResponse(text="ACGTRYN" + "C" * 25)
        if url.endswith("/info/data"):
            return FakeResponse(json_data={"releases": [62]})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(ss, "_requests_get", fake_get)
    rec = ss.fetch_ensembl_gene("AT4G18197", "Arabidopsis thaliana")
    assert rec.assembly == "TAIR10"
    assert rec.annotation_release == "Ensembl release 62"
    assert rec.source_record_version.endswith(".2")
    assert len(rec.segments) == 2
    assert rec.ambiguity_count == 3
    assert set(rec.ambiguity_codes) == {"N", "R", "Y"}
    assert any("ambiguous" in w.lower() for w in rec.warnings)


def test_fetch_ensembl_gene_error_branches(monkeypatch):
    monkeypatch.setattr(ss, "_requests_get", lambda *a, **k: FakeResponse(json_data={"Transcript": []}))
    with pytest.raises(ValueError, match="no transcript"):
        ss.fetch_ensembl_gene("X", "rice")


def test_fetch_gene_routes_and_rejects_unknown(monkeypatch):
    sentinel = ss.GeneSequenceRecord("G", "O", "mock", "A.1", "d", [("s", "A" * 23)])
    monkeypatch.setattr(ss, "fetch_ncbi_gene", lambda gene, organism: sentinel)
    monkeypatch.setattr(ss, "fetch_ensembl_gene", lambda gene, organism: sentinel)
    assert ss.fetch_gene("G", "O", "NCBI") is sentinel
    assert ss.fetch_gene("G", "O", "Ensembl") is sentinel
    with pytest.raises(ValueError, match="NCBI or Ensembl"):
        ss.fetch_gene("G", "O", "Other")


def test_requests_get_retry_success_and_failure(monkeypatch):
    calls = {"n": 0}

    def flaky(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] < 2:
            raise RuntimeError("temporary")
        return FakeResponse(text="ok")

    monkeypatch.setattr(ss.requests, "get", flaky)
    monkeypatch.setattr(ss.time, "sleep", lambda *_: None)
    assert ss._requests_get("https://example.test", retries=2).text == "ok"
    assert calls["n"] == 2

    monkeypatch.setattr(ss.requests, "get", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("down")))
    with pytest.raises(RuntimeError, match="Network request failed"):
        ss._requests_get("https://example.test", retries=1)
