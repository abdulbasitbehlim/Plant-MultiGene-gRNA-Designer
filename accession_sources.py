#!/usr/bin/env python3
"""Direct accession-ID retrieval for the Plant MultiGene gRNA Designer.

This module adds a third input path without changing the existing gene-symbol/ID
lookup or manual FASTA workflows. NCBI accessions are fetched directly from
nuccore/RefSeq; Ensembl stable gene or transcript IDs are resolved through the
Ensembl REST API.
"""
from __future__ import annotations

from io import StringIO
from typing import List, Tuple
import re

from Bio import SeqIO

from plant_multiguide import clean_dna
from sequence_sources import (
    ENSEMBL,
    NCBI_EMAIL,
    NCBI_EUTILS,
    NCBI_TOOL,
    GeneSequenceRecord,
    _ambiguity_fields,
    _append_ambiguity_warning,
    _ensembl_release,
    _requests_get,
)


def _gene_label_from_genbank(rec, fallback: str) -> str:
    for feature_type in ("gene", "CDS"):
        for feature in rec.features:
            if feature.type == feature_type:
                values = feature.qualifiers.get("gene", [])
                if values:
                    return str(values[0])
    return fallback


def _segments_from_ncbi_record(rec) -> tuple[List[Tuple[str, str]], List[str]]:
    raw_seq = str(rec.seq)
    seq = clean_dna(raw_seq)
    warnings: List[str] = []
    cds_features = [f for f in rec.features if f.type == "CDS"]
    exon_features = [f for f in rec.features if f.type == "exon"]
    segments: List[Tuple[str, str]] = []

    if len(cds_features) > 1:
        warnings.append(
            "This accession contains multiple CDS annotations. The first CDS was selected; "
            "verify that the record represents the intended gene before experimental use."
        )

    if cds_features:
        cds = cds_features[0]
        c0, c1 = int(cds.location.start), int(cds.location.end)
        for idx, exon in enumerate(exon_features, start=1):
            e0, e1 = int(exon.location.start), int(exon.location.end)
            a, b = max(c0, e0), min(c1, e1)
            if b > a and b - a >= 23:
                segments.append((f"coding_exon_{idx}:{a+1}-{b}", seq[a:b]))
        if not segments:
            piece = seq[c0:c1]
            if piece:
                segments = [(f"spliced_CDS:{c0+1}-{c1}", piece)]
                warnings.append(
                    "Exon boundaries were unavailable for this accession; scanning uses the "
                    "spliced CDS, so genomic mapping should be checked for exon-junction guides."
                )
    else:
        if len(seq) > 500_000:
            raise ValueError(
                "This accession resolves to a large nucleotide record without a clear CDS. "
                "Use a gene/transcript-specific accession or reviewed FASTA instead."
            )
        segments = [("accession_sequence", seq)]
        warnings.append(
            "No CDS feature was found in this accession; the complete accession sequence is "
            "scanned and the intended coding/genomic context must be reviewed."
        )

    if not segments or not any(seq_piece for _, seq_piece in segments):
        raise ValueError("The accession did not contain a usable nucleotide sequence.")
    return segments, warnings


def fetch_ncbi_accession(accession: str) -> GeneSequenceRecord:
    requested = accession.strip()
    if not requested:
        raise ValueError("Enter an NCBI nucleotide/RefSeq accession ID.")

    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
    response = _requests_get(
        f"{NCBI_EUTILS}/efetch.fcgi",
        params={**common, "db": "nuccore", "id": requested, "rettype": "gb", "retmode": "text"},
    )
    records = list(SeqIO.parse(StringIO(response.text), "genbank"))
    if not records:
        raise ValueError(f"NCBI could not retrieve a readable nucleotide record for '{requested}'.")

    rec = records[0]
    segments, warnings = _segments_from_ncbi_record(rec)
    raw_seq = str(rec.seq)
    amb_count, amb_codes = _ambiguity_fields(raw_seq)
    _append_ambiguity_warning(warnings, amb_count, amb_codes)

    gene = _gene_label_from_genbank(rec, requested)
    organism = str(rec.annotations.get("organism", "unknown"))
    record_date = str(rec.annotations.get("date", "unknown"))
    warnings.insert(0, f"Sequence retrieved directly from NCBI using accession '{requested}'.")

    return GeneSequenceRecord(
        gene=gene,
        organism=organism,
        source="NCBI accession",
        accession=rec.id,
        description=rec.description,
        segments=segments,
        warnings=warnings,
        assembly="not resolved from direct accession lookup",
        annotation_release=f"RefSeq/GenBank record date {record_date}",
        source_record_version=rec.id,
        ambiguity_count=amb_count,
        ambiguity_codes=amb_codes,
    )


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
    gene_label = str(data.get("display_name") or data.get("id") or "unknown")
    return tx, gene_label


def fetch_ensembl_accession(accession: str) -> GeneSequenceRecord:
    requested = accession.strip()
    if not requested:
        raise ValueError("Enter an Ensembl stable gene or transcript ID.")

    stable_id = re.sub(r"\.\d+$", "", requested)
    headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": NCBI_TOOL}
    data = _requests_get(
        f"{ENSEMBL}/lookup/id/{stable_id}", params={"expand": 1}, headers=headers
    ).json()
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
        raw_seq = _requests_get(
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
        raise ValueError("No Ensembl exon sequence long enough for SpCas9 target discovery was retrieved.")

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
    description = str(data.get("description") or f"Ensembl record resolved from {requested}")

    return GeneSequenceRecord(
        gene=gene_label,
        organism=organism,
        source="Ensembl accession",
        accession=requested,
        description=description,
        segments=segments,
        warnings=warnings,
        assembly=str(data.get("assembly_name") or tx.get("assembly_name") or "unknown"),
        annotation_release=f"Ensembl release {_ensembl_release(headers)}",
        source_record_version=versioned_tx,
        ambiguity_count=ambiguity_count,
        ambiguity_codes=amb_codes,
    )


def fetch_accession(accession: str, source: str) -> GeneSequenceRecord:
    source_l = source.lower()
    if source_l.startswith("ncbi"):
        return fetch_ncbi_accession(accession)
    if source_l.startswith("ensembl"):
        return fetch_ensembl_accession(accession)
    raise ValueError("Accession source must be NCBI or Ensembl.")
