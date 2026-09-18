#!/usr/bin/env python3
"""Sequence retrieval helpers for plant gene-family guide design.

Network retrieval is deliberately separated from guide design so that the core can
be tested offline and users can always fall back to reviewed FASTA sequences.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from io import StringIO
from typing import Dict, List, Sequence, Tuple
from datetime import datetime, timezone
import hashlib
import os
import re
import time
import requests
from Bio import SeqIO

from plant_multiguide import clean_dna, ambiguity_summary, MAX_INPUT_BP

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "plant.multigene.grna@example.com")
NCBI_TOOL = "Plant_MultiGene_gRNA_Designer"
ENSEMBL = "https://rest.ensembl.org"

SPECIES_ALIASES = {
    "arabidopsis": "arabidopsis_thaliana",
    "arabidopsis thaliana": "arabidopsis_thaliana",
    "rice": "oryza_sativa",
    "oryza sativa": "oryza_sativa",
    "maize": "zea_mays",
    "corn": "zea_mays",
    "zea mays": "zea_mays",
    "tomato": "solanum_lycopersicum",
    "solanum lycopersicum": "solanum_lycopersicum",
    "soybean": "glycine_max",
    "glycine max": "glycine_max",
    "sorghum": "sorghum_bicolor",
    "sorghum bicolor": "sorghum_bicolor",
    "wheat": "triticum_aestivum",
    "triticum aestivum": "triticum_aestivum",
    "barley": "hordeum_vulgare",
    "hordeum vulgare": "hordeum_vulgare",
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
        return sum(len(seq) for _, seq in self.segments)

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
            "They are normalized to N; candidate spacers containing ambiguity are skipped, "
            "and ambiguity is never silently counted as a mismatch. An unresolved base is "
            "accepted only at the degenerate N position of an otherwise resolved NGG/CCN PAM."
        )


def _ensembl_release(headers: Dict[str, str]) -> str:
    try:
        data = _requests_get(f"{ENSEMBL}/info/data", headers=headers).json()
        releases = data.get("releases", [])
        return str(releases[0]) if releases else "unknown"
    except Exception:
        return "unknown"


def _ncbi_genomic_record(gene_id: str, common: Dict[str, str]) -> str:
    try:
        data = _requests_get(
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
    if len(raw) > 2 * MAX_INPUT_BP:
        raise ValueError("FASTA text exceeds the interactive input limit.")
    if not any(line.lstrip().startswith(">") for line in raw.splitlines()):
        return {"sequence_1": clean_dna(raw)}
    raw = "\n".join(line.strip() for line in raw.splitlines() if line.strip())
    if not raw.startswith(">"):
        raise ValueError("FASTA sequence text must follow a named header.")
    records: Dict[str, str] = {}
    for i, rec in enumerate(SeqIO.parse(StringIO(raw), "fasta"), start=1):
        name = rec.id or f"sequence_{i}"
        if name in records:
            raise ValueError(f"Duplicate FASTA identifier: {name}. Use unique headers.")
        if not rec.id or not str(rec.seq):
            raise ValueError("Every FASTA record needs a name and a nonempty sequence.")
        records[name] = clean_dna(str(rec.seq))
    if not records:
        raise ValueError("No readable FASTA records were found.")
    return records


def manual_records(raw: str, gene_names: Sequence[str] | None = None, organism: str = "manual",
                   group_segments: bool = False) -> Dict[str, GeneSequenceRecord]:
    raw = "\n".join(line.strip() for line in raw.splitlines())
    seqs = parse_multifasta(raw)
    raw_by_header: Dict[str, str] = {}
    if any(line.lstrip().startswith(">") for line in raw.splitlines()):
        for i, rec in enumerate(SeqIO.parse(StringIO(raw), "fasta"), start=1):
            raw_by_header[rec.id or f"sequence_{i}"] = str(rec.seq)
    else:
        raw_by_header["sequence_1"] = raw
    out: Dict[str, GeneSequenceRecord] = {}
    names = list(gene_names or [])
    pairs = list(seqs.items())
    if names and (len(names) != len(pairs) or len(set(names)) != len(names)):
        raise ValueError("Gene labels must be unique and match the number of FASTA records.")
    if names and len(pairs) == len(names):
        iterable = [(gene, header, seq) for gene, (header, seq) in zip(names, pairs)]
    else:
        iterable = [(header, header, seq) for header, seq in pairs]
    for gene, header, seq in iterable:
        if group_segments:
            if names or header.count("|") != 1 or not all(header.split("|")):
                raise ValueError("Grouped segments require GeneID|SegmentID headers and no replacement gene labels.")
            gene = header.split("|")[0]
        count, codes = _ambiguity_fields(raw_by_header.get(header, seq))
        warnings = ["Manual sequence boundaries are user supplied; confirm genomic/exonic context, intended genotype/cultivar, and assembly before experimental use."]
        _append_ambiguity_warning(warnings, count, codes)
        if gene in out:
            previous = out[gene]
            previous.segments.append((header, seq))
            previous.ambiguity_count += count
            previous.ambiguity_codes = tuple(sorted(set(previous.ambiguity_codes) | set(codes)))
            previous.warnings.extend(warnings[1:])
            previous.sequence_sha256 = ""
            previous.__post_init__()
            continue
        out[gene] = GeneSequenceRecord(
            gene=gene, organism=organism, source="Manual FASTA", accession=header,
            description=f"Manual sequence: {header}", segments=[(header, seq)], warnings=warnings,
            assembly="user-supplied / unspecified", annotation_release="not applicable",
            source_record_version="user-supplied", ambiguity_count=count, ambiguity_codes=codes,
        )
    return out


def extract_coding_segments(rec) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Keep CDS parts and exons independent, including pieces shorter than 23 nt."""
    seq = clean_dna(str(rec.seq))
    cds = [f for f in rec.features if f.type == "CDS"]
    exons = [f for f in rec.features if f.type == "exon"]
    if len(cds) > 1:
        raise ValueError("Multiple CDS annotations found. Use a gene-specific accession or reviewed segments.")
    if not seq or len(seq) > MAX_INPUT_BP:
        raise ValueError("Nucleotide record is empty or exceeds the interactive length limit.")
    if not cds:
        return [("accession_sequence", seq)], ["No CDS feature was found; full sequence requires coding and genomic context review."]
    parts = cds[0].location.parts
    intervals = []
    for part in parts:
        start, end = int(part.start), int(part.end)
        if exons:
            for exon in exons:
                for ep in exon.location.parts:
                    a, b = max(start, int(ep.start)), min(end, int(ep.end))
                    if b > a:
                        intervals.append((a, b))
        else:
            intervals.append((start, end))
    if not intervals:
        raise ValueError("CDS and exon annotations do not overlap; inspect the source record.")
    # Do not concatenate adjacent exons or fall back to a spliced CDS for short exons.
    prefix = "coding_exon" if exons else "coding_part" if len(parts) > 1 else "spliced_CDS"
    segments = [(f"{prefix}_{i}:{a+1}-{b}", seq[a:b])
                for i, (a, b) in enumerate(sorted(set(intervals)), 1)]
    warnings = []
    if not exons and len(parts) == 1:
        warnings.append("Exon boundaries were unavailable; spliced CDS targeting is provisional until genomic mapping excludes junction-spanning candidates.")
    return segments, warnings


def validate_record_set(records: Dict[str, GeneSequenceRecord]) -> None:
    """Reject known duplicate gene labels and mixed species, without guessing homology."""
    genes = [r.gene.casefold() for r in records.values()]
    if len(set(genes)) != len(genes):
        raise ValueError("Two records resolve to the same gene. Combine reviewed exons or choose distinct genes.")
    resolved = [(r.source.split()[0], r.source_record_version) for r in records.values()
                if r.source != "Manual FASTA" and r.source_record_version not in {"unknown", "user-supplied"}]
    if len(set(resolved)) != len(resolved):
        raise ValueError("Two inputs resolve to the same source record. Use distinct genes rather than aliases.")
    organisms = {normalize_species_name(r.organism) for r in records.values()
                 if r.organism.lower() not in {"unknown", "manual"}}
    if len(organisms) > 1:
        raise ValueError("Records resolve to different organisms. Use one plant organism per design.")

def _requests_get(url: str, *, params=None, headers=None, timeout: int = 30, retries: int = 3):
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
    query = f"{gene}[Gene Name] AND {organism}[Organism] AND alive[prop]"
    s = _requests_get(
        f"{NCBI_EUTILS}/esearch.fcgi",
        params={**common, "db": "gene", "term": query, "retmax": 5, "retmode": "json"},
    )
    ids = s.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError(f"NCBI Gene could not resolve '{gene}' in '{organism}'.")
    if len(ids) != 1:
        raise ValueError(f"NCBI found multiple genes for '{gene}'. Use a unique accession or reviewed FASTA.")

    link = _requests_get(
        f"{NCBI_EUTILS}/elink.fcgi",
        params={**common, "dbfrom": "gene", "db": "nuccore", "id": ids[0],
                "linkname": "gene_nuccore_refseqrna", "retmode": "json"},
    )
    dbs = link.json().get("linksets", [{}])[0].get("linksetdbs", [])
    nuccore_ids = dbs[0].get("links", []) if dbs else []
    if not nuccore_ids:
        raise ValueError(f"No RefSeq RNA record was linked to NCBI GeneID {ids[0]}.")

    gb = _requests_get(
        f"{NCBI_EUTILS}/efetch.fcgi",
        params={**common, "db": "nuccore", "id": ",".join(nuccore_ids[:25]),
                "rettype": "gb", "retmode": "text"},
    )
    records = list(SeqIO.parse(StringIO(gb.text), "genbank"))
    if not records:
        raise ValueError("NCBI returned no readable nucleotide records.")

    def record_rank(rec):
        has_cds = any(f.type == "CDS" for f in rec.features)
        prefix = 0 if rec.id.upper().startswith("NM_") else 1 if rec.id.upper().startswith("XM_") else 2
        return (0 if has_cds else 1, prefix, -len(rec.seq))

    rec = sorted(records, key=record_rank)[0]
    segments, warnings = extract_coding_segments(rec)
    amb_count, amb_codes = _ambiguity_fields(str(rec.seq))
    _append_ambiguity_warning(warnings, amb_count, amb_codes)
    genomic_record = _ncbi_genomic_record(ids[0], common)
    record_date = str(rec.annotations.get("date", "unknown"))
    return GeneSequenceRecord(
        gene=gene, organism=organism, source="NCBI RefSeq",
        accession=rec.id, description=rec.description, segments=segments, warnings=warnings,
        assembly=genomic_record, annotation_release=f"RefSeq/GenBank record date {record_date}",
        source_record_version=rec.id, ambiguity_count=amb_count, ambiguity_codes=amb_codes,
    )


def fetch_ensembl_gene(gene: str, organism: str) -> GeneSequenceRecord:
    species = normalize_species_name(organism)
    headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": NCBI_TOOL}
    if re.match(r"^ENS[A-Z0-9]*G\d+", gene, re.I):
        url = f"{ENSEMBL}/lookup/id/{gene}"
    else:
        url = f"{ENSEMBL}/lookup/symbol/{species}/{gene}"
    r = _requests_get(url, params={"expand": 1}, headers=headers)
    data = r.json()
    transcripts = data.get("Transcript", [])
    if not transcripts:
        raise ValueError("Ensembl returned the gene but no transcript/exon annotation.")
    canonical = str(data.get("canonical_transcript", "")).split(".")[0]
    tx = next(
        (t for t in transcripts if t.get("is_canonical") or str(t.get("id", "")).split(".")[0] == canonical),
        max(transcripts, key=lambda t: abs(int(t.get("end", 0)) - int(t.get("start", 0)))),
    )
    exons = tx.get("Exon", [])
    if not exons:
        raise ValueError("The selected Ensembl transcript has no expanded exon records.")
    segments: List[Tuple[str, str]] = []
    ambiguity_count = 0
    ambiguity_codes_set = set()
    for i, exon in enumerate(exons, start=1):
        exon_id = exon.get("id")
        if not exon_id:
            continue
        seq_r = _requests_get(
            f"{ENSEMBL}/sequence/id/{exon_id}",
            headers={"Accept": "text/plain", "User-Agent": NCBI_TOOL},
        )
        raw_seq = seq_r.text
        acount, acodes = _ambiguity_fields(raw_seq)
        ambiguity_count += acount
        ambiguity_codes_set.update(acodes)
        seq = clean_dna(raw_seq)
        if len(seq) >= 23:
            segments.append((f"exon_{i}:{exon_id}", seq))
    if not segments:
        raise ValueError("No Ensembl exon sequence long enough for SpCas9 target discovery was retrieved.")
    warnings = [
        "Ensembl mode scans individual transcript exons to avoid exon-junction guides. Exonic sequence is not automatically restricted to CDS; review coding context for knockout experiments."
    ]
    amb_codes = tuple(sorted(ambiguity_codes_set))
    _append_ambiguity_warning(warnings, ambiguity_count, amb_codes)
    release = _ensembl_release(headers)
    tx_id = str(tx.get("id", data.get("id", gene)))
    tx_version = tx.get("version")
    versioned_tx = f"{tx_id}.{tx_version}" if tx_version not in (None, "") else tx_id
    return GeneSequenceRecord(
        gene=gene, organism=organism, source="Ensembl REST", accession=versioned_tx,
        description=data.get("description", "") or f"Ensembl record for {gene}",
        segments=segments, warnings=warnings, assembly=str(data.get("assembly_name", "unknown")),
        annotation_release=f"Ensembl release {release}", source_record_version=versioned_tx,
        ambiguity_count=ambiguity_count, ambiguity_codes=amb_codes,
    )


def fetch_gene(gene: str, organism: str, source: str = "NCBI") -> GeneSequenceRecord:
    source_l = source.lower()
    if source_l.startswith("ncbi"):
        return fetch_ncbi_gene(gene, organism)
    if source_l.startswith("ensembl"):
        return fetch_ensembl_gene(gene, organism)
    raise ValueError("Source must be NCBI or Ensembl.")
