#!/usr/bin/env python3
"""Core algorithms for Plant MultiGene gRNA Designer.

The module is intentionally UI-independent. It discovers SpCas9 (20 nt + NGG)
protospacers in one or more sequence segments, finds exact guide sequences shared
by every requested gene, and can additionally propose mismatch-aware consensus
guides for homologous gene families.

Mismatch-aware scores are transparent ranking proxies, not calibrated editing
probabilities and not a reimplementation of CRISPys/CFD.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple
import re

DNA = frozenset("ACGTN")
IUPAC_AMBIGUOUS = frozenset("NRYWSKMBDHVX")
MAX_GENES = 20
RECOMMENDED_MAX_GENES = 12
MAX_PAIR_COMPARISONS = 5_000_000


@dataclass(frozen=True)
class TargetSite:
    gene: str
    segment_id: str
    spacer: str
    pam: str
    strand: str
    start: int  # 1-based within segment, lower genomic/segment coordinate
    end: int    # inclusive
    gc_percent: float
    sequence_score: float


@dataclass(frozen=True)
class GeneMatch:
    gene: str
    target_spacer: str
    pam: str
    segment_id: str
    strand: str
    start: int
    end: int
    mismatches: int
    mismatch_positions: Tuple[int, ...]
    seed_mismatches: int
    compatibility_proxy: float


@dataclass
class SharedGuide:
    spacer: str
    design_type: str  # Exact shared or Mismatch-aware consensus
    matches: List[GeneMatch] = field(default_factory=list)
    gc_percent: float = 0.0
    sequence_score: float = 0.0
    average_compatibility: float = 0.0
    minimum_compatibility: float = 0.0
    exact_gene_count: int = 0
    worst_mismatch_count: int = 0
    notes: List[str] = field(default_factory=list)

    @property
    def genes_covered(self) -> int:
        return len({m.gene for m in self.matches})


@dataclass(frozen=True)
class PanelHit:
    contig: str
    spacer: str
    pam: str
    strand: str
    start: int
    mismatches: int
    seed_mismatches: int


@dataclass(frozen=True)
class SearchDiagnostics:
    genes: int
    total_pam_sites: int
    unique_seed_spacers: int
    upper_pair_comparisons: int
    strategy: str
    recommended_gene_limit: int = RECOMMENDED_MAX_GENES
    hard_gene_limit: int = MAX_GENES
    hard_pair_comparison_limit: int = MAX_PAIR_COMPARISONS


def clean_dna(raw: str) -> str:
    if not raw:
        return ""
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip() and not ln.lstrip().startswith(">")]
    seq = "".join(lines).upper().replace("U", "T")
    seq = re.sub(r"[\s\d]", "", seq)
    seq = re.sub(r"[RYSWKMBDHVX]", "N", seq)
    bad = sorted(set(seq) - DNA)
    if bad:
        raise ValueError("Unsupported sequence character(s): " + ", ".join(bad))
    return seq


def ambiguity_summary(raw: str) -> Dict[str, object]:
    """Describe IUPAC ambiguity without treating ambiguous symbols as mismatches.

    Ambiguity codes are normalized to N by :func:`clean_dna`. Candidate scanning
    skips any spacer containing N. An unresolved base is allowed only at the
    degenerate N position of an otherwise resolved NGG/CCN PAM; unresolved G/CC
    positions never create a PAM.
    """
    body = "".join(
        ln.strip() for ln in (raw or "").splitlines()
        if ln.strip() and not ln.lstrip().startswith(">")
    ).upper().replace("U", "T")
    body = re.sub(r"[\s\d]", "", body)
    codes = sorted({c for c in body if c in IUPAC_AMBIGUOUS})
    positions = tuple(i + 1 for i, c in enumerate(body) if c in IUPAC_AMBIGUOUS)
    return {"count": len(positions), "codes": codes, "positions": positions}


def reverse_complement(seq: str) -> str:
    return clean_dna(seq).translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def gc_percent(seq: str) -> float:
    seq = clean_dna(seq)
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def sequence_quality_score(spacer: str) -> float:
    """Transparent 0-100 sequence-quality ranking; not an activity probability."""
    s = clean_dna(spacer)
    if len(s) != 20 or "N" in s:
        return 0.0
    gc = gc_percent(s)
    score = 60.0
    if 40 <= gc <= 60:
        score += 20
    elif 30 <= gc <= 70:
        score += 10
    else:
        score -= 15
    if "TTTT" in s:
        score -= 20
    if re.search(r"A{5,}|C{5,}|G{5,}|T{5,}", s):
        score -= 12
    if s[-1] == "G":
        score += 4
    if s[-2] in "AG":
        score += 3
    return round(max(0.0, min(100.0, score)), 1)


def scan_spcas9(sequence: str, gene: str, segment_id: str = "segment") -> List[TargetSite]:
    """Enumerate both-strand 20 nt + NGG sites without crossing segment boundaries."""
    seq = clean_dna(sequence)
    sites: List[TargetSite] = []
    # Plus-strand PAM: [20nt][NGG]
    for pam_start in range(20, len(seq) - 2):
        pam = seq[pam_start:pam_start + 3]
        if len(pam) == 3 and pam[1:] == "GG":
            spacer = seq[pam_start - 20:pam_start]
            if "N" not in spacer:
                sites.append(TargetSite(
                    gene=gene, segment_id=segment_id, spacer=spacer, pam=pam,
                    strand="+", start=pam_start - 20 + 1, end=pam_start,
                    gc_percent=round(gc_percent(spacer), 1),
                    sequence_score=sequence_quality_score(spacer),
                ))
    # Reverse-strand PAM appears as CCN on the supplied strand: [CCN][20nt]
    for pam_start in range(0, len(seq) - 22):
        pam_raw = seq[pam_start:pam_start + 3]
        if len(pam_raw) == 3 and pam_raw[:2] == "CC":
            downstream = seq[pam_start + 3:pam_start + 23]
            if len(downstream) == 20 and "N" not in downstream:
                spacer = reverse_complement(downstream)
                pam = reverse_complement(pam_raw)
                sites.append(TargetSite(
                    gene=gene, segment_id=segment_id, spacer=spacer, pam=pam,
                    strand="-", start=pam_start + 4, end=pam_start + 23,
                    gc_percent=round(gc_percent(spacer), 1),
                    sequence_score=sequence_quality_score(spacer),
                ))
    return sites


def scan_gene_segments(gene_segments: Mapping[str, Sequence[Tuple[str, str]]]) -> Dict[str, List[TargetSite]]:
    """Scan independent exon/CDS/genomic segments; guides never span segment junctions."""
    out: Dict[str, List[TargetSite]] = {}
    for gene, segments in gene_segments.items():
        gene_sites: List[TargetSite] = []
        for segment_id, seq in segments:
            gene_sites.extend(scan_spcas9(seq, gene, segment_id))
        out[gene] = gene_sites
    return out


def mismatch_positions(guide: str, target: str) -> Tuple[int, ...]:
    g, t = clean_dna(guide), clean_dna(target)
    if len(g) != 20 or len(t) != 20:
        raise ValueError("Guide and target spacers must both be 20 nt.")
    return tuple(i + 1 for i, (a, b) in enumerate(zip(g, t)) if a != b)


def seed_mismatch_count(positions: Iterable[int]) -> int:
    # PAM-proximal 8 nt represented as guide positions 13-20.
    return sum(1 for p in positions if p >= 13)


def compatibility_proxy(guide: str, target: str) -> float:
    """Explainable mismatch-tolerance proxy; intentionally not called CFD or efficiency."""
    pos = mismatch_positions(guide, target)
    score = 100.0
    for p in pos:
        if p >= 13:
            score -= 22.0
        elif p >= 9:
            score -= 12.0
        else:
            score -= 6.0
    if len(pos) >= 4:
        score -= 8.0 * (len(pos) - 3)
    return round(max(0.0, score), 1)


def _weighted_distance(a: str, b: str) -> Tuple[int, int, int]:
    pos = mismatch_positions(a, b)
    seed = seed_mismatch_count(pos)
    weighted = sum(3 if p >= 13 else 2 if p >= 9 else 1 for p in pos)
    return weighted, len(pos), seed


def _to_match(gene: str, guide: str, site: TargetSite) -> GeneMatch:
    pos = mismatch_positions(guide, site.spacer)
    return GeneMatch(
        gene=gene, target_spacer=site.spacer, pam=site.pam, segment_id=site.segment_id,
        strand=site.strand, start=site.start, end=site.end,
        mismatches=len(pos), mismatch_positions=pos,
        seed_mismatches=seed_mismatch_count(pos),
        compatibility_proxy=compatibility_proxy(guide, site.spacer),
    )


def _build_shared(guide: str, design_type: str, matches: List[GeneMatch], notes: List[str] | None = None) -> SharedGuide:
    comps = [m.compatibility_proxy for m in matches]
    return SharedGuide(
        spacer=guide,
        design_type=design_type,
        matches=matches,
        gc_percent=round(gc_percent(guide), 1),
        sequence_score=sequence_quality_score(guide),
        average_compatibility=round(sum(comps) / len(comps), 1) if comps else 0.0,
        minimum_compatibility=min(comps) if comps else 0.0,
        exact_gene_count=sum(m.mismatches == 0 for m in matches),
        worst_mismatch_count=max((m.mismatches for m in matches), default=0),
        notes=notes or [],
    )


def exact_shared_guides(sites_by_gene: Mapping[str, Sequence[TargetSite]]) -> List[SharedGuide]:
    genes = list(sites_by_gene)
    if len(genes) < 2:
        return []
    per_gene: Dict[str, Dict[str, List[TargetSite]]] = {}
    for gene, sites in sites_by_gene.items():
        d: Dict[str, List[TargetSite]] = defaultdict(list)
        for s in sites:
            d[s.spacer].append(s)
        per_gene[gene] = d
    common = set(per_gene[genes[0]])
    for gene in genes[1:]:
        common &= set(per_gene[gene])
    results: List[SharedGuide] = []
    for spacer in common:
        matches = []
        for gene in genes:
            site = max(per_gene[gene][spacer], key=lambda x: x.sequence_score)
            matches.append(_to_match(gene, spacer, site))
        results.append(_build_shared(
            spacer, "Exact shared", matches,
            ["Same 20-nt spacer is PAM-compatible in every requested gene."],
        ))
    return sorted(results, key=_rank_key, reverse=True)


def _consensus(spacers: Sequence[str], preferred: str) -> str:
    cols = []
    for i in range(20):
        c = Counter(s[i] for s in spacers)
        top_n = max(c.values())
        tied = {b for b, n in c.items() if n == top_n}
        cols.append(preferred[i] if preferred[i] in tied else sorted(tied)[0])
    return "".join(cols)


def estimate_search_diagnostics(
    sites_by_gene: Mapping[str, Sequence[TargetSite]],
    max_pair_comparisons: int = MAX_PAIR_COMPARISONS,
) -> SearchDiagnostics:
    genes = len(sites_by_gene)
    total_sites = sum(len(v) for v in sites_by_gene.values())
    seeds = len({s.spacer for sites in sites_by_gene.values() for s in sites})
    # Two complete nearest-neighbour passes are the worst case: seed selection
    # and consensus re-optimization. No early stop is used for valid seeds.
    upper = 2 * seeds * total_sites
    return SearchDiagnostics(
        genes=genes, total_pam_sites=total_sites, unique_seed_spacers=seeds,
        upper_pair_comparisons=upper,
        strategy=(
            "Seed-exhaustive over every distinct observed PAM-compatible spacer; "
            "nearest PAM-compatible site selected per gene; majority consensus then "
            "re-optimized against every gene; all surviving unique proposals are ranked "
            "before max_results truncation. This is not exhaustive over all 4^20 synthetic spacers."
        ),
        hard_pair_comparison_limit=max_pair_comparisons,
    )


def _enforce_search_limits(
    sites_by_gene: Mapping[str, Sequence[TargetSite]],
    max_pair_comparisons: int = MAX_PAIR_COMPARISONS,
) -> SearchDiagnostics:
    d = estimate_search_diagnostics(sites_by_gene, max_pair_comparisons)
    if d.genes > MAX_GENES:
        raise ValueError(
            f"Mismatch-aware mode supports at most {MAX_GENES} genes per run. "
            f"Use multiple batches or a dedicated large-library workflow for {d.genes} genes."
        )
    if d.upper_pair_comparisons > max_pair_comparisons:
        raise ValueError(
            "Mismatch-aware search workload is too large for the interactive tool: "
            f"estimated worst-case {d.upper_pair_comparisons:,} pair comparisons exceeds "
            f"the {max_pair_comparisons:,} guard. Reduce genes/sequence length, use stricter "
            "input regions, or run a dedicated scalable multi-target workflow."
        )
    return d


def mismatch_aware_guides(
    sites_by_gene: Mapping[str, Sequence[TargetSite]],
    max_mismatches_per_gene: int = 2,
    max_seed_mismatches_per_gene: int = 1,
    min_compatibility: float = 55.0,
    max_results: int = 100,
    max_pair_comparisons: int = MAX_PAIR_COMPARISONS,
) -> List[SharedGuide]:
    """Find consensus guides with a PAM-compatible near-match in every gene.

    Search semantics are explicit: every distinct *observed* PAM-compatible spacer
    is evaluated as a seed (no first-hit/greedy stopping). For each seed, the nearest
    PAM-compatible target in every gene is selected, a majority consensus is formed,
    and the consensus is re-optimized against every gene. All surviving unique
    consensus proposals are globally ranked before ``max_results`` is applied.

    This is seed-exhaustive, not mathematically exhaustive over all possible 20-mers;
    therefore the top-ranked result is best among generated proposals, not a proof of
    the global optimum over the 4^20 sequence space.
    """
    genes = list(sites_by_gene)
    if len(genes) < 2 or any(not sites_by_gene[g] for g in genes):
        return []
    _enforce_search_limits(sites_by_gene, max_pair_comparisons=max_pair_comparisons)

    seed_spacers = list(dict.fromkeys(s.spacer for g in genes for s in sites_by_gene[g]))
    proposals: Dict[str, SharedGuide] = {}

    for seed in seed_spacers:
        selected: List[TargetSite] = []
        feasible = True
        for gene in genes:
            best = min(
                sites_by_gene[gene],
                key=lambda site: (_weighted_distance(seed, site.spacer), -site.sequence_score),
            )
            _, mm, seed_mm = _weighted_distance(seed, best.spacer)
            if mm > max_mismatches_per_gene or seed_mm > max_seed_mismatches_per_gene:
                feasible = False
                break
            selected.append(best)
        if not feasible:
            continue

        guide = _consensus([site.spacer for site in selected], preferred=seed)
        matches: List[GeneMatch] = []
        for gene in genes:
            best = min(
                sites_by_gene[gene],
                key=lambda site: (_weighted_distance(guide, site.spacer), -site.sequence_score),
            )
            match = _to_match(gene, guide, best)
            if (
                match.mismatches > max_mismatches_per_gene
                or match.seed_mismatches > max_seed_mismatches_per_gene
                or match.compatibility_proxy < min_compatibility
            ):
                feasible = False
                break
            matches.append(match)
        if not feasible:
            continue

        notes = [
            "Consensus guide has a PAM-compatible target in every requested gene.",
            "Search is exhaustive across observed PAM-compatible seed spacers, with no early stop; it is not exhaustive over all possible synthetic 20-mers.",
            "Mismatch compatibility is a transparent prioritization proxy, not a calibrated cleavage probability.",
        ]
        obj = _build_shared(guide, "Mismatch-aware consensus", matches, notes)
        previous = proposals.get(guide)
        if previous is None or _rank_key(obj) > _rank_key(previous):
            proposals[guide] = obj

    exact_spacers = {g.spacer for g in exact_shared_guides(sites_by_gene)}
    results = [g for spacer, g in proposals.items() if spacer not in exact_spacers]
    results.sort(key=_rank_key, reverse=True)
    return results[:max_results]

def _rank_key(g: SharedGuide) -> Tuple[float, ...]:
    gc_pref = -abs(g.gc_percent - 50.0)
    return (
        float(g.genes_covered),
        float(g.exact_gene_count),
        g.minimum_compatibility,
        g.average_compatibility,
        -float(g.worst_mismatch_count),
        g.sequence_score,
        gc_pref,
    )


def design_shared_guides(
    gene_segments: Mapping[str, Sequence[Tuple[str, str]]],
    include_mismatch_aware: bool = True,
    max_mismatches_per_gene: int = 2,
    max_seed_mismatches_per_gene: int = 1,
    min_compatibility: float = 55.0,
    max_results: int = 100,
    max_pair_comparisons: int = MAX_PAIR_COMPARISONS,
) -> Tuple[Dict[str, List[TargetSite]], List[SharedGuide]]:
    if len(gene_segments) < 2:
        raise ValueError("Enter at least two genes/sequences for a shared-guide design.")
    if len(gene_segments) > MAX_GENES:
        raise ValueError(f"This interactive tool supports at most {MAX_GENES} genes per run.")
    sites = scan_gene_segments(gene_segments)
    exact = exact_shared_guides(sites)
    consensus = []
    if include_mismatch_aware:
        consensus = mismatch_aware_guides(
            sites,
            max_mismatches_per_gene=max_mismatches_per_gene,
            max_seed_mismatches_per_gene=max_seed_mismatches_per_gene,
            min_compatibility=min_compatibility,
            max_results=max_results,
            max_pair_comparisons=max_pair_comparisons,
        )
    combined = exact + consensus
    combined.sort(key=lambda g: (g.design_type == "Exact shared",) + _rank_key(g), reverse=True)
    return sites, combined[:max_results]



def evaluate_candidate_guide(
    guide: str,
    sites_by_gene: Mapping[str, Sequence[TargetSite]],
) -> SharedGuide:
    """Evaluate an arbitrary 20-nt spacer against the currently scanned genes.

    The closest PAM-compatible site in each gene is retained even when it exceeds
    design thresholds; the validation layer can then explain exactly why the
    candidate passes, needs review, or fails.
    """
    spacer = clean_dna(guide)
    if len(spacer) != 20 or "N" in spacer:
        raise ValueError("Candidate guide must be exactly 20 resolved DNA bases (A/C/G/T).")
    matches: List[GeneMatch] = []
    missing = []
    for gene, sites in sites_by_gene.items():
        if not sites:
            missing.append(gene)
            continue
        best = min(
            sites,
            key=lambda site: (_weighted_distance(spacer, site.spacer), -site.sequence_score),
        )
        matches.append(_to_match(gene, spacer, best))
    exact = bool(matches) and not missing and all(m.mismatches == 0 for m in matches)
    notes = [
        "Custom validation candidate: closest PAM-compatible target retained for each gene.",
        "This evaluation does not create a new PAM site; it compares the supplied spacer with PAM-compatible sites already found in each gene.",
    ]
    if missing:
        notes.append("No PAM-compatible site was available for: " + ", ".join(missing))
    return _build_shared(
        spacer,
        "Exact shared" if exact else "Mismatch-aware consensus",
        matches,
        notes,
    )

def guide_summary_row(g: SharedGuide) -> Dict[str, object]:
    return {
        "Spacer (20 nt)": g.spacer,
        "Design type": g.design_type,
        "Genes covered": g.genes_covered,
        "Exact genes": g.exact_gene_count,
        "Worst mismatches": g.worst_mismatch_count,
        "Minimum compatibility": g.minimum_compatibility,
        "Average compatibility": g.average_compatibility,
        "GC%": g.gc_percent,
        "Sequence quality": g.sequence_score,
    }


def match_rows(g: SharedGuide) -> List[Dict[str, object]]:
    rows = []
    for m in g.matches:
        rows.append({
            "Gene": m.gene,
            "Target spacer": m.target_spacer,
            "PAM": m.pam,
            "Segment": m.segment_id,
            "Strand": m.strand,
            "Start": m.start,
            "End": m.end,
            "Mismatches": m.mismatches,
            "Mismatch positions": ",".join(map(str, m.mismatch_positions)) or "Exact",
            "Seed mismatches": m.seed_mismatches,
            "Compatibility proxy": m.compatibility_proxy,
        })
    return rows


def screen_reference_panel(
    guide: str,
    panel: Mapping[str, str],
    max_mismatches: int = 3,
    max_hits: int = 250,
) -> List[PanelHit]:
    """PAM-aware near-match screen for a supplied FASTA panel (not a whole-genome claim)."""
    hits: List[PanelHit] = []
    for contig, seq in panel.items():
        for site in scan_spcas9(seq, gene=contig, segment_id=contig):
            pos = mismatch_positions(guide, site.spacer)
            if len(pos) <= max_mismatches:
                hits.append(PanelHit(
                    contig=contig, spacer=site.spacer, pam=site.pam,
                    strand=site.strand, start=site.start,
                    mismatches=len(pos), seed_mismatches=seed_mismatch_count(pos),
                ))
    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches, h.contig, h.start))
    return hits[:max_hits]
