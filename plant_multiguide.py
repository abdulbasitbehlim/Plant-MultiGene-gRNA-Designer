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
import hashlib
import json
import re

DNA = frozenset("ACGTN")
IUPAC_AMBIGUOUS = frozenset("NRYWSKMBDHVX")
MAX_GENES = 20
RECOMMENDED_MAX_GENES = 12
MAX_PAIR_COMPARISONS = 5_000_000
MAX_INPUT_BP = 500_000


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
class PanelScreen:
    hits: Tuple[PanelHit, ...]
    total_hits: int
    scanned_sites: int
    panel_bp: int
    panel_sha256: str
    mismatch_radius: int
    hit_limit: int

    @property
    def truncated(self) -> bool:
        return self.total_hits > len(self.hits)

    def __iter__(self):
        return iter(self.hits)


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
    if sum(len(seq) for segments in gene_segments.values() for _, seq in segments) > MAX_INPUT_BP:
        raise ValueError(f"Input exceeds the {MAX_INPUT_BP:,} base interactive limit. Use reviewed target regions.")
    out: Dict[str, List[TargetSite]] = {}
    for gene, segments in gene_segments.items():
        if not gene.strip() or not segments:
            raise ValueError("Each gene needs a nonempty identifier and at least one segment.")
        if len({name for name, _ in segments}) != len(segments):
            raise ValueError(f"Duplicate segment identifiers in gene {gene}.")
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


def _check_thresholds(max_mm: int, max_seed_mm: int, min_compat: float, max_results: int = 1) -> None:
    if not isinstance(max_mm, int) or not 0 <= max_mm <= 20:
        raise ValueError("Mismatch limit must be an integer from 0 to 20.")
    if not isinstance(max_seed_mm, int) or not 0 <= max_seed_mm <= 8:
        raise ValueError("Seed mismatch limit must be an integer from 0 to 8.")
    if not 0 <= min_compat <= 100:
        raise ValueError("Minimum compatibility must be from 0 to 100.")
    if not isinstance(max_results, int) or max_results < 1:
        raise ValueError("max_results must be a positive integer.")


def _site_key(site: TargetSite) -> tuple:
    return (-site.sequence_score, site.segment_id, site.start, site.end, site.strand, site.pam)


def _select_match(guide: str, gene: str, sites: Sequence[TargetSite], max_mm: int,
                  max_seed_mm: int, min_compat: float, fallback: bool = False) -> GeneMatch | None:
    """Choose among feasible sites first; a closer invalid site must not hide one."""
    choices = [(_to_match(gene, guide, site), site) for site in sites]
    eligible = [(m, s) for m, s in choices if m.mismatches <= max_mm
                and m.seed_mismatches <= max_seed_mm and m.compatibility_proxy >= min_compat]
    pool = eligible or (choices if fallback else [])
    if not pool:
        return None
    return min(pool, key=lambda x: (-x[0].compatibility_proxy, x[0].mismatches,
                                    x[0].seed_mismatches, _site_key(x[1])))[0]


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
    genes = sorted(sites_by_gene)
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
    for spacer in sorted(common):
        matches = []
        for gene in genes:
            site = min(per_gene[gene][spacer], key=_site_key)
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
            "all constraints applied before selecting each per-gene site; observed feasible spacers retained; majority consensus then "
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
    feasible PAM-compatible target in every gene is selected. The observed spacer
    is retained, a majority consensus is formed, and the consensus is independently
    checked against every gene. Unique proposals are ranked before ``max_results``.

    This is seed-exhaustive, not mathematically exhaustive over all possible 20-mers;
    therefore the top-ranked result is best among generated proposals, not a proof of
    the global optimum over the 4^20 sequence space.
    """
    _check_thresholds(max_mismatches_per_gene, max_seed_mismatches_per_gene, min_compatibility, max_results)
    genes = sorted(sites_by_gene)
    if len(genes) < 2 or any(not sites_by_gene[g] for g in genes):
        return []
    _enforce_search_limits(sites_by_gene, max_pair_comparisons=max_pair_comparisons)
    seed_spacers = sorted({s.spacer for g in genes for s in sites_by_gene[g]})
    proposals: Dict[str, SharedGuide] = {}

    def matches_for(spacer):
        matches = [_select_match(spacer, gene, sites_by_gene[gene], max_mismatches_per_gene,
                                 max_seed_mismatches_per_gene, min_compatibility) for gene in genes]
        return matches if all(m is not None for m in matches) else []

    def retain(spacer, matches):
        if matches:
            proposals[spacer] = _build_shared(spacer, "Mismatch-aware consensus", matches, [
                "All target sites satisfy the selected mismatch, seed and compatibility constraints.",
                "Observed feasible spacers and seed-derived consensus proposals are retained; synthetic search is not exhaustive.",
                "Compatibility is an uncalibrated ranking proxy, not an editing probability.",
            ])

    for seed in seed_spacers:
        matches = matches_for(seed)
        if not matches:
            continue
        retain(seed, matches)
        consensus = _consensus([m.target_spacer for m in matches], preferred=seed)
        if consensus != seed and consensus not in proposals:
            retain(consensus, matches_for(consensus))
    exact_spacers = {g.spacer for g in exact_shared_guides(sites_by_gene)}
    results = [g for spacer, g in proposals.items() if spacer not in exact_spacers]
    results.sort(key=_rank_key, reverse=True)
    return results[:max_results]

def _rank_key(g: SharedGuide) -> tuple:
    gc_pref = -abs(g.gc_percent - 50.0)
    return (
        float(g.genes_covered),
        float(g.exact_gene_count),
        g.minimum_compatibility,
        g.average_compatibility,
        -float(g.worst_mismatch_count),
        g.sequence_score,
        gc_pref,
        g.spacer,  # deterministic final tie break
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
    _check_thresholds(max_mismatches_per_gene, max_seed_mismatches_per_gene, min_compatibility, max_results)
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
    max_mismatches_per_gene: int = 2,
    max_seed_mismatches_per_gene: int = 1,
    min_compatibility: float = 55.0,
) -> SharedGuide:
    """Evaluate an arbitrary 20-nt spacer against the currently scanned genes.

    The best proxy-scored site satisfying all thresholds is retained per gene.
    If none is feasible, the best proxy-scored site is retained for diagnosis;
    validation explains the threshold failures.
    """
    _check_thresholds(max_mismatches_per_gene, max_seed_mismatches_per_gene, min_compatibility)
    spacer = clean_dna(guide)
    if len(spacer) != 20 or "N" in spacer:
        raise ValueError("Candidate guide must be exactly 20 resolved DNA bases (A/C/G/T).")
    matches: List[GeneMatch] = []
    missing = []
    for gene, sites in sorted(sites_by_gene.items()):
        if not sites:
            missing.append(gene)
            continue
        matches.append(_select_match(spacer, gene, sites, max_mismatches_per_gene,
                                     max_seed_mismatches_per_gene, min_compatibility, fallback=True))
    exact = bool(matches) and not missing and all(m.mismatches == 0 for m in matches)
    notes = [
        "Custom validation candidate: best feasible PAM-compatible target retained per gene; if none is feasible, the best proxy match is shown for diagnosis.",
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


def screen_reference_panel_report(
    guide: str, panel: Mapping[str, str], max_mismatches: int = 3, max_hits: int = 250,
) -> PanelScreen:
    """Complete bounded local scan; capped display never hides total hit count."""
    _check_thresholds(max_mismatches, 8, 0, max_hits)
    guide = clean_dna(guide)
    if len(guide) != 20 or "N" in guide:
        raise ValueError("Panel guide must be 20 resolved DNA bases.")
    if not panel or any(not name.strip() or not clean_dna(seq) for name, seq in panel.items()):
        raise ValueError("Provide a nonempty reference panel with nonempty named sequences.")
    normalized = {name: clean_dna(seq) for name, seq in sorted(panel.items())}
    panel_bp = sum(map(len, normalized.values()))
    if panel_bp > MAX_INPUT_BP:
        raise ValueError(f"Panel exceeds {MAX_INPUT_BP:,} bases; use an external genome-wide workflow.")
    digest = hashlib.sha256(json.dumps(normalized, sort_keys=True).encode()).hexdigest()
    hits: List[PanelHit] = []
    scanned_sites = 0
    for contig, seq in normalized.items():
        sites = scan_spcas9(seq, gene=contig, segment_id=contig)
        scanned_sites += len(sites)
        for site in sites:
            pos = mismatch_positions(guide, site.spacer)
            if len(pos) <= max_mismatches:
                hits.append(PanelHit(contig, site.spacer, site.pam, site.strand, site.start,
                                     len(pos), seed_mismatch_count(pos)))
    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches, h.contig, h.start, h.strand))
    return PanelScreen(tuple(hits[:max_hits]), len(hits), scanned_sites, panel_bp,
                       digest, max_mismatches, max_hits)


def screen_reference_panel(guide: str, panel: Mapping[str, str], max_mismatches: int = 3,
                           max_hits: int = 250) -> List[PanelHit]:
    """Compatibility API; use the report API for total counts and provenance."""
    return list(screen_reference_panel_report(guide, panel, max_mismatches, max_hits).hits)


def suggest_exact_guide_set(sites_by_gene: Mapping[str, Sequence[TargetSite]], max_guides: int = 5) -> dict:
    """Deterministic greedy exact set cover; not a minimum-size or efficacy guarantee."""
    _check_thresholds(0, 0, 0, max_guides)
    by_spacer = defaultdict(dict)
    for gene, sites in sorted(sites_by_gene.items()):
        for site in sorted(sites, key=_site_key):
            by_spacer[site.spacer].setdefault(gene, site)
    uncovered = set(sites_by_gene)
    chosen = []
    while uncovered and by_spacer and len(chosen) < max_guides:
        spacer = min(by_spacer, key=lambda s: (-len(set(by_spacer[s]) & uncovered),
                                              -sequence_quality_score(s), s))
        covered = set(by_spacer[spacer])
        if not covered & uncovered:
            break
        matches = [_to_match(g, spacer, site) for g, site in sorted(by_spacer.pop(spacer).items())]
        chosen.append(_build_shared(spacer, "Exact guide set member", matches,
                                   ["Exact subset coverage in a greedy multi-guide fallback; not a single all-gene guide."]))
        uncovered -= covered
    return {"guides": chosen, "uncovered_genes": sorted(uncovered),
            "complete": not uncovered, "method": "greedy exact set cover; minimum size not guaranteed"}
