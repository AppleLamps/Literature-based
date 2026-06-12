"""ABC discovery engine.

Implements Swanson's open-discovery model over a MeSH co-occurrence substrate:

    A  --(co-occurs in A's literature)-->  B  --(co-occurs in B's literature)-->  C

where C is *not* directly linked to A. Each surviving C is a candidate hypothesis
A--C bridged by one or more intermediates B.

The engine is substrate-agnostic in spirit: edges carry a co-occurrence weight
and an optional predicate/direction. With the live MeSH substrate predicates are
unavailable, so direction is left as ``None`` and the directional-coherence bonus
defaults off (documented in DATA.md / DECISIONS.md).
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

from .eutils import EutilsClient


@dataclass
class ALiterature:
    a_query: str
    a_name: str
    corpus_size: int
    fetched_pmids: List[str]
    # B-term -> co-occurrence count with A (number of A-articles carrying B)
    b_counts: Dict[str, int]
    # B-term -> set of representative PMIDs for the A-B edge
    b_pmids: Dict[str, List[str]]
    # every MeSH descriptor that co-occurs with A in the fetched corpus
    a_cooccurring: Set[str]


@dataclass
class Bridge:
    b_term: str
    ab_count: int
    bc_count: int
    ab_pmids: List[str]
    bc_pmids: List[str]
    b_productivity: int = 0  # number of distinct C-candidates this B yields (for IDF)
    b_tree: str = ""         # primary MeSH tree letter of B (semantic-type weight)
    direction: Optional[str] = None  # reserved for predicate-aware substrates


@dataclass
class Candidate:
    c_term: str
    bridges: List[Bridge] = field(default_factory=list)
    score: float = 0.0
    flags: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_bridges(self) -> int:
        return len(self.bridges)

    @property
    def b_terms(self) -> List[str]:
        return [b.b_term for b in self.bridges]


def _stop_set(cfg: Dict[str, Any]) -> Set[str]:
    return {d.lower() for d in cfg["mesh_stoplist"]["descriptors"]}


def _tally_mesh(
    records: Dict[str, Dict[str, Any]],
    stop: Set[str],
    only_major: bool,
    exclude: Set[str],
) -> Tuple[Dict[str, int], Dict[str, List[str]]]:
    """Tally descriptor -> document count, plus example PMIDs per descriptor."""
    counts: Dict[str, int] = defaultdict(int)
    pmids: Dict[str, List[str]] = defaultdict(list)
    exclude_lc = {e.lower() for e in exclude}
    for pmid, rec in records.items():
        seen: Set[str] = set()
        for m in rec.get("mesh", []):
            name = m["name"]
            lc = name.lower()
            if lc in stop or lc in exclude_lc:
                continue
            if only_major and not m.get("major"):
                continue
            if name in seen:
                continue
            seen.add(name)
            counts[name] += 1
            if len(pmids[name]) < 8:
                pmids[name].append(pmid)
    return dict(counts), dict(pmids)


def gather_a_literature(
    client: EutilsClient,
    cfg: Dict[str, Any],
    a_query: str,
    a_name: str,
    *,
    maxdate: Optional[str] = None,
    mindate: Optional[str] = None,
) -> ALiterature:
    """Phase 4.1: gather A's literature and the B-terms linked to A."""
    corpus_cfg = cfg["corpus"]
    field = corpus_cfg["a_literature_query_field"]
    term = f'"{a_query}"{field}' if field else a_query
    total, pmids = client.esearch(
        term, retmax=corpus_cfg["a_max_pmids"], mindate=mindate, maxdate=maxdate
    )
    records = client.efetch_mesh(pmids)
    stop = _stop_set(cfg)
    counts, pmid_map = _tally_mesh(
        records, stop, corpus_cfg["only_major_topic"], exclude={a_name, a_query}
    )
    a_cooccurring = set(counts.keys())
    # Keep only B-terms above the minimum A-B co-occurrence threshold.
    min_ab = corpus_cfg["min_ab_cooccurrence"]
    b_counts = {b: c for b, c in counts.items() if c >= min_ab}
    return ALiterature(
        a_query=a_query,
        a_name=a_name,
        corpus_size=total,
        fetched_pmids=pmids,
        b_counts=b_counts,
        b_pmids={b: pmid_map.get(b, []) for b in b_counts},
        a_cooccurring=a_cooccurring,
    )


def select_b_terms(a_lit: ALiterature, cfg: Dict[str, Any]) -> List[str]:
    top_n = cfg["corpus"]["b_top_n"]
    ranked = sorted(a_lit.b_counts.items(), key=lambda kv: kv[1], reverse=True)
    return [b for b, _ in ranked[:top_n]]


def expand_b_and_collect_c(
    client: EutilsClient,
    cfg: Dict[str, Any],
    a_lit: ALiterature,
    b_terms: List[str],
    *,
    maxdate: Optional[str] = None,
    mindate: Optional[str] = None,
    mesh_tree=None,
    progress=None,
) -> List[Candidate]:
    """Phases 4.2-4.3: expand each B to its C-terms and collect C not linked to A."""
    corpus_cfg = cfg["corpus"]
    stop = _stop_set(cfg)
    field = corpus_cfg["a_literature_query_field"]
    min_bc = corpus_cfg["min_bc_cooccurrence"]

    # Concepts that are directly linked to A are not hidden -> excluded as C.
    directly_linked = {x.lower() for x in a_lit.a_cooccurring}
    directly_linked.add(a_lit.a_name.lower())
    directly_linked.add(a_lit.a_query.lower())

    # Semantic-type filter: keep only C-candidates in the allowed MeSH branches
    # (default tree D = Chemicals & Drugs, i.e. candidate substances/therapies).
    tree_filter_on = bool(cfg["corpus"].get("candidate_tree_filter")) and mesh_tree is not None \
        and getattr(mesh_tree, "available", False)
    whitelist = cfg["corpus"].get("candidate_tree_whitelist", ["D"])

    def c_allowed(name: str) -> bool:
        if not tree_filter_on:
            return True
        return mesh_tree.in_branches(name, whitelist, unknown_ok=False)

    # c_term -> list of Bridge
    candidates: Dict[str, List[Bridge]] = defaultdict(list)
    b_productivity: Dict[str, int] = {}

    for idx, b in enumerate(b_terms):
        if progress:
            progress(idx + 1, len(b_terms), b)
        term = f'"{b}"{field}' if field else b
        _, b_pmids = client.esearch(
            term, retmax=corpus_cfg["b_max_pmids"], mindate=mindate, maxdate=maxdate
        )
        if not b_pmids:
            continue
        records = client.efetch_mesh(b_pmids)
        c_counts, c_pmid_map = _tally_mesh(records, stop, corpus_cfg["only_major_topic"], exclude={b})
        kept = 0
        for c_term, bc_count in c_counts.items():
            if bc_count < min_bc:
                continue
            if c_term.lower() in directly_linked:
                continue  # already connected to A -> not a hidden link
            if not c_allowed(c_term):
                continue  # wrong semantic type (not a candidate substance/therapy)
            kept += 1
            bridge = Bridge(
                b_term=b,
                ab_count=a_lit.b_counts.get(b, 0),
                bc_count=bc_count,
                ab_pmids=a_lit.b_pmids.get(b, [])[: cfg["output"]["pmids_per_edge"]],
                bc_pmids=c_pmid_map.get(c_term, [])[: cfg["output"]["pmids_per_edge"]],
            )
            candidates[c_term].append(bridge)
        b_productivity[b] = kept

    # A B-term that bridges to few C's is a specific, informative intermediate; a
    # hub B that links to almost everything carries little signal. Record each
    # bridge's parent-B productivity so scoring can apply an IDF-style weight, and
    # the B-term's primary MeSH tree so scoring can prefer mechanistic
    # intermediates (physiological processes) over disease co-membership.
    def b_tree_letter(name: str) -> str:
        if mesh_tree is None or not getattr(mesh_tree, "available", False):
            return ""
        t = mesh_tree.trees(name)
        return t[0][0] if t else ""

    tree_cache = {b: b_tree_letter(b) for b in b_terms}
    for brs in candidates.values():
        for br in brs:
            br.b_productivity = b_productivity.get(br.b_term, 0)
            br.b_tree = tree_cache.get(br.b_term, "")

    return [Candidate(c_term=c, bridges=brs) for c, brs in candidates.items()]
