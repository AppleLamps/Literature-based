"""End-to-end orchestration for the LBD pipeline."""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from . import discovery, scoring
from .cascade import Cascade, CascadeOutcome
from .config import load_config, repo_root
from .discovery import ALiterature, Candidate
from .eutils import EutilsClient
from .europepmc import EuropePmcClient


@dataclass
class DiscoveryResult:
    a_input: str
    a_name: str
    a_query: str
    corpus_size: int
    fetched: int
    n_b_terms: int
    b_terms_used: List[str]
    n_candidates: int
    ranked: List[Candidate]
    outcomes: Dict[str, CascadeOutcome] = field(default_factory=dict)


def build_clients(cfg: Dict[str, Any], cache_dir: str, with_embedder: bool = True):
    os.makedirs(cache_dir, exist_ok=True)
    eutils = EutilsClient(cfg, os.path.join(cache_dir, "eutils_cache.sqlite"))
    europepmc = EuropePmcClient(cfg)
    embedder = None
    if with_embedder:
        from .embeddings import Embedder
        embedder = Embedder(quiet=True)
    return eutils, europepmc, embedder


def build_mesh_tree(cfg: Dict[str, Any], cache_dir: str):
    from .mesh_tree import MeshTree
    return MeshTree(cache_dir, year=cfg["corpus"].get("mesh_year", 2025))


def _progress(stage: str):
    def fn(i, n, b):
        sys.stderr.write(f"\r  [{stage}] {i}/{n}  expanding B: {b[:40]:<40}")
        sys.stderr.flush()
        if i == n:
            sys.stderr.write("\n")
    return fn


def run_discovery(
    cfg: Dict[str, Any],
    a_input: str,
    *,
    eutils: EutilsClient,
    europepmc: EuropePmcClient,
    embedder=None,
    mesh_tree=None,
    maxdate: Optional[str] = None,
    mindate: Optional[str] = None,
    run_cascade: bool = True,
    this_year: int = 2026,
    cascade_top_k: Optional[int] = None,
    verbose: bool = True,
) -> DiscoveryResult:
    """Core discovery for one A-term. Used by both validation and live runs."""
    resolved = eutils.resolve_mesh(a_input)
    a_name = resolved["name"] if resolved else a_input
    a_query = a_name  # query the canonical MeSH descriptor
    if verbose:
        print(f"  A-term '{a_input}' -> MeSH descriptor '{a_name}'", file=sys.stderr)

    a_lit = discovery.gather_a_literature(
        eutils, cfg, a_query, a_name, maxdate=maxdate, mindate=mindate, mesh_tree=mesh_tree
    )
    if verbose:
        print(f"  A corpus: {a_lit.corpus_size} records "
              f"(fetched {len(a_lit.fetched_pmids)}), {len(a_lit.b_counts)} B-terms above threshold",
              file=sys.stderr)

    b_terms = discovery.select_b_terms(a_lit, cfg)
    candidates = discovery.expand_b_and_collect_c(
        eutils, cfg, a_lit, b_terms,
        maxdate=maxdate, mindate=mindate, mesh_tree=mesh_tree,
        progress=_progress("B->C") if verbose else None,
    )
    scoring.rank_candidates(candidates, cfg)
    ranked = scoring.refine_scores_with_specificity(
        candidates, cfg, eutils,
        maxdate=maxdate, mindate=mindate,
        field=cfg["corpus"]["a_literature_query_field"],
        progress=_progress("specificity") if verbose else None,
    )
    if verbose:
        print(f"  Collected {len(ranked)} C-candidates "
              f"({sum(1 for c in candidates if c.n_bridges >= 2)} with >=2 bridges)", file=sys.stderr)

    result = DiscoveryResult(
        a_input=a_input, a_name=a_name, a_query=a_query,
        corpus_size=a_lit.corpus_size, fetched=len(a_lit.fetched_pmids),
        n_b_terms=len(b_terms), b_terms_used=b_terms,
        n_candidates=len(ranked), ranked=ranked,
    )

    if run_cascade:
        cascade = Cascade(cfg, eutils, europepmc, embedder, mesh_tree=mesh_tree)
        top_k = cascade_top_k or cfg["output"]["max_hypotheses"]
        if verbose:
            print(f"  Running novelty/false-positive cascade on top {top_k} candidates...",
                  file=sys.stderr)
        for i, cand in enumerate(ranked[:top_k]):
            outcome = cascade.evaluate(a_name, a_query, cand, this_year)
            result.outcomes[cand.c_term] = outcome
            if verbose:
                sys.stderr.write(
                    f"\r    cascade {i+1}/{min(top_k, len(ranked))}: "
                    f"{cand.c_term[:30]:<30} -> {outcome.disposition}        "
                )
                sys.stderr.flush()
        if verbose:
            sys.stderr.write("\n")
    return result


def find_target_rank(result: DiscoveryResult, target_aliases: List[str]) -> Optional[Dict[str, Any]]:
    """Locate a target C-term in the ranked list by exact MeSH-descriptor match."""
    aliases = {a.lower() for a in target_aliases}
    for rank, cand in enumerate(result.ranked, start=1):
        cl = cand.c_term.lower()
        if cl in aliases:
            return {
                "rank": rank,
                "c_term": cand.c_term,
                "score": cand.score,
                "n_bridges": cand.n_bridges,
                "bridges": [b.b_term for b in sorted(cand.bridges, key=lambda x: x.ab_count*x.bc_count, reverse=True)],
                "total_candidates": len(result.ranked),
            }
    return None
