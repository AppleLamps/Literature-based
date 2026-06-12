"""Bridge scoring for ABC candidates.

A candidate C is scored by how strongly and how redundantly it is bridged to A:
  * strength of each A-B and B-C edge (co-occurrence counts, log-compressed so a
    handful of giant hub edges do not swamp many independent moderate bridges),
  * number of independent B-paths (diversity bonus),
  * optional directional coherence bonus (only when a predicate-aware substrate
    such as SemMedDB is live; neutral on the MeSH co-occurrence substrate).
"""
from __future__ import annotations

import math
from typing import Any, Dict, List

from .discovery import Candidate


def edge_weight(count: int, transform: str) -> float:
    if transform == "log1p":
        return math.log1p(count)
    if transform == "sqrt":
        return math.sqrt(count)
    return float(count)


def _idf(productivity: int, total_candidates: int) -> float:
    """IDF-style weight for an intermediate B by how many C's it links to."""
    return math.log((1.0 + total_candidates) / (1.0 + productivity))


def score_candidate(cand: Candidate, cfg: Dict[str, Any], total_candidates: int = 0) -> float:
    s = cfg["scoring"]
    t = s["edge_transform"]
    use_idf = s.get("use_idf_weighting", True)
    use_tree = s.get("use_b_tree_weighting", True)
    tree_w = s.get("b_tree_weights", {})
    tree_default = tree_w.get("default", 1.0)
    total = 0.0
    for b in cand.bridges:
        contrib = edge_weight(b.ab_count, t) * edge_weight(b.bc_count, t)
        if use_idf and total_candidates > 0:
            contrib *= _idf(b.b_productivity, total_candidates)
        if use_tree and b.b_tree:
            contrib *= tree_w.get(b.b_tree, tree_default)
        if b.direction in ("coherent", "consistent"):
            contrib *= 1.0 + s["directional_bonus"]
        total += contrib
    # Reward independent, diverse bridges.
    total *= 1.0 + 0.10 * math.log1p(cand.n_bridges)
    cand.score = total
    cand.flags["n_bridges"] = cand.n_bridges
    cand.flags["sum_ab"] = sum(b.ab_count for b in cand.bridges)
    cand.flags["sum_bc"] = sum(b.bc_count for b in cand.bridges)
    return total


def rank_candidates(candidates, cfg: Dict[str, Any]):
    total = len(candidates)
    for c in candidates:
        score_candidate(c, cfg, total_candidates=total)
        c.flags["raw_score"] = c.score
    return sorted(candidates, key=lambda c: c.score, reverse=True)


def refine_scores_with_specificity(
    candidates, cfg: Dict[str, Any], client, *, maxdate=None, mindate=None, field="[MeSH Terms]",
    progress=None,
):
    """Second pass: down-weight C-candidates by their global prevalence.

    Raw co-occurrence scoring rewards hub concepts (common anatomy/physiology)
    that co-occur with everything. We divide each candidate's raw bridge score by
    ``log1p(global_frequency_of_C) ** alpha`` so that a specific concept reached
    through several real bridges beats a generic hub reached through many trivial
    ones. The global frequency is measured under the same date window as the run,
    so the validation freeze and live runs stay consistent.

    Only candidates with at least ``specificity_min_bridges`` bridges are refined
    (singletons cannot survive the plausibility floor anyway); among those, the
    top ``specificity_refine_top_k`` by raw score get a marginal-frequency lookup.
    """
    s = cfg["scoring"]
    if not s.get("specificity_refine", True):
        return sorted(candidates, key=lambda c: c.score, reverse=True)
    alpha = s["specificity_alpha"]
    min_b = s["specificity_min_bridges"]
    eligible = [c for c in candidates if c.n_bridges >= min_b]
    eligible.sort(key=lambda c: c.flags.get("raw_score", c.score), reverse=True)
    shortlist = eligible[: s["specificity_refine_top_k"]]

    fetched: List[int] = []
    for i, c in enumerate(shortlist):
        clause = f'"{c.c_term}"{field}' if field else c.c_term
        n_c = client.esearch_count(clause, mindate=mindate, maxdate=maxdate)
        c.flags["c_marginal"] = n_c
        fetched.append(n_c)
        if progress:
            progress(i + 1, len(shortlist), c.c_term)

    # Estimate a background marginal for candidates we did not look up, so every
    # candidate is scored on one consistent scale (no refined/unrefined cliff).
    median_nc = 1000
    if fetched:
        sf = sorted(fetched)
        median_nc = sf[len(sf) // 2] or 1000

    for c in candidates:
        n_c = c.flags.get("c_marginal", median_nc)
        denom = math.log1p(max(n_c, 1)) ** alpha
        raw = c.flags.get("raw_score", c.score)
        c.score = raw / denom if denom > 0 else 0.0

    return sorted(candidates, key=lambda c: c.score, reverse=True)
