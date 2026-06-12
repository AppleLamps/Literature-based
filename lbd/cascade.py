"""Novelty and false-positive cascade (Phase 5).

Every candidate is presumed REJECTED until it clears, in order, a sequence of
gates. Each gate returns one of pass / fail / warn with a human-readable reason.
A candidate SURVIVES only if no gate returns ``fail``.

Gate order (cheap-to-expensive, so we prune before spending network calls):
  1. plausibility_floor   - enough independent bridges and minimum bridge score
  2. triviality           - name containment / near-synonym / MeSH-tree hierarchy
  3. real_gap             - A and C must not already co-occur in PubMed
  4. recency_scoop        - A+C not recently proposed (PubMed recent window)
  5. preprint_scoop       - A+C not proposed in bioRxiv/medRxiv (Europe PMC PPR)
  6. mechanistic_coherence- predicate-direction consistency (substrate-dependent)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .discovery import Candidate
from .eutils import EutilsClient
from .europepmc import EuropePmcClient


@dataclass
class GateResult:
    name: str
    status: str  # "pass" | "fail" | "warn"
    reason: str
    detail: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CascadeOutcome:
    disposition: str  # "SURVIVED" | "REJECTED"
    gates: List[GateResult] = field(default_factory=list)

    def as_flags(self) -> Dict[str, str]:
        return {g.name: g.status for g in self.gates}

    def fail_reason(self) -> str:
        for g in self.gates:
            if g.status == "fail":
                return f"{g.name}: {g.reason}"
        return ""


def _tree_numbers(client: EutilsClient, name: str) -> List[str]:
    resolved = client.resolve_mesh(name)
    if resolved:
        return resolved.get("tree_numbers", [])
    return []


def _is_hierarchical(a_trees: List[str], c_trees: List[str]) -> Optional[str]:
    """Return a reason string if A and C are parent/child in any MeSH tree."""
    for at in a_trees:
        for ct in c_trees:
            if at == ct:
                return f"identical MeSH tree node {at}"
            if ct.startswith(at + ".") or at.startswith(ct + "."):
                return f"MeSH hierarchy {at} <-> {ct}"
    return None


class Cascade:
    def __init__(
        self,
        cfg: Dict[str, Any],
        eutils: EutilsClient,
        europepmc: EuropePmcClient,
        embedder=None,
    ):
        self.cfg = cfg
        self.eutils = eutils
        self.europepmc = europepmc
        self.embedder = embedder
        self.field = cfg["corpus"]["a_literature_query_field"]
        # Tree numbers are resolved through the (cached) E-utilities client, but
        # the parse/regex work still repeats per concept. Memoize within a run so
        # the constant A-term is resolved once, not once per candidate.
        self._tree_cache: Dict[str, List[str]] = {}

    def _tree_numbers_cached(self, name: str) -> List[str]:
        cached = self._tree_cache.get(name)
        if cached is None:
            cached = _tree_numbers(self.eutils, name)
            self._tree_cache[name] = cached
        return cached

    # ---- individual gates ------------------------------------------------
    def gate_plausibility(self, cand: Candidate) -> GateResult:
        s = self.cfg["scoring"]
        if cand.n_bridges < s["min_bridge_terms"]:
            return GateResult(
                "plausibility_floor", "fail",
                f"only {cand.n_bridges} bridge(s); need >= {s['min_bridge_terms']}",
                {"n_bridges": cand.n_bridges},
            )
        if cand.score < s["min_bridge_score"]:
            return GateResult(
                "plausibility_floor", "fail",
                f"bridge score {cand.score:.3f} < floor {s['min_bridge_score']}",
                {"score": cand.score},
            )
        return GateResult(
            "plausibility_floor", "pass",
            f"{cand.n_bridges} independent bridges, score {cand.score:.2f}",
            {"n_bridges": cand.n_bridges, "score": cand.score},
        )

    def gate_triviality(self, a_name: str, cand: Candidate) -> GateResult:
        a_lc, c_lc = a_name.lower(), cand.c_term.lower()
        if a_lc in c_lc or c_lc in a_lc:
            return GateResult(
                "triviality", "fail",
                f"name containment between '{a_name}' and '{cand.c_term}'",
            )
        if self.embedder and self.embedder.available:
            sim = self.embedder.similarity(a_name, cand.c_term)
            cand.flags["embed_sim"] = round(sim, 3) if sim is not None else None
            if sim is not None and sim >= 0.97:
                return GateResult(
                    "triviality", "fail",
                    f"near-synonym (embedding cosine {sim:.3f} >= 0.97)",
                    {"embed_sim": sim},
                )
        if self.cfg["cascade"]["triviality_check_mesh_tree"]:
            a_trees = self._tree_numbers_cached(a_name)
            c_trees = self._tree_numbers_cached(cand.c_term)
            cand.flags["c_tree_numbers"] = c_trees
            reason = _is_hierarchical(a_trees, c_trees)
            if reason:
                return GateResult("triviality", "fail", f"trivial MeSH hierarchy: {reason}",
                                  {"a_trees": a_trees, "c_trees": c_trees})
        return GateResult("triviality", "pass", "distinct concept, not a hierarchy/synonym")

    def gate_real_gap(self, a_query: str, cand: Candidate) -> GateResult:
        f = self.field
        a_clause = f'"{a_query}"{f}' if f else a_query
        c_clause = f'"{cand.c_term}"{f}' if f else cand.c_term
        query = f"{a_clause} AND {c_clause}"
        n = self.eutils.esearch_count(query)
        cand.flags["ac_cooccurrence"] = n
        limit = self.cfg["cascade"]["gap_max_cooccurrence"]
        if n > limit:
            return GateResult(
                "real_gap", "fail",
                f"A and C already co-occur in {n} PubMed record(s) (limit {limit}) -> known, not hidden",
                {"cooccurrence": n, "query": query},
            )
        return GateResult("real_gap", "pass", f"A and C co-occur in {n} records (<= {limit}): genuine gap",
                          {"cooccurrence": n})

    def gate_recency_scoop(self, a_query: str, cand: Candidate, this_year: int) -> GateResult:
        f = self.field
        win = self.cfg["cascade"]["recency_window_years"]
        # The recent A+C count is a date-restricted subset of the all-dates count
        # already measured by gate_real_gap. If that total is zero, the recent
        # count is necessarily zero too, so skip the redundant esearch.
        if cand.flags.get("ac_cooccurrence") == 0:
            cand.flags["recent_pubmed_hits"] = 0
            return GateResult("recency_scoop", "pass",
                              f"no recent PubMed A+C link in last {win} years (A+C never co-occur)")
        a_clause = f'"{a_query}"{f}' if f else a_query
        c_clause = f'"{cand.c_term}"{f}' if f else cand.c_term
        mindate = f"{this_year - win}/01/01"
        n = self.eutils.esearch_count(f"{a_clause} AND {c_clause}", mindate=mindate, maxdate=f"{this_year}/12/31")
        cand.flags["recent_pubmed_hits"] = n
        if n > self.cfg["cascade"]["recency_max_hits"]:
            return GateResult("recency_scoop", "fail",
                              f"{n} PubMed record(s) link A+C in the last {win} years -> already proposed",
                              {"recent_hits": n})
        return GateResult("recency_scoop", "pass", f"no recent PubMed A+C link in last {win} years")

    def gate_preprint_scoop(self, a_query: str, cand: Candidate) -> GateResult:
        if not self.cfg["cascade"]["preprint_scoop_check"]:
            return GateResult("preprint_scoop", "warn", "preprint check disabled in config")
        # Europe PMC free-text within title/abstract; keep concepts quoted.
        query = f'"{a_query}" AND "{cand.c_term}"'
        res = self.europepmc.search_count(query, sources=["PPR"])
        if res["count"] < 0:
            return GateResult("preprint_scoop", "warn", "Europe PMC unreachable; preprint scoop unverified")
        cand.flags["preprint_hits"] = res["count"]
        if res["count"] > self.cfg["cascade"]["recency_max_hits"]:
            titles = "; ".join(s["title"][:80] for s in res["sample"][:2])
            return GateResult("preprint_scoop", "fail",
                              f"{res['count']} preprint(s) already link A+C: {titles}",
                              {"preprint_hits": res["count"], "sample": res["sample"]})
        return GateResult("preprint_scoop", "pass", "no bioRxiv/medRxiv preprint links A+C")

    def gate_coherence(self, cand: Candidate) -> GateResult:
        # Predicate direction is only available on a predicate-aware substrate
        # (SemMedDB). On the MeSH co-occurrence substrate we cannot verify the
        # sign of each edge, so we PASS-with-warning rather than silently claim
        # coherence. Honest by design (see DATA.md).
        directed = [b for b in cand.bridges if b.direction is not None]
        if not directed:
            return GateResult("mechanistic_coherence", "warn",
                              "predicate direction unavailable on MeSH substrate; chain not sign-verified")
        incoherent = [b for b in directed if b.direction == "contradictory"]
        if incoherent:
            return GateResult("mechanistic_coherence", "fail",
                              f"{len(incoherent)} contradictory predicate edge(s)")
        return GateResult("mechanistic_coherence", "pass", "predicate chain directionally consistent")

    # ---- driver ----------------------------------------------------------
    def evaluate(self, a_name: str, a_query: str, cand: Candidate, this_year: int) -> CascadeOutcome:
        gates: List[GateResult] = []

        g = self.gate_plausibility(cand)
        gates.append(g)
        if g.status == "fail":
            return CascadeOutcome("REJECTED", gates)

        g = self.gate_triviality(a_name, cand)
        gates.append(g)
        if g.status == "fail":
            return CascadeOutcome("REJECTED", gates)

        g = self.gate_real_gap(a_query, cand)
        gates.append(g)
        if g.status == "fail":
            return CascadeOutcome("REJECTED", gates)

        g = self.gate_recency_scoop(a_query, cand, this_year)
        gates.append(g)
        if g.status == "fail":
            return CascadeOutcome("REJECTED", gates)

        g = self.gate_preprint_scoop(a_query, cand)
        gates.append(g)
        if g.status == "fail":
            return CascadeOutcome("REJECTED", gates)

        gates.append(self.gate_coherence(cand))  # warn-only on MeSH substrate

        return CascadeOutcome("SURVIVED", gates)
