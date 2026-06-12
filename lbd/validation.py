"""Validation gate (Phase 3): rediscover Swanson's findings under a date freeze.

No credit for new hypotheses until the known ones come back. For each historical
case we restrict the corpus to papers published before the discovery date and
require the target C-term to surface among the top-ranked candidates for the
given A-term.
"""
from __future__ import annotations

import os
from datetime import datetime, timezone
from typing import Any, Dict, List

from . import pipeline
from .config import repo_root

# Canonical Swanson benchmarks. Targets list MeSH aliases that count as a hit.
CASES: List[Dict[str, Any]] = [
    {
        "name": "Raynaud's syndrome -> dietary fish oil",
        "a_input": "Raynaud's syndrome",
        "maxdate": "1985/11/30",
        "freeze_label": "before 1986 (cutoff 1985/11/30)",
        "targets": [
            "Fish Oils", "Dietary Fats, Unsaturated", "Eicosapentaenoic Acid",
            "Docosahexaenoic Acids", "Fatty Acids, Omega-3", "Fatty Acids, Unsaturated",
        ],
        "expected_bridges": ["Blood Viscosity", "Platelet Aggregation", "Vascular Resistance", "Blood Platelets"],
        "reference": "Swanson DR. Fish oil, Raynaud's syndrome, and undiscovered public knowledge. Perspect Biol Med 1986.",
    },
    {
        "name": "Migraine -> magnesium",
        "a_input": "Migraine",
        "maxdate": "1987/12/31",
        "freeze_label": "before 1988 (cutoff 1987/12/31)",
        "targets": ["Magnesium"],
        "expected_bridges": [
            "Cortical Spreading Depression", "Platelet Aggregation", "Vasospasm",
            "Calcium Channel Blockers", "Calcium", "Serotonin", "Vascular Resistance",
        ],
        "reference": "Swanson DR. Migraine and magnesium: eleven neglected connections. Perspect Biol Med 1988.",
    },
]


def run_validation(cfg: Dict[str, Any], cache_dir: str, verbose: bool = True) -> Dict[str, Any]:
    eutils, europepmc, embedder = pipeline.build_clients(cfg, cache_dir, with_embedder=False)
    mesh_tree = pipeline.build_mesh_tree(cfg, cache_dir)
    top_n = cfg["validation"]["top_n"]
    results = []
    all_pass = True
    for case in CASES:
        if verbose:
            print(f"\n=== Validation: {case['name']} | freeze {case['freeze_label']} ===")
        res = pipeline.run_discovery(
            cfg, case["a_input"],
            eutils=eutils, europepmc=europepmc, embedder=None, mesh_tree=mesh_tree,
            maxdate=case["maxdate"], run_cascade=False, verbose=verbose,
        )
        hit = pipeline.find_target_rank(res, case["targets"])
        passed = hit is not None and hit["rank"] <= top_n
        all_pass = all_pass and passed
        recovered_bridges = []
        if hit:
            recovered_bridges = [b for b in hit["bridges"]
                                 if any(eb.lower() in b.lower() or b.lower() in eb.lower()
                                        for eb in case["expected_bridges"])]
        results.append({
            "case": case, "result": res, "hit": hit,
            "passed": passed, "recovered_bridges": recovered_bridges,
        })
        if verbose:
            if hit:
                print(f"  -> target '{hit['c_term']}' at RANK {hit['rank']}/{hit['total_candidates']} "
                      f"(score {hit['score']:.2f}); PASS={passed}")
                print(f"  -> recovered expected bridges: {recovered_bridges}")
            else:
                print(f"  -> target NOT FOUND among {res.n_candidates} candidates; PASS=False")
    return {"all_pass": all_pass, "results": results, "top_n": top_n}


def write_validation_md(val: Dict[str, Any], path: str) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    L: List[str] = []
    L.append("# VALIDATION.md - Swanson rediscovery under publication freeze\n")
    L.append(f"_Generated {now}._\n")
    L.append("The pipeline is not trusted for novel hypotheses until it blind-recovers known "
             "discoveries. For each case the PubMed corpus is frozen to papers published before "
             "the discovery date, and the historical C-term must appear in the top "
             f"{val['top_n']} ranked candidates for the A-term.\n")
    overall = "PASS" if val["all_pass"] else "FAIL"
    L.append(f"## Overall: **{overall}**\n")
    L.append("| Case | Freeze | Target found | Rank | Score | Recovered bridges | Result |")
    L.append("|---|---|---|---|---|---|---|")
    for r in val["results"]:
        c = r["case"]
        hit = r["hit"]
        if hit:
            rankcell = f"{hit['rank']}/{hit['total_candidates']}"
            tgt = hit["c_term"]
            score = f"{hit['score']:.2f}"
            bridges = ", ".join(r["recovered_bridges"][:4]) or "(none matched expected)"
        else:
            rankcell, tgt, score, bridges = "-", "NOT FOUND", "-", "-"
        res_txt = "PASS" if r["passed"] else "FAIL"
        L.append(f"| {c['name']} | {c['freeze_label']} | {tgt} | {rankcell} | {score} | {bridges} | {res_txt} |")
    L.append("")
    for r in val["results"]:
        c = r["case"]
        res = r["result"]
        hit = r["hit"]
        L.append(f"### {c['name']}\n")
        L.append(f"- A-term input: `{c['a_input']}` -> MeSH descriptor `{res.a_name}`")
        L.append(f"- Freeze: {c['freeze_label']}")
        L.append(f"- Reference: {c['reference']}")
        L.append(f"- Frozen A-corpus size: {res.corpus_size} PubMed records "
                 f"(fetched {res.fetched}); {res.n_b_terms} intermediate B-terms expanded; "
                 f"{res.n_candidates} C-candidates collected.")
        if hit:
            L.append(f"- **Target recovered: `{hit['c_term']}` at rank {hit['rank']} "
                     f"of {hit['total_candidates']}** (bridge score {hit['score']:.2f}, "
                     f"{hit['n_bridges']} independent bridges).")
            L.append(f"- Top bridging intermediates B (A->B->C): "
                     f"{', '.join(hit['bridges'][:10])}")
            L.append(f"- Expected bridges recovered: {', '.join(r['recovered_bridges']) or 'none of the listed set'}")
        else:
            L.append("- **Target not recovered.**")
        L.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
