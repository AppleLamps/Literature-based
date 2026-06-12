"""Output generation (Phase 6): hypotheses.csv, per-survivor briefs, REPORT.md."""
from __future__ import annotations

import os
from typing import Any, Dict, List

import pandas as pd

from .cascade import CascadeOutcome
from .discovery import Candidate


def _pubmed_url(pmid: str) -> str:
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def candidate_row(a_name: str, cand: Candidate, outcome: CascadeOutcome, cfg: Dict[str, Any]) -> Dict[str, Any]:
    top_b = sorted(cand.bridges, key=lambda b: (b.ab_count * b.bc_count), reverse=True)
    top_b = top_b[: cfg["output"]["top_b_per_hypothesis"]]
    flags = outcome.as_flags()
    return {
        "A": a_name,
        "C": cand.c_term,
        "disposition": outcome.disposition,
        "bridge_score": round(cand.score, 4),
        "n_bridges": cand.n_bridges,
        "bridging_B_terms": " | ".join(b.b_term for b in top_b),
        "predicate_chains": " ; ".join(
            f"A-[{b.ab_count}]-{b.b_term}-[{b.bc_count}]-C" for b in top_b
        ),
        "AC_cooccurrence": cand.flags.get("ac_cooccurrence", ""),
        "recent_pubmed_hits": cand.flags.get("recent_pubmed_hits", ""),
        "preprint_hits": cand.flags.get("preprint_hits", ""),
        "embed_sim": cand.flags.get("embed_sim", ""),
        "gap_check": flags.get("real_gap", ""),
        "flag_plausibility": flags.get("plausibility_floor", ""),
        "flag_triviality": flags.get("triviality", ""),
        "flag_recency": flags.get("recency_scoop", ""),
        "flag_preprint": flags.get("preprint_scoop", ""),
        "flag_coherence": flags.get("mechanistic_coherence", ""),
        "reject_reason": outcome.fail_reason(),
    }


def write_hypotheses_csv(rows: List[Dict[str, Any]], path: str) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if df.empty:
        df.to_csv(path, index=False)
        return df
    # Rank by survival, then bridge score.
    df["_surv"] = (df["disposition"] == "SURVIVED").astype(int)
    df = df.sort_values(["_surv", "bridge_score"], ascending=[False, False]).drop(columns="_surv")
    df.to_csv(path, index=False)
    return df


def write_brief(a_name: str, cand: Candidate, outcome: CascadeOutcome, cfg: Dict[str, Any], path: str) -> None:
    top_b = sorted(cand.bridges, key=lambda b: (b.ab_count * b.bc_count), reverse=True)
    top_b = top_b[: cfg["output"]["top_b_per_hypothesis"]]
    n_edge = cfg["output"]["pmids_per_edge"]

    lines: List[str] = []
    lines.append(f"# Hypothesis brief: {a_name} -> {cand.c_term}\n")
    lines.append("## Proposed hypothesis (one sentence)\n")
    lines.append(
        f"**{cand.c_term}** may be mechanistically relevant to **{a_name}**, "
        f"linked indirectly through {cand.n_bridges} shared intermediate concept(s) "
        f"that connect the two otherwise-separate literatures.\n"
    )

    lines.append("## Mechanistic chain through the bridging B-terms\n")
    lines.append("Each row is an A->B->C path. PMIDs are representative records for each edge.\n")
    lines.append("| Intermediate B | A-B co-occ | B-C co-occ | A-B PMIDs | B-C PMIDs |")
    lines.append("|---|---|---|---|---|")
    for b in top_b:
        ab = ", ".join(f"[{p}]({_pubmed_url(p)})" for p in b.ab_pmids[:n_edge]) or "-"
        bc = ", ".join(f"[{p}]({_pubmed_url(p)})" for p in b.bc_pmids[:n_edge]) or "-"
        lines.append(f"| {b.b_term} | {b.ab_count} | {b.bc_count} | {ab} | {bc} |")
    lines.append("")

    lines.append("## Evidence that A-C is currently unstudied\n")
    ac = cand.flags.get("ac_cooccurrence", "n/a")
    rec = cand.flags.get("recent_pubmed_hits", "n/a")
    pre = cand.flags.get("preprint_hits", "n/a")
    lines.append(f"- Direct PubMed co-occurrence of A and C (all years): **{ac}** records.")
    lines.append(f"- A+C records in the recent PubMed window: **{rec}**.")
    lines.append(f"- A+C bioRxiv/medRxiv preprints (Europe PMC): **{pre}**.")
    sim = cand.flags.get("embed_sim")
    if sim is not None:
        lines.append(f"- Concept embedding cosine similarity (PubMedBERT): **{sim}** "
                     f"(low/moderate argues against a trivial synonymy).")
    lines.append("")

    lines.append("## Single most decisive test\n")
    bnames = ", ".join(b.b_term for b in top_b[:3])
    lines.append(
        f"Run a focused study measuring whether modulating **{cand.c_term}** changes the "
        f"intermediate mechanism(s) shared with **{a_name}** ({bnames}), then whether that "
        f"propagates to an {a_name} endpoint. Concretely: assemble the cohort/model where "
        f"both A and the strongest bridge ({top_b[0].b_term}) are measurable, apply or stratify "
        f"by {cand.c_term}, and test for the predicted change in {a_name} severity. If the "
        f"public-data shortcut is preferred, query a cohort/biobank (e.g. UK Biobank, NHANES, "
        f"or GEO expression sets) for the A--{cand.c_term} association conditioned on {top_b[0].b_term}.\n"
    )

    lines.append("## Cascade record\n")
    for g in outcome.gates:
        lines.append(f"- **{g.name}**: {g.status.upper()} - {g.reason}")
    lines.append("")
    lines.append("> This is a computationally generated hypothesis, not a finding. "
                 "It requires expert appraisal and experimental validation.\n")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def safe_filename(a_name: str, c_term: str) -> str:
    def clean(s: str) -> str:
        return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_")[:60]
    return f"{clean(a_name)}__{clean(c_term)}.md"
