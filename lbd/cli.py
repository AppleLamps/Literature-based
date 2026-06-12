"""Command-line interface for the LBD pipeline.

Canonical commands:
    .venv/bin/python -m lbd validate
    .venv/bin/python -m lbd run --a-term "<concept>"

Windows:
    .venv\\Scripts\\python.exe -m lbd validate
    .venv\\Scripts\\python.exe -m lbd run --a-term "<concept>"
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from typing import Any, Dict, List

from . import output, pipeline, report, validation
from .config import load_config, repo_root


def _cache_dir(args) -> str:
    return args.cache_dir or os.path.join(repo_root(), "cache")


def cmd_validate(args) -> int:
    cfg = load_config(args.config)
    cache_dir = _cache_dir(args)
    val = validation.run_validation(cfg, cache_dir, verbose=not args.quiet)
    out_path = os.path.join(repo_root(), "VALIDATION.md")
    validation.write_validation_md(val, out_path)
    print(f"\nValidation {'PASSED' if val['all_pass'] else 'FAILED'}. Wrote {out_path}")
    for r in val["results"]:
        hit = r["hit"]
        if hit:
            print(f"  - {r['case']['name']}: rank {hit['rank']}/{hit['total_candidates']} "
                  f"({hit['c_term']}) -> {'PASS' if r['passed'] else 'FAIL'}")
        else:
            print(f"  - {r['case']['name']}: NOT FOUND -> FAIL")
    return 0 if val["all_pass"] else 1


def _validation_summary_line(cfg, cache_dir) -> str:
    """Read VALIDATION.md if present; else note it was not run."""
    vpath = os.path.join(repo_root(), "VALIDATION.md")
    if os.path.exists(vpath):
        with open(vpath, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("## Overall:"):
                    return f"See VALIDATION.md. {line.strip()}"
        return "See VALIDATION.md."
    return "Validation not run in this session (run `python -m lbd validate`)."


def cmd_run(args) -> int:
    cfg = load_config(args.config)
    # CLI overrides for common knobs.
    if args.max_hypotheses is not None:
        cfg["output"]["max_hypotheses"] = args.max_hypotheses
    if args.b_top_n is not None:
        cfg["corpus"]["b_top_n"] = args.b_top_n
    if args.a_max_pmids is not None:
        cfg["corpus"]["a_max_pmids"] = args.a_max_pmids

    cache_dir = _cache_dir(args)
    eutils, europepmc, embedder = pipeline.build_clients(cfg, cache_dir, with_embedder=not args.no_embeddings)
    mesh_tree = pipeline.build_mesh_tree(cfg, cache_dir)
    this_year = datetime.now(timezone.utc).year if args.this_year is None else args.this_year

    runs: List[pipeline.DiscoveryResult] = []
    rows: List[Dict[str, Any]] = []
    survivors_by_a: Dict[str, list] = {}

    results_dir = os.path.join(repo_root(), "results")
    briefs_dir = os.path.join(results_dir, "briefs")
    os.makedirs(briefs_dir, exist_ok=True)

    for a_term in args.a_term:
        print(f"\n=== Discovery run: A = '{a_term}' ===", file=sys.stderr)
        res = pipeline.run_discovery(
            cfg, a_term,
            eutils=eutils, europepmc=europepmc, embedder=embedder, mesh_tree=mesh_tree,
            run_cascade=not args.no_cascade, this_year=this_year,
            cascade_top_k=cfg["output"]["max_hypotheses"],
            verbose=not args.quiet,
        )
        runs.append(res)
        survivors: list = []
        for cand in res.ranked[: cfg["output"]["max_hypotheses"]]:
            outcome = res.outcomes.get(cand.c_term)
            if outcome is None:
                continue
            rows.append(output.candidate_row(res.a_name, cand, outcome, cfg))
            if outcome.disposition == "SURVIVED":
                survivors.append((cand, outcome))
        survivors_by_a[res.a_name] = survivors
        # Write a brief per survivor.
        for cand, outcome in survivors:
            fname = output.safe_filename(res.a_name, cand.c_term)
            output.write_brief(res.a_name, cand, outcome, cfg, os.path.join(briefs_dir, fname))
        print(f"  -> {len(survivors)} survivor(s) for {res.a_name}", file=sys.stderr)

    csv_path = os.path.join(repo_root(), "hypotheses.csv")
    df = output.write_hypotheses_csv(rows, csv_path)
    report_path = os.path.join(repo_root(), "REPORT.md")
    report.write_report(runs, survivors_by_a, cfg, report_path,
                        _validation_summary_line(cfg, cache_dir))

    total_surv = sum(len(v) for v in survivors_by_a.values())
    print(f"\nWrote {csv_path} ({len(df)} rows), {report_path}, and "
          f"{total_surv} brief(s) in results/briefs/.")
    print(f"Network calls this run: {eutils.network_calls} (rest served from cache).")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="lbd", description="Literature-based discovery (Swanson ABC) pipeline.")
    p.add_argument("--config", help="Path to a JSON config override file.")
    p.add_argument("--cache-dir", help="Directory for the on-disk response cache.")
    p.add_argument("--quiet", action="store_true", help="Suppress progress output.")
    sub = p.add_subparsers(dest="command", required=True)

    pv = sub.add_parser("validate", help="Run the Swanson rediscovery validation gate.")
    pv.set_defaults(func=cmd_validate)

    pr = sub.add_parser("run", help="Run discovery for one or more A-terms.")
    pr.add_argument("--a-term", action="append", required=True,
                    help="Disease/phenotype A-term. Repeat for multiple.")
    pr.add_argument("--max-hypotheses", type=int, dest="max_hypotheses",
                    help="Cascade and report at most this many top candidates.")
    pr.add_argument("--b-top-n", type=int, dest="b_top_n", help="Number of B-intermediates to expand.")
    pr.add_argument("--a-max-pmids", type=int, dest="a_max_pmids", help="Max A-literature PMIDs to fetch.")
    pr.add_argument("--no-cascade", action="store_true", help="Skip the novelty cascade (debug).")
    pr.add_argument("--no-embeddings", action="store_true", help="Skip loading the embedding model.")
    pr.add_argument("--this-year", type=int, dest="this_year", help="Override 'current year' for recency.")
    pr.set_defaults(func=cmd_run)
    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
