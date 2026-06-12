"""Europe PMC client, used for the preprint scoop check.

NCBI E-utilities does not index bioRxiv/medRxiv. Europe PMC does (source code
``PPR``) alongside PubMed, and exposes a free REST search. We use it only in the
novelty cascade to ask "has anyone recently linked A and C, including in
preprints?". Documented in DATA.md.
"""
from __future__ import annotations

import json
import time
from typing import Any, Dict, Optional

import requests


class EuropePmcClient:
    def __init__(self, cfg: Dict[str, Any]):
        self.cfg = cfg["europepmc"]
        self.base = self.cfg["base_url"].rstrip("/")
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "lbd-pipeline (mailto:lucasmillen@gmail.com)"})

    def search_count(self, query: str, *, sources: Optional[list] = None) -> Dict[str, Any]:
        """Return {'count': int, 'sample': [titles]} for a Europe PMC query.

        ``sources`` filters by Europe PMC source code, e.g. ['PPR'] for preprints.
        """
        q = query
        if sources:
            src_clause = " OR ".join(f"SRC:{s}" for s in sources)
            q = f"({query}) AND ({src_clause})"
        params = {
            "query": q,
            "format": "json",
            "resultType": "lite",
            "pageSize": 10,
        }
        last_err = None
        for attempt in range(self.cfg["max_retries"]):
            try:
                resp = self.session.get(f"{self.base}/search", params=params, timeout=self.cfg["timeout"])
                if resp.status_code == 200:
                    data = json.loads(resp.text)
                    results = data.get("resultList", {}).get("result", [])
                    sample = [
                        {
                            "title": r.get("title", ""),
                            "source": r.get("source", ""),
                            "year": r.get("pubYear", ""),
                            "id": r.get("id", ""),
                        }
                        for r in results[:5]
                    ]
                    return {"count": int(data.get("hitCount", 0)), "sample": sample}
                raise requests.HTTPError(f"HTTP {resp.status_code}")
            except (requests.RequestException, requests.HTTPError) as exc:
                last_err = exc
                time.sleep(1.0 * (2 ** attempt))
        # Network failure should not crash the cascade; report unknown.
        return {"count": -1, "sample": [], "error": str(last_err)}
