"""NCBI E-utilities client.

A thin, polite, cached wrapper over the Entrez E-utilities REST endpoints.

Design choices (logged in DECISIONS.md):
  * On-disk SQLite response cache keyed by endpoint + sorted params, so a rerun
    of the same A-term costs zero network calls. LBD is inherently iterative;
    caching is what makes "run a new A-term with one command" practical.
  * Token-bucket rate limiting at ``requests_per_second`` (3/s anonymous,
    10/s with an NCBI API key) to stay inside NCBI's usage policy.
  * Exponential backoff on 429/5xx and transient network errors.
"""
from __future__ import annotations

import hashlib
import json
import re
import sqlite3
import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple

import requests


class RateLimiter:
    def __init__(self, rate_per_second: float):
        self.min_interval = 1.0 / max(rate_per_second, 0.1)
        self._last = 0.0

    def wait(self) -> None:
        now = time.monotonic()
        delta = now - self._last
        if delta < self.min_interval:
            time.sleep(self.min_interval - delta)
        self._last = time.monotonic()


class ResponseCache:
    def __init__(self, path: str):
        self.path = path
        self._conn = sqlite3.connect(path)
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS cache (key TEXT PRIMARY KEY, body TEXT, ts REAL)"
        )
        self._conn.commit()

    @staticmethod
    def make_key(endpoint: str, params: Dict[str, Any]) -> str:
        # Exclude volatile identity params from the cache key.
        relevant = {k: v for k, v in params.items() if k not in ("api_key", "email", "tool")}
        blob = endpoint + "|" + json.dumps(relevant, sort_keys=True)
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    def get(self, key: str) -> Optional[str]:
        row = self._conn.execute("SELECT body FROM cache WHERE key = ?", (key,)).fetchone()
        return row[0] if row else None

    def put(self, key: str, body: str) -> None:
        self._conn.execute(
            "INSERT OR REPLACE INTO cache (key, body, ts) VALUES (?, ?, ?)",
            (key, body, time.time()),
        )
        self._conn.commit()


class EutilsClient:
    def __init__(self, cfg: Dict[str, Any], cache_path: str):
        self.cfg = cfg["eutils"]
        self.base = self.cfg["base_url"].rstrip("/")
        self.limiter = RateLimiter(self.cfg["requests_per_second"])
        self.cache = ResponseCache(cache_path)
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": f"{self.cfg['tool']} (mailto:{self.cfg['email']})"})
        self.network_calls = 0

    # ---- low level -------------------------------------------------------
    def _identity(self) -> Dict[str, Any]:
        ident = {"tool": self.cfg["tool"], "email": self.cfg["email"]}
        if self.cfg.get("api_key"):
            ident["api_key"] = self.cfg["api_key"]
        return ident

    @staticmethod
    def _looks_malformed(endpoint: str, body: str) -> bool:
        """Detect empty/error responses we must not cache."""
        if not body or len(body.strip()) < 2:
            return True
        if endpoint.startswith("esearch"):
            try:
                res = json.loads(body).get("esearchresult", {})
            except json.JSONDecodeError:
                return True
            if "ERROR" in res or "error" in res:
                return True
        return False

    def _request(self, endpoint: str, params: Dict[str, Any], method: str = "GET", force: bool = False) -> str:
        key = ResponseCache.make_key(endpoint, params)
        if not force:
            cached = self.cache.get(key)
            if cached is not None:
                return cached
        full_params = dict(params)
        full_params.update(self._identity())
        url = f"{self.base}/{endpoint}"
        last_err: Optional[Exception] = None
        for attempt in range(self.cfg["max_retries"]):
            self.limiter.wait()
            try:
                if method == "POST":
                    resp = self.session.post(url, data=full_params, timeout=self.cfg["timeout"])
                else:
                    resp = self.session.get(url, params=full_params, timeout=self.cfg["timeout"])
                self.network_calls += 1
                if resp.status_code == 200:
                    if self._looks_malformed(endpoint, resp.text):
                        # Empty/error payload: retry rather than cache the poison.
                        raise requests.HTTPError("malformed/empty E-utilities response")
                    self.cache.put(key, resp.text)
                    return resp.text
                if resp.status_code in (429, 500, 502, 503, 504):
                    raise requests.HTTPError(f"HTTP {resp.status_code}")
                resp.raise_for_status()
            except (requests.RequestException, requests.HTTPError) as exc:
                last_err = exc
                sleep_s = self.cfg["backoff_base_seconds"] * (2 ** attempt)
                time.sleep(sleep_s)
        raise RuntimeError(f"E-utilities request failed after retries: {endpoint} {params}: {last_err}")

    # ---- high level ------------------------------------------------------
    def esearch(
        self,
        term: str,
        *,
        db: str = "pubmed",
        retmax: int = 100000,
        mindate: Optional[str] = None,
        maxdate: Optional[str] = None,
        datetype: str = "pdat",
    ) -> Tuple[int, List[str]]:
        """Return (total_count, pmids). Dates as YYYY/MM/DD."""
        params: Dict[str, Any] = {
            "db": db,
            "term": term,
            "retmax": retmax,
            "retmode": "json",
            "retstart": 0,
        }
        if mindate or maxdate:
            params["datetype"] = datetype
            params["mindate"] = mindate or "1700/01/01"
            params["maxdate"] = maxdate or "2100/12/31"
        body = self._request("esearch.fcgi", params)
        data = json.loads(body)
        result = data.get("esearchresult", {})
        count = int(result.get("count", 0))
        ids = result.get("idlist", [])
        return count, ids

    def esearch_count(self, term: str, **kwargs) -> int:
        kwargs.setdefault("retmax", 0)
        count, _ = self.esearch(term, **kwargs)
        return count

    def efetch_mesh(self, pmids: List[str]) -> Dict[str, Dict[str, Any]]:
        """Fetch MeSH headings for PMIDs.

        Returns ``{pmid: {"mesh": [{"ui","name","major"}...], "year": int|None}}``.
        Batched via POST to respect URL length limits.
        """
        out: Dict[str, Dict[str, Any]] = {}
        batch = self.cfg["efetch_batch_size"]
        for i in range(0, len(pmids), batch):
            chunk = pmids[i : i + batch]
            params = {
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "xml",
                "rettype": "medline",
            }
            body = self._request("efetch.fcgi", params, method="POST")
            out.update(self._parse_mesh_xml(body))
        return out

    @staticmethod
    def _parse_mesh_xml(xml_text: str) -> Dict[str, Dict[str, Any]]:
        out: Dict[str, Dict[str, Any]] = {}
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return out
        for article in root.findall(".//PubmedArticle"):
            pmid_el = article.find(".//PMID")
            if pmid_el is None or not pmid_el.text:
                continue
            pmid = pmid_el.text.strip()
            mesh: List[Dict[str, Any]] = []
            for heading in article.findall(".//MeshHeadingList/MeshHeading"):
                desc = heading.find("DescriptorName")
                if desc is None or not desc.text:
                    continue
                major = desc.get("MajorTopicYN", "N") == "Y"
                # A heading is "major" if the descriptor or any qualifier is major.
                if not major:
                    for qual in heading.findall("QualifierName"):
                        if qual.get("MajorTopicYN", "N") == "Y":
                            major = True
                            break
                mesh.append({"ui": desc.get("UI", ""), "name": desc.text.strip(), "major": major})
            year = None
            ystr = article.find(".//DateCompleted/Year")
            if ystr is None:
                ystr = article.find(".//PubDate/Year")
            if ystr is not None and ystr.text and ystr.text.isdigit():
                year = int(ystr.text)
            out[pmid] = {"mesh": mesh, "year": year}
        return out

    # ---- MeSH term resolution -------------------------------------------
    _MESH_TERM_RE = re.compile(r'"([^"]+)"\[MeSH Terms\]', re.IGNORECASE)

    @staticmethod
    def _extract_qt(body: str) -> str:
        try:
            return json.loads(body).get("esearchresult", {}).get("querytranslation", "")
        except json.JSONDecodeError:
            return ""

    def resolve_mesh(self, term: str) -> Optional[Dict[str, Any]]:
        """Resolve a free-text concept to a canonical MeSH descriptor.

        Uses PubMed's Automatic Term Mapping (the same translation a normal
        PubMed search applies), which maps e.g. "Migraine" -> "migraine
        disorders"[MeSH Terms]. This is far more reliable than picking the first
        db=mesh hit. Returns ``{"name","ui","uid","tree_numbers"}`` or ``None``.
        """
        params = {"db": "pubmed", "term": term, "retmax": 0, "retmode": "json"}
        body = self._request("esearch.fcgi", params)
        qt = self._extract_qt(body)
        if not qt:  # poisoned/empty cache or transient: force a fresh fetch
            body = self._request("esearch.fcgi", params, force=True)
            qt = self._extract_qt(body)
        match = self._MESH_TERM_RE.search(qt)
        descriptor = match.group(1) if match else term
        # Canonicalise capitalisation and fetch UI + tree numbers from db=mesh.
        _, ids = self.esearch(descriptor, db="mesh", retmax=1)
        uid = ids[0] if ids else None
        name, ui, tree_numbers = descriptor, "", []
        if uid:
            name, ui, tree_numbers = self._fetch_mesh_record(uid, fallback_name=descriptor)
        return {"name": name, "ui": ui, "uid": uid, "tree_numbers": tree_numbers}

    def _fetch_mesh_record(self, uid: str, fallback_name: str) -> Tuple[str, str, List[str]]:
        body = self._request("efetch.fcgi", {"db": "mesh", "id": uid, "retmode": "text"})
        name = fallback_name
        ui = ""
        trees: List[str] = []
        lines = body.splitlines()
        for line in lines:
            s = line.strip()
            if not s:
                continue
            # First non-empty line is like "1: Raynaud Disease".
            head = s.split(":", 1)
            if len(head) == 2 and head[0].strip().isdigit() and head[1].strip():
                name = head[1].strip()
            break
        for line in lines:
            s = line.strip()
            if s.startswith("Tree Number(s):"):
                payload = s.split(":", 1)[1]
                trees = [t.strip() for t in payload.split(",") if t.strip()]
            elif s.startswith("Unique ID:"):
                ui = s.split(":", 1)[1].strip()
        return name, ui, trees
