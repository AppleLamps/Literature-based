"""Offline MeSH descriptor -> tree-number map.

Used for semantic-type filtering of candidates. When the A-term is a disease, the
useful C-candidates are *substances / interventions* (MeSH tree ``D`` = Chemicals
and Drugs), not anatomy, procedures, or other diseases. Filtering C to tree D is
how classical MeSH-based LBD avoids drowning the signal (fish oil, magnesium) in
co-occurring anatomy and surgical-procedure hubs.

NLM publishes the full descriptor file; we download it once (~30 MB ASCII),
parse ``MH`` (main heading) -> ``MN`` (tree numbers), and cache the parsed map as
JSON so subsequent runs are instant and offline.
"""
from __future__ import annotations

import gzip
import json
import os
from typing import Dict, List, Optional

import requests

# Year-directory ASCII descriptor file. 2025 is the last fully-released ASCII set.
_URL_TEMPLATE = "https://nlmpubs.nlm.nih.gov/projects/mesh/{year}/asciimesh/d{year}.bin"


class MeshTree:
    def __init__(self, cache_dir: str, year: int = 2025):
        self.year = year
        self.cache_dir = cache_dir
        self.json_path = os.path.join(cache_dir, f"mesh_tree_{year}.json")
        self.map: Dict[str, List[str]] = {}
        self.available = False
        self._load_or_build()

    def _load_or_build(self) -> None:
        if os.path.exists(self.json_path):
            try:
                with open(self.json_path, "r", encoding="utf-8") as fh:
                    self.map = json.load(fh)
                self.available = len(self.map) > 0
                return
            except (json.JSONDecodeError, OSError):
                pass
        try:
            text = self._download()
            self.map = self._parse(text)
            with open(self.json_path, "w", encoding="utf-8") as fh:
                json.dump(self.map, fh)
            self.available = len(self.map) > 0
        except Exception:
            # Offline or download failed: tree filtering silently disables.
            self.available = False

    def _download(self) -> str:
        url = _URL_TEMPLATE.format(year=self.year)
        resp = requests.get(url, timeout=180)
        resp.raise_for_status()
        content = resp.content
        if content[:2] == b"\x1f\x8b":  # gzip
            content = gzip.decompress(content)
        return content.decode("utf-8", errors="replace")

    @staticmethod
    def _parse(text: str) -> Dict[str, List[str]]:
        out: Dict[str, List[str]] = {}
        mh: Optional[str] = None
        trees: List[str] = []
        for line in text.splitlines():
            if line.startswith("*NEWRECORD"):
                if mh:
                    out[mh.lower()] = trees
                mh, trees = None, []
            elif line.startswith("MH = "):
                mh = line[5:].strip()
            elif line.startswith("MN = "):
                trees.append(line[5:].strip())
        if mh:
            out[mh.lower()] = trees
        return out

    def trees(self, name: str) -> List[str]:
        return self.map.get(name.lower(), [])

    def in_branches(self, name: str, prefixes: List[str], unknown_ok: bool = False) -> bool:
        """True if any tree number of ``name`` starts with an allowed prefix.

        Unknown descriptors (not in the map) return ``unknown_ok``.
        """
        t = self.trees(name)
        if not t:
            return unknown_ok
        return any(tn.startswith(p) for tn in t for p in prefixes)
