"""Local sentence-embedding loader.

Used as an optional semantic signal in scoring and the triviality check: the
cosine similarity between the A and C concept strings. Very high similarity flags
near-synonyms (a triviality smell); moderate similarity is mild evidence the two
concepts are mechanistically relatable rather than random.

Model preference order (all run locally, CPU is fine):
  1. A biomedical model (PubMedBERT / SPECTER2 family) for in-domain text.
  2. Fallback to all-MiniLM-L6-v2 if the biomedical model cannot be downloaded.

Everything is lazy: if no model can be loaded (offline), the pipeline still runs
and simply skips the embedding-derived signals.
"""
from __future__ import annotations

from typing import List, Optional

# Tried in order. First that loads wins.
PREFERRED_MODELS = [
    "pritamdeka/S-PubMedBert-MS-MARCO",   # PubMedBERT sentence model, biomedical
    "allenai/specter2_base",              # SPECTER2 scientific-paper model
    "sentence-transformers/all-MiniLM-L6-v2",  # general fallback
]


class Embedder:
    def __init__(self, model_name: Optional[str] = None, quiet: bool = True):
        self.model = None
        self.model_name = None
        self._load(model_name, quiet)

    def _load(self, model_name: Optional[str], quiet: bool) -> None:
        try:
            from sentence_transformers import SentenceTransformer
        except Exception:
            return
        candidates = [model_name] if model_name else PREFERRED_MODELS
        for name in candidates:
            if not name:
                continue
            try:
                self.model = SentenceTransformer(name)
                self.model_name = name
                if not quiet:
                    print(f"[embeddings] loaded {name}")
                return
            except Exception as exc:  # download failure, offline, etc.
                if not quiet:
                    print(f"[embeddings] could not load {name}: {exc}")
                continue

    @property
    def available(self) -> bool:
        return self.model is not None

    def similarity(self, a: str, b: str) -> Optional[float]:
        if not self.available:
            return None
        import numpy as np

        emb = self.model.encode([a, b], normalize_embeddings=True)
        return float(np.dot(emb[0], emb[1]))

    def embed(self, texts: List[str]):
        if not self.available:
            return None
        return self.model.encode(texts, normalize_embeddings=True)


if __name__ == "__main__":
    emb = Embedder(quiet=False)
    if emb.available:
        print("model:", emb.model_name)
        print("sim(Raynaud Disease, Fish Oils) =", emb.similarity("Raynaud Disease", "Fish Oils"))
        print("sim(Migraine Disorders, Magnesium) =", emb.similarity("Migraine Disorders", "Magnesium"))
    else:
        print("no embedding model available; pipeline will skip embedding signals")
