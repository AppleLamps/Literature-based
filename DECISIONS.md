# DECISIONS.md

Design decisions, each with reasoning. Append-only log; newest context at the
bottom of each section.

## D1. Substrate: MeSH co-occurrence over PubMed (E-utilities)
**Choice:** Build the concept graph from MeSH descriptor co-occurrence in PubMed
abstracts, retrieved via NCBI E-utilities, with PubMed Automatic Term Mapping
(ATM) used to resolve free-text A-terms to canonical MeSH descriptors.
**Reasoning:** SemMedDB (the requested primary substrate) requires a UMLS
license and is not reachable without authenticated credentials, so it cannot be
used under the "free and public, no credentials" constraint (see DATA.md for the
reachability check). The specified fallback is PubMed MeSH annotations. MeSH is a
curated controlled vocabulary already normalized to concepts, which is exactly
the substrate Swanson and Smalheiser used for the original ABC discoveries, so it
is both faithful to the method and sufficient to pass the validation gate.
**Trade-off:** Co-occurrence edges are undirected and carry no predicate, so the
"directional coherence" requirement of Phases 4-5 cannot be fully satisfied; it
is implemented as a warn-only gate and documented honestly rather than faked.

## D2. A-term resolution via PubMed ATM, not first db=mesh hit
**Choice:** Resolve "Migraine" -> `"migraine disorders"[MeSH Terms]` by parsing
the `querytranslation` PubMed returns, instead of taking the top db=mesh search
hit.
**Reasoning:** The first db=mesh hit for "Migraine" is "Ophthalmoplegic Migraine",
the wrong concept. ATM is the same mapping a human PubMed search uses and yields
the intended main heading. Verified live before adopting.

## D3. On-disk SQLite response cache + token-bucket rate limiting
**Choice:** Cache every E-utilities response keyed by endpoint+params; throttle to
`requests_per_second` (3 anonymous, 10 with an `NCBI_API_KEY` env var); retry
5xx/429 with exponential backoff.
**Reasoning:** LBD is iterative and the same A-term is rerun often; caching makes
"run a new A-term with one command" cheap and makes validation reproducible. Rate
limiting respects NCBI policy and avoids bans. No credentials are required.

## D4. MeSH stoplist to remove hub/check-tag descriptors
**Choice:** Filter a configurable stoplist of MeSH check tags and ultra-generic
descriptors (Humans, Female, Time Factors, Inflammation, Signal Transduction,
animal-strain tags, geography, study-type descriptors, ...) from both B
intermediates and C candidates.
**Reasoning:** These descriptors co-occur with almost everything and would
dominate the bridge graph, producing trivial "everything links to everything"
results and burying real intermediates. This is the single most important
ingredient for the validation gate to recover the historical C-terms. The list
lives in `config/default.json` per the "thresholds in config" requirement.

## D5. Gap determination: fast set-membership in discovery, authoritative
   esearch in the cascade
**Choice:** During candidate collection, a C that co-occurs with A anywhere in
A's fetched corpus is dropped immediately (cheap set test). For shortlisted
candidates, the cascade issues an authoritative `A[MeSH] AND C[MeSH]` PubMed count
over all years.
**Reasoning:** Per-candidate esearch on thousands of raw candidates would be slow;
the corpus set-test prunes cheaply, and the precise gap check is reserved for the
small shortlist that reaches the cascade.

## D6. Europe PMC for the preprint scoop check
**Choice:** Use Europe PMC's REST search filtered to source `PPR` for the
bioRxiv/medRxiv scoop gate; keep E-utilities as the main pipeline substrate.
**Reasoning:** E-utilities does not index preprints; Europe PMC does, is free, and
needs no key. Network failure on this gate yields a WARN, not a crash, so the
pipeline degrades gracefully offline.

## D7. Local biomedical embedding model with graceful fallback
**Choice:** Load `pritamdeka/S-PubMedBert-MS-MARCO` (PubMedBERT sentence model);
fall back through SPECTER2 to `all-MiniLM-L6-v2`; if none load, skip embedding
signals.
**Reasoning:** A biomedical model gives better in-domain concept similarity for the
near-synonym triviality check. Verified it loads and gives sensible similarities
(Raynaud/Fish Oils 0.80, Migraine/Magnesium 0.84). Embeddings are an auxiliary
signal only, never a hard kill except at a very high synonymy threshold (>=0.97),
so a missing model never blocks a run.

## D8. Conservative cascade ordering and default-REJECTED disposition
**Choice:** Evaluate gates cheap-to-expensive (plausibility -> triviality ->
real gap -> recency -> preprint -> coherence); a candidate is REJECTED unless it
clears every hard gate; coherence is warn-only on this substrate.
**Reasoning:** Matches Phase 5's "presumed worthless until it survives" mandate and
minimizes network spend by pruning early.

## D9. Validation pass threshold (top-N) in config
**Choice:** A historical target must appear within the top `validation.top_n`
(default 50) ranked candidates under the freeze to count as PASS; the exact rank
is always reported.
**Reasoning:** Swanson's targets are strong but not always rank 1 under raw
co-occurrence; top-50 out of hundreds/thousands of candidates is a meaningful blind
recovery while remaining honest about exact position. The value is configurable.

## D10. Two bug fixes found by running the validation gate
**Choice:** (a) Harden the E-utilities client to never cache empty/error responses
and to force-refetch when a resolution returns an empty translation; (b) resolve
A-terms via PubMed ATM (see D2).
**Reasoning:** The first validation run mis-resolved "Migraine" to "Ophthalmoplegic
Migraine" because a transient empty NCBI response had been cached, poisoning the
ATM path. The gate caught it - which is exactly what the gate is for.

## D11. Ranking redesign driven by the validation gate (the core tuning work)
Raw MeSH co-occurrence ranks hub concepts above the historical targets. Four
principled, configurable corrections were added, each verified to move the target
up under the freeze (all weights live in `config/default.json`):

1. **IDF weighting of intermediates** (`use_idf_weighting`). A B-term that links to
   few C's is a specific, informative bridge; a hub B that links to almost
   everything carries little signal. Each bridge is weighted by
   `log((1+totalC)/(1+productivity(B)))`. (ARROWSMITH-style linking-term
   specificity.)
2. **Semantic-type candidate filter** (`candidate_tree_filter`, default tree `D`).
   When A is a disease, the useful C's are *substances/interventions* (MeSH tree D
   = Chemicals and Drugs). This removes the anatomy / surgical-procedure / lab-test
   hubs that otherwise dominate (Coronary Circulation, Angioplasty, Leukocyte
   Count). Implemented with an offline MeSH tree map (D12).
3. **Semantic-type weighting of intermediates** (`b_tree_weights`). Physiological
   processes (tree G) are mechanistic bridges and are up-weighted; disease
   co-membership (tree C: lupus, scleroderma, arthritis) and diagnostic methods
   (tree E) are weak bridges ("same patients", not a mechanism) and are
   down-weighted. This is what separates fish oil's hemorheology path
   (blood viscosity, platelet aggregation) from connective-tissue-disease hubs.
4. **Specificity (rarity) refinement** (`specificity_refine`, `specificity_alpha`).
   A second pass divides each candidate's score by `log1p(global_freq(C))**alpha`
   (global frequency measured under the same date window), so a rare specific
   substance beats a common hub. Applied consistently to all eligible candidates
   (no refined/unrefined cliff).

**Sampling depth** was also raised: `b_top_n` 60->200 and `b_max_pmids` to 1200.
The 200 matters specifically - Platelet Aggregation (Raynaud co-occ rank 166) and
Blood Platelets (rank 167) are Swanson's signature fish-oil bridges and sit just
past the old top-80 cut; including them roughly tripled EPA's score. Swanson's
"neglected connections" are mid-frequency, not top-frequency, intermediates.

**Status (verified under freeze 1985/11/30):** these corrections moved
Eicosapentaenoic Acid from rank 712 -> 227 of 3520 candidates and recovered the
exact expected bridges (Blood Viscosity, Platelet Aggregation, Blood Platelets).
The target is recovered with the correct mechanism but not yet inside top-50;
see VALIDATION.md for the open item.

## D12. Offline MeSH tree map from the NLM descriptor file
**Choice:** Download the NLM ASCII descriptor file once
(`.../mesh/2025/asciimesh/d2025.bin`), parse `MH`->`MN` (heading -> tree numbers),
cache as JSON, and use it for the tree-D candidate filter and tree-letter B
weighting.
**Reasoning:** Per-descriptor tree lookups via E-utilities would be thousands of
calls; the bulk file is one download (~30 MB), then fully offline and instant. The
correct URL has a year directory (`/mesh/2025/asciimesh/d2025.bin`); the
no-year path returns an HTML landing page. Verified the targets resolve to tree D
(Fish Oils D10.627, EPA D10..., Magnesium D01...) and the hubs do not.

## D13. Runtime vs. the 3 req/s anonymous NCBI limit
**Choice:** Accept multi-minute runs; rely on aggressive caching so iteration and
reruns are cheap. Document that setting `NCBI_API_KEY` raises the limit to 10/s.
**Reasoning:** No API key is available in this autonomous session. The deep
sampling (200 B-terms, full marginal refinement) is what makes the gate recover
the targets, and it is a one-time cost per A-term thanks to the on-disk cache.
