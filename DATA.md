# DATA.md - Data substrate

## What is live in this build

**Live substrate: MeSH-term co-occurrence over PubMed, via NCBI E-utilities.**

- Concepts = MeSH descriptors (a curated, pre-normalized controlled vocabulary).
- Edges = co-occurrence of two descriptors on the same PubMed record.
- Free-text A-terms are resolved to canonical MeSH descriptors through PubMed
  Automatic Term Mapping (e.g. `Migraine` -> `migraine disorders[MeSH Terms]`).
- Auxiliary endpoints:
  - **Europe PMC REST** (free, no key): used only for the preprint scoop check
    (`SRC:PPR` = bioRxiv/medRxiv) in the novelty cascade.
  - **Local PubMedBERT sentence model** (`pritamdeka/S-PubMedBert-MS-MARCO`):
    concept similarity for the near-synonym triviality check.
  - **NLM MeSH descriptor file** (`https://nlmpubs.nlm.nih.gov/projects/mesh/2025/asciimesh/d2025.bin`,
    free, no key): downloaded once (~30 MB ASCII), parsed to a descriptor ->
    tree-number map, cached as JSON. Drives the semantic-type candidate filter
    (keep C in tree D = Chemicals & Drugs) and the intermediate-term tree
    weighting. Note the URL needs the year directory; the no-year path serves an
    HTML page, not the file.

No credentials are used anywhere. No paid services. All compute is local except
the public PubMed / Europe PMC HTTP endpoints.

## Why not SemMedDB (the requested primary substrate)

SemMedDB (SemRep semantic predications over all of PubMed) is the ideal substrate
because it provides **typed, directed predicates** (e.g. `B CAUSES A`,
`C TREATS A`), which would let the mechanistic-coherence gate verify the *sign* of
each edge.

It was **not usable here** because it is gated behind a free **UMLS Metathesaurus
License** and authenticated download:

- SemMedDB is distributed as MySQL dump files from the NLM (Lister Hill /
  `lhncbc.nlm.nih.gov/ii/tools/SemRep_SemMedDB_SKR.html`), and downloads require
  signing in with UMLS Terminology Services (UTS) credentials.
- There is no anonymous public REST API for SemMedDB predications.
- The constraint for this project is "free and public data, no credentials," and
  the task explicitly says: if SemMedDB is not reachable without credentials, log
  that and fall back to PubMed MeSH / PubTator annotations.

**Reachability check performed:** NCBI E-utilities `einfo.fcgi` returned HTTP 200
(reachable). No unauthenticated SemMedDB predication endpoint exists, so the
documented fallback (PubMed MeSH co-occurrence) is the live substrate.

## Consequence for the method

| Capability | SemMedDB (unavailable) | MeSH co-occurrence (live) |
|---|---|---|
| Concept normalization | UMLS CUIs | MeSH descriptors (also normalized) |
| Edge existence | predication present | co-occurrence on a record |
| Edge **direction / predicate** | yes (SUBJECT-PREDICATE-OBJECT) | **no** (undirected) |
| Mechanistic-sign coherence gate | enforceable | **warn-only, not verified** |

The pipeline is written substrate-agnostically: a `Bridge` carries an optional
`direction`, and the scoring/cascade already consume it. If a UMLS-licensed user
later drops in a SemMedDB loader that populates `direction`, the directional bonus
and the coherence gate activate with no other code changes.

## PubTator3 (considered, not wired in)

PubTator3 offers free bioconcept annotations and a relations endpoint. It was
considered as an enrichment layer but not adopted for the core run because MeSH
co-occurrence is sufficient to pass the Swanson validation gate and keeps the
substrate single and auditable. It is the natural next enrichment if finer
entity typing (gene/chemical/disease) is wanted.

## Corpus sizes observed

Recorded per run in `REPORT.md` and `VALIDATION.md` (A-corpus size, number of
B-intermediates expanded, number of C-candidates collected).
