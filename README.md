# LBD - Literature-Based Discovery (Swanson ABC) pipeline

Surfaces **hidden, testable connections** between biomedical concepts: pairs that
are each well studied and logically linkable through a shared intermediate, but
have never been studied together because they sit in separate literatures.

This is Swanson's **ABC model**:

```
A  --relates to-->  B  --relates to-->  C        with A--C unknown
```

If `A` relates to `B`, and `B` relates to `C`, but `A` and `C` have never been
connected in the literature, then **A--C is a candidate hypothesis**. The engine
discovers such candidates, then runs every one through a novelty / false-positive
cascade so that what you read has already survived the obvious ways to be wrong.

The engine is first **proven** by blind-rediscovering two of Swanson's own
historical findings under a strict publication-date freeze (see VALIDATION.md):

- **Raynaud's syndrome -> dietary fish oil** (freeze: before 1986)
- **Migraine -> magnesium** (freeze: before 1988)

## Substrate

Live substrate is **MeSH-term co-occurrence over PubMed** via NCBI E-utilities
(SemMedDB requires a UMLS license and is unavailable without credentials - see
DATA.md). Preprint scoop checks use Europe PMC; concept similarity uses a local
PubMedBERT model. Everything is free, public, and local; no credentials required.

## Install

Requires Python 3.11+ and ~3 GB disk for the local model/torch.

### Linux / macOS
```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

### Windows (PowerShell)
```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

(To install without the pinned file: `pip install requests pandas numpy scipy
networkx scikit-learn sentence-transformers`.)

### Optional: faster, higher-rate NCBI access
Set a free NCBI API key to raise the rate limit from 3 to 10 req/s:
```bash
export NCBI_API_KEY=xxxxxxxx      # Windows: $env:NCBI_API_KEY="xxxxxxxx"
```

## Prove it works (validation gate)

Run the Swanson rediscovery benchmark first. It writes `VALIDATION.md`.

```bash
.venv/bin/python -m lbd validate                 # Linux/macOS
.venv\Scripts\python.exe -m lbd validate         # Windows
```

This must pass (both targets recovered in the top-N under the freeze) before any
new hypothesis is trustworthy.

## Run discovery on a new A-term (one command)

```bash
.venv/bin/python -m lbd run --a-term "<concept>"
```
Windows:
```powershell
.venv\Scripts\python.exe -m lbd run --a-term "<concept>"
```

Examples:
```bash
.venv/bin/python -m lbd run --a-term "Alzheimer Disease"
.venv/bin/python -m lbd run --a-term "Psoriasis" --a-term "Endometriosis"
```

### Outputs
- `hypotheses.csv` - every scored candidate, ranked by survival then bridge score,
  with bridging B-terms, predicate-chain summary, gap-check result, and all
  cascade flags.
- `results/briefs/<A>__<C>.md` - a one-page brief per survivor: the one-sentence
  hypothesis, the mechanistic chain with representative PMIDs per edge, the
  evidence that A--C is unstudied, and the single most decisive experiment.
- `REPORT.md` - what was searched, corpus sizes, candidate/survivor counts,
  per-survivor confidence, and the explicit "these are hypotheses, not findings"
  statement.

## Useful flags

```
--config FILE          JSON overrides, deep-merged over config/default.json
--b-top-n N            number of B-intermediates to expand (default 60)
--a-max-pmids N        cap A-literature PMIDs fetched (default 1500)
--max-hypotheses N     how many top candidates to run through the cascade
--no-cascade           skip the cascade (raw candidate inspection)
--no-embeddings        skip loading the local model
--cache-dir DIR        on-disk response cache location (default ./cache)
--quiet                suppress progress
```

## How it works (pipeline stages)

1. **Gather A-literature** - resolve A to a MeSH descriptor (via PubMed automatic
   term mapping), fetch its PubMed records, tally co-occurring descriptors as
   B-intermediates (with link counts).
2. **Expand B -> C** - for each top B (default top 200), fetch its literature and
   tally its co-occurring descriptors as C-candidates. When A is a disease, C is
   filtered to substances/interventions (MeSH tree D = Chemicals & Drugs) via an
   offline MeSH tree map, so anatomy/procedure/disease hubs are excluded.
3. **Collect the gap** - keep only C that are *not* directly linked to A.
4. **Score** - each A-B-C bridge is weighted by (a) the edge strengths
   (log-compressed), (b) the IDF of the intermediate B (a B that links to few C's
   is informative; a hub B is not), and (c) the B-term's semantic type
   (physiological-process bridges up-weighted, disease-comembership/method bridges
   down-weighted). A second pass divides by the candidate's global rarity so a
   specific substance beats a common hub.
5. **Cascade** (every candidate REJECTED until it clears all):
   plausibility -> triviality -> real-gap -> recency -> preprint scoop -> coherence.
6. **Output** - ranked CSV, per-survivor briefs, REPORT.md.

All thresholds and weights live in `config/default.json` (corpus depth, the MeSH
stoplist, the tree-D candidate filter, the per-tree bridge weights, the
specificity exponent, and the cascade limits); nothing is hard-coded.

## Repository layout

```
lbd/
  cli.py          command-line entry point (python -m lbd ...)
  config.py       config loader (deep-merge overrides)
  eutils.py       cached, rate-limited NCBI E-utilities client
  europepmc.py    Europe PMC client (preprint scoop check)
  embeddings.py   local PubMedBERT / fallback sentence model
  mesh_tree.py    offline MeSH descriptor -> tree-number map (semantic filter)
  discovery.py    ABC engine: gather A, expand B, collect C
  scoring.py      bridge scoring (IDF + tree weight + rarity) + ranking
  cascade.py      novelty / false-positive cascade
  validation.py   Swanson rediscovery harness
  pipeline.py     orchestration shared by validate + run
  output.py       hypotheses.csv + per-survivor briefs
  report.py       REPORT.md
config/default.json   all thresholds, weights, and the MeSH stoplist
```

## Important caveats

These are **hypotheses generated from literature structure, not findings**. The
MeSH co-occurrence substrate is undirected, so mechanistic *direction* is not
verified (the coherence gate is warn-only). A "gap" means absence of evidence in
the indexed literature, not proof the link is novel or correct. Every survivor
needs expert appraisal and experimental validation. See REPORT.md and DATA.md.
