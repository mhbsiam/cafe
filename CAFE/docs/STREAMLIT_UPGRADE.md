# Streamlit 1.58 Upgrade Notes & Feature Roadmap

This document accompanies the dependency bump from **Streamlit 1.41.0 → 1.58.0** and
records (a) what changed for CAFE, and (b) the newer Streamlit capabilities worth adopting
in follow-up work. No feature below is implemented yet — this is a prioritized backlog.

## Why the scientific stack did not move

The upgrade was deliberately conservative. The scientific stack is **frozen at its
numpy-1.x ceiling**:

- `scanpy` **1.11+ and 1.12 all require `numpy>=2` and `pandas>=2.3`**. `scanpy==1.10.4`
  (our pin) is already the newest release compatible with `numpy==1.26.4`.
- `plotly` 6 pulls in `kaleido` v1, which **removes the `engine=` argument** used by
  ~10 `fig.write_image(..., engine='kaleido')` call sites in `bin/Tools/Visualization.py`.

So `numpy`, `scanpy`, `anndata`, `plotly`, and `kaleido` stay on their current major lines.
Moving them is a separate **numpy-2 migration** (see roadmap at the bottom).

## What Streamlit 1.58 changes for CAFE

The app is already modern (`st.navigation` / `st.Page`, 36 `@st.fragment` render
functions, no deprecated APIs), so **no code changes were required** to run on 1.58.

**Server change (main validation item):** 1.58 replaces Tornado with **Starlette/Uvicorn**
by default (tornado is no longer a dependency). CAFE only configures
`--server.maxUploadSize` (3000 MB desktop / 2000 MB web) and a theme, so nothing needed
migration — but **large-CSV upload over the new server is the #1 thing to smoke-test**.

**Free wins from 1.58:**
- Fragment apps no longer crash on stale auto-reruns (#15130) — matters for our 36
  fragments.
- Tables / dataframes / data-editors no longer overscroll (#15309).
- `st.selectbox` no longer hides its first option at exactly seven selections (#14997).
- Faster startup (reduced external-IP-lookup timeout).

**Behavioral changes to re-check:**
- Widget/icon sizing is more consistent (#15056, #15098) — we ship custom CSS in
  `bin/theme.py`; do a visual pass for spacing/layout regressions.
- `st.multiselect` disables "Select all" for very large option lists (#15301) — our marker
  multiselects can carry 20–100+ options.
- `st.markdown` now blocks `javascript:`/`vbscript:` links — low risk; glance at any custom
  HTML markdown.

## Feature adoption backlog (prioritized)

### 1. `@st.fragment(parallel=True)` — recommended first (new in 1.58)
Lets per-tab fragment computations run **concurrently** for snappier tab switching. CAFE
already groups its 36 `@st.fragment` render functions under 8 `st.tabs()` blocks:

| Location (`bin/Tools/Visualization.py`) | Tab group |
|---|---|
| line 320 | UMAP · Clusters / Marker Expression / Selected Clusters / Single Marker / Per Sample |
| line 2254 | Cluster Heatmaps / Composition Bar Charts / Cell % & Count / Compare Two Clusters |
| line 2704 | Sankey · Clusters↔Groups / Marker Threshold / Cells vs Clusters / Cell Count by Group |
| line 2721 | Dot Plot+Dendrogram / Expression Heatmap / Histogram / Group×Cluster / Marker Profile ×2 |
| line 3544 | Group Violin / Marker Bar ×2 / KDE / Relative Abundance / Correlation Matrix |
| line 3798 | Compare 3 Clusters / Find Clusters by Marker / One vs All / Cluster Frequency / Violin+Bar |
| line 3819 | Group Comparison / Marker-hi / Effect Size / Diff. Abundance / Normality / Group Corr. |

Also the per-group loop inside `render_tab8` (`Visualization.py:2502`) is a candidate.

**Caveat:** each parallel fragment must read `st.session_state.adata` locally and avoid
shared mutable state / ordering assumptions between fragments.

### 2. `st.pagination` (new in 1.58)
Page large tables currently rendered all at once: expression tables
(`Visualization.py:141, 156, 169`), normality-test results (~`4461`), rank-genes output
(~`2758`). Hundreds–thousands of rows today.

### 3. `@st.cache_data` / `@st.cache_resource`
Replace ad-hoc `st.session_state` memoization on expensive scanpy ops
(`bin/Tools/Data_Processing.py` PCA/UMAP/Leiden ~446/726/734; `Visualization.py:184`
rank_genes + dendrogram) with proper decorator caching keyed on parameters.

### 4. `st.bottom`
Pin the hottest of the 31 scattered `st.download_button` / apply controls into a
persistent footer so users don't scroll up to re-download.

### 5. `type` param on `st.expander` / `st.status` (new in 1.58)
Cheap visual polish — a more compact style for CAFE's expanders/status blocks
(e.g. `Data_Processing.py:219, 257`).

## Roadmap: the numpy-2 migration (unlocks the rest)

To modernize the scientific stack (scanpy 1.12, anndata 0.13, plotly 6, kaleido 1) the team
must migrate to `numpy>=2`. Concrete work involved:

1. Bump `numpy>=2`, `pandas>=2.3`, then `scanpy` 1.12.x, `anndata` 0.13.x.
2. Rework all `fig.write_image(..., engine='kaleido')` call sites for kaleido v1
   (the `engine=`/`scope` API is removed).
3. Re-validate `umap-learn` / `numba` under numpy 2.
4. Full manual re-run of the QC → PCA → UMAP → Leiden → Visualization pipeline plus all
   export paths (no automated CI covers this).
