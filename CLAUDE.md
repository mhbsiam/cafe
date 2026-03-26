# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Layout

This repo has two independent components:

- **`CAFE/`** — The Python/Streamlit application (the actual tool)
- **`src/`** — Astro/Starlight documentation site (deployed to GitHub Pages)

---

## CAFE Application

### Running the App

```bash
cd CAFE

# Preferred: using Pixi (manages its own isolated environment)
pixi run cafe

# Alternative: direct Python (requires conda env activated)
conda activate cafe
streamlit run cafe.py
```

There are three entry point variants:
- `cafe.py` — standard desktop use (3 GB upload limit)
- `web.py` — web-hosted deployment (2 GB upload limit)
- `cafe_hpc.py` — HPC cluster use

All three delegate to `bin/run.py` which defines the Streamlit multi-page navigation.

### Environment Setup

```bash
cd CAFE

# Pixi (recommended)
pixi install

# OR Conda
conda env create -f cafe.yaml
conda activate cafe
```

### No Test Suite

There are no automated tests. Validation is done manually by running the app.

---

## Documentation Site

```bash
# Install dependencies
pnpm install

# Dev server
pnpm dev

# Build static site
pnpm build
```

Deployed automatically to `https://mhbsiam.github.io/cafe` via `.github/workflows/deploy.yml` on push to `main`.

---

## Architecture

CAFE is a **no-code bioinformatics web app** for analyzing spectral flow cytometry (SFCM) data. Users upload CSV files exported from FlowJo, run dimensionality reduction / clustering, and download publication-ready figures.

### Application Flow

1. `cafe.py` → sets Streamlit server config → delegates to `bin/run.py`
2. `bin/run.py` → defines navigation pages using `st.navigation`
3. Pages load independently; each is a self-contained Streamlit page module

### Key Modules (`CAFE/bin/`)

| Path | Role |
|------|------|
| `App/CAFE.py` | Welcome/landing page, CSV format instructions |
| `Tools/Data_Processing.py` | CSV ingestion, preprocessing, session state |
| `Tools/Visualization.py` | Core analysis: UMAP, clustering, all plots (~4300 lines) |
| `Adv/Merge_Clusters.py` | Post-clustering merge operations |
| `Adv/Cluster_Annotation.py` | Manual cluster labeling |
| `Adv/Cluster_Evaluation.py` | Cluster quality metrics |
| `Adv/Selective_Clustering.py` | Subset re-clustering |
| `Adv/Downsampling.py` | Data downsampling utilities |

### Core Stack

- **Streamlit** — entire UI layer; session state (`st.session_state`) is the primary data store between pages
- **Scanpy / AnnData** — single-cell analysis backbone; data is stored as AnnData objects
- **Pandas / NumPy / SciPy** — data processing
- **Matplotlib / Seaborn / Plotly** — figure generation (Kaleido for static export)
- **igraph** — graph-based clustering (Leiden/Louvain)

### Data Flow

CSV upload → Pandas DataFrame → AnnData object (stored in `st.session_state`) → analysis steps mutate the AnnData → figures rendered inline → ZIP download of results

### Streamlit Conventions

- All inter-page state lives in `st.session_state`
- Each page checks for required session keys at the top and shows an error if prerequisites are missing (user must run Data Processing before Visualization)
- The `.streamlit/config.toml` sets the theme (red primary color, sans-serif font)
