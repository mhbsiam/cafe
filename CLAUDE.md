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
# Easiest: install as a global command (pipx / uv), then run from anywhere
pipx install "git+https://github.com/mhbsiam/cafe"   # or: pipx install .
cafe

cd CAFE

# Dev: using Pixi (manages its own isolated environment)
pixi run cafe

# Dev alternative: direct Python (requires conda env activated)
conda activate cafe
python cafe.py
```

The app is a proper installable package (`pyproject.toml` at the **repo root**,
package source under `CAFE/src/cafe_app/`). Installing it exposes three console
scripts:
- `cafe` → `cafe_app.launcher:main` — standard desktop use (3 GB upload limit)
- `cafe-web` → `cafe_app.launcher:main_web` — web-hosted deployment (2 GB upload limit)
- `cafe-hpc` → `cafe_app.hpc:main` — HPC cluster CLI

The desktop/web launchers delegate to `src/cafe_app/run.py` / `run_w.py`, which
define the Streamlit multi-page navigation. The root `CAFE/cafe.py` and `web.py`
are thin dev wrappers (used by `pixi run`) that call the same launcher.

Dependencies are declared in the root `pyproject.toml`, kept in sync with the
`[pypi-dependencies]` in `CAFE/pixi.toml`.

### Environment Setup

```bash
cd CAFE

# Pixi (recommended)
pixi install

# OR Conda
conda env create -f cafe.yaml
conda activate cafe
```

### Tests

A pytest suite lives in `CAFE/tests/`. Run it with:

```bash
cd CAFE && pixi run test    # or: pytest tests/ -v
```

`tests/conftest.py` puts `src/cafe_app/tools/` on `sys.path` so tests can
`from utils import ...` the same way the app pages do.

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

1. `cafe_app.launcher:main` (or dev wrapper `cafe.py`) → sets Streamlit server config → delegates to `run.py`
2. `run.py` → defines navigation pages using `st.navigation`
3. Pages load independently; each is a self-contained Streamlit page module

Page modules are still loaded by relative path (`st.Page("app/CAFE.py")`) and
import shared helpers as top-level modules (`from theme import ...`,
`from utils import ...`) — `run.py` puts `tools/` on `sys.path`. Image assets
are resolved via `theme.IMG_DIR` (absolute, `__file__`-relative) so they work
regardless of the current working directory.

### Key Modules (`CAFE/src/cafe_app/`)

| Path | Role |
|------|------|
| `launcher.py` | Console entry points (`cafe`, `cafe-web`); bakes in theme config |
| `run.py` / `run_w.py` | Streamlit multi-page navigation (desktop / web) |
| `hpc.py` | HPC command-line pipeline (`cafe-hpc`) |
| `app/CAFE.py` | Welcome/landing page (desktop build) |
| `app/landing_content.py` | Shared landing page content used by desktop and web builds |
| `web/CAFE.py` | Welcome/landing page (hosted build) |
| `theme.py` | Shared design tokens, CSS, page-header helpers, and `IMG_DIR` |
| `tools/Data_Processing.py` | CSV ingestion, preprocessing, session state |
| `tools/Visualization.py` | Core analysis: UMAP, clustering, all plots (~4300 lines) |
| `tools/_plot_export.py` | Deferred download button for large figures |
| `tools/utils.py` | Shared CSV, batch correction, and statistical helpers |
| `advanced/Merge_Clusters.py` | Post-clustering merge operations |
| `advanced/Cluster_Annotation.py` | Manual cluster labeling |
| `advanced/Cluster_Evaluation.py` | Cluster quality metrics |
| `advanced/Selective_Clustering.py` | Subset re-clustering |
| `advanced/Downsampling.py` | Data downsampling utilities |
| `img/` | Bundled image assets (logo, workflow diagrams) |

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
- The theme (teal `#0f5070` primary color, Source Sans 3 font, shared design tokens) is set by `CAFE/src/cafe_app/theme.py`; `CAFE/.streamlit/config.toml` applies for `pixi run`/`conda` launches, and the same values are passed as config flags by `launcher.py` for the installed `cafe` command
