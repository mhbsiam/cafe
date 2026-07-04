import os
import sys
import warnings

# scanpy 1.10.4 still reads anndata's now-deprecated `__version__` attribute
# (anndata >= 0.12), which emits a noisy FutureWarning three times on every
# launch. Silence just that message — not all FutureWarnings — until scanpy
# ships a fix. Must run before any page module imports scanpy.
warnings.filterwarnings(
    "ignore",
    message=r"`__version__` is deprecated",
    category=FutureWarning,
)

import streamlit as st

# Streamlit only puts this entry script's directory (bin/) on sys.path. The
# shared utils/_plot_export modules live in bin/Tools, so add that once here —
# every page can then `from utils import ...` without repeating the path shim.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Tools'))

from theme import apply_theme

# A2: deployment-mode flag — 'desktop' for local app, 'web' for hosted deployment
st.session_state.setdefault('cafe_deployment', 'desktop')

page1 = st.Page("App/CAFE.py", default=True)
page2 = st.Page("Tools/Data_Processing.py")
page3 = st.Page("Tools/Visualization.py")

page4 = st.Page("Adv/Merge_Clusters.py")
page5 = st.Page("Adv/Cluster_Annotation.py")

page6 = st.Page("Adv/Cluster_Evaluation.py")
page7 = st.Page("Adv/Selective_Clustering.py")
page8 = st.Page("Adv/Downsampling.py")

pg = st.navigation(
    {
        "Welcome": [page1],
        "Tools": [page2, page3],
        "Advanced": [page4, page5, page6, page7, page8],
    }
)

# Track page transitions so a page can reset transient UI state on a fresh
# visit. st.session_state persists across navigation, so a page can't otherwise
# tell "the user just navigated here" apart from "my own widget triggered a
# rerun" — both re-run this script. This exposes `_page_changed` (True only on
# the first run after arriving at a different page) for pages that need it.
st.session_state['_page_changed'] = st.session_state.get('_active_page') != pg.url_path
st.session_state['_active_page'] = pg.url_path

pg.run()

apply_theme()
