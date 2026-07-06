import os
import sys

import streamlit as st

# Streamlit only puts this entry script's directory (the package root) on
# sys.path. The shared utils/_plot_export modules live in ``tools/``, so add
# that once here — every page can then `from utils import ...` without
# repeating the path shim.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools'))

from theme import apply_theme

# A2: deployment-mode flag — 'web' for hosted deployment
st.session_state.setdefault('cafe_deployment', 'web')

page1 = st.Page("web/CAFE.py", default=True)
page2 = st.Page("tools/Data_Processing.py")
page3 = st.Page("tools/Visualization.py")

page4 = st.Page("advanced/Merge_Clusters.py")
page5 = st.Page("advanced/Cluster_Annotation.py")

page6 = st.Page("advanced/Cluster_Evaluation.py")
page7 = st.Page("advanced/Selective_Clustering.py")
page8 = st.Page("advanced/Downsampling.py")

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
