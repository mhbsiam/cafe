import os
import sys
import warnings

# Silence scanpy's noisy anndata `__version__` FutureWarning; must run before scanpy is imported.
warnings.filterwarnings(
    "ignore",
    message=r"`__version__` is deprecated",
    category=FutureWarning,
)

import streamlit as st

# Put tools/ on sys.path so pages can `from utils import ...` without a per-page shim.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools'))

from theme import apply_theme
from session_store import ensure_token, restore_session_adata, render_clear_control

# Refresh-persistence: mint a per-session token (survives refresh via the URL) and, if the browser was
# refreshed, restore the loaded AnnData from disk before any page runs.
ensure_token()
restore_session_adata()

page1 = st.Page("app/CAFE.py", default=True)
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

# Expose _page_changed (True only on the first run after navigating to a different page).
st.session_state['_page_changed'] = st.session_state.get('_active_page') != pg.url_path
st.session_state['_active_page'] = pg.url_path

# Render the sidebar Clear control before pg.run() so pages that call st.stop() (e.g. every Data Processing
# step) can't skip it.
render_clear_control()

pg.run()

apply_theme()
