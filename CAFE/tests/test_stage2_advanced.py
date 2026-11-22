"""Stage 2 fixes in the advanced pages and celltype_gmm (posterior monotonicity, sorting, merge validation, downsampling)."""
import ast
import importlib.util
import os

import numpy as np
import pandas as pd
import pytest

CAFE_APP = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'cafe_app')
ADVANCED = os.path.join(CAFE_APP, 'advanced')
WORKTREE = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import celltype_gmm  # noqa: E402  (conftest puts this worktree's tools/ on sys.path)
from utils import safe_sort_clusters  # noqa: E402


def _assert_in_worktree(module):
    assert os.path.realpath(module.__file__).startswith(WORKTREE)


def _load_functions(page, names):
    """Exec only the named top-level functions of a Streamlit page (the page body needs a running app)."""
    path = os.path.join(ADVANCED, page)
    tree = ast.parse(open(path).read())
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    ns = {'__file__': path}
    exec('import os, time\nimport numpy as np\nimport pandas as pd\n'
         'from joblib import Parallel, delayed\nfrom scipy.stats import gaussian_kde\n'
         'from sklearn.decomposition import PCA', ns)
    exec(compile(tree, path, 'exec'), ns)
    return ns


# ── G1: monotone GMM posterior ────────────────────────────────────────────────

def test_gmm_posterior_is_zero_below_negative_mean_and_one_above_positive_mean():
    _assert_in_worktree(celltype_gmm)
    rng = np.random.default_rng(0)
    # Negative N(0, 0.3), wider positive N(3, 1.5), equal weights; plus probes at -2 and 6.
    x = np.concatenate([rng.normal(0, 0.3, 5000), rng.normal(3, 1.5, 5000), [-2.0, 6.0]])
    proba, info = celltype_gmm.fit_marker_gmm(x)

    gmm, pos = info["model"]["gmm"], info["model"]["pos"]
    raw = gmm.predict_proba(np.array([[-2.0]]))[0, pos]
    assert raw > 0.5  # the unclamped posterior calls this very dim cell positive (the bug)

    assert proba[-2] == 0.0
    assert proba[-1] == 1.0
    assert 0 < info["cutoff"] < 3
    neg_mean, pos_mean = info["means"]
    assert neg_mean < info["cutoff"] < pos_mean

    # The pooled-scoring path agrees with the per-fit path.
    np.testing.assert_array_equal(celltype_gmm._score_with_model(x, info["model"]), proba)

    # Monotone non-decreasing in x.
    grid = np.linspace(-4, 8, 2001)
    assert np.all(np.diff(celltype_gmm._score_with_model(grid, info["model"])) >= -1e-12)


# ── M2: mixed cluster labels ──────────────────────────────────────────────────

def test_safe_sort_clusters_mixed_labels():
    assert safe_sort_clusters(["10", "Tcells", "2", "0"]) == ["0", "10", "2", "Tcells"]
    assert safe_sort_clusters(["10", "2", "0"]) == ["0", "2", "10"]


def test_cluster_annotation_uses_safe_sort_not_isdigit():
    src = open(os.path.join(ADVANCED, 'Cluster_Annotation.py')).read()
    assert 'isdigit' not in src
    assert src.count('safe_sort_clusters(') >= 3  # two fixed sites + the existing purity table


# ── M3: merge target colliding with an existing cluster ───────────────────────

def test_merge_target_equal_to_other_existing_cluster_is_rejected():
    ns = _load_functions('Merge_Clusters.py', {'_validate_merge_groups'})
    validate = ns['_validate_merge_groups']
    clusters = ["0", "1", "2"]

    errors = validate([{"target": "1", "sources": ["0", "2"]}], clusters)
    assert any("Pick another name" in e for e in errors)

    # Renaming into one of the group's own sources, or to a fresh name, stays valid.
    assert validate([{"target": "0", "sources": ["0", "2"]}], clusters) == []
    assert validate([{"target": "T cells", "sources": ["0", "2"]}], clusters) == []


def test_merge_drops_stability_scores_of_the_old_clusters():
    src = open(os.path.join(ADVANCED, 'Merge_Clusters.py')).read()
    assert "drop(columns='leiden_stability'" in src


# ── D2-D5: downsampling ───────────────────────────────────────────────────────

@pytest.fixture
def downsampling():
    return _load_functions('Downsampling.py', {
        'run_pca', 'kde_fft', 'downsample_uniform', 'downsample_density_weighted',
        'downsample_inverse_density', 'process_file', 'process_directory_parallel',
    })


def _write_inputs(folder):
    rng = np.random.default_rng(0)
    for name in ("S1_Ctrl.csv", "S2_Ctrl.CSV"):
        pd.DataFrame(rng.normal(size=(300, 3)), columns=["CD3", "CD4", "CD8"]).to_csv(
            os.path.join(folder, name), index=False)


def test_downsampling_keeps_names_and_reads_uppercase_csv(tmp_path, downsampling):
    src, out = tmp_path / "in", tmp_path / "out"
    src.mkdir()
    out.mkdir()
    _write_inputs(src)

    _, skipped = downsampling['process_directory_parallel'](
        str(src), str(out), "uniform", 5, 20, 100, 1)

    assert sorted(os.listdir(out)) == ["S1_Ctrl.csv", "S2_Ctrl.CSV"]
    assert skipped == []
    assert len(pd.read_csv(out / "S2_Ctrl.CSV")) == 100


def test_downsampling_refuses_same_folder(tmp_path, downsampling):
    _write_inputs(tmp_path)
    before = {f: (tmp_path / f).read_bytes() for f in os.listdir(tmp_path)}
    with pytest.raises(ValueError, match="different from the input"):
        downsampling['process_directory_parallel'](
            str(tmp_path), str(tmp_path / "."), "uniform", 5, 20, 100, 1)
    assert {f: (tmp_path / f).read_bytes() for f in os.listdir(tmp_path)} == before


def test_downsampling_errors_when_no_csv(tmp_path, downsampling):
    src, out = tmp_path / "in", tmp_path / "out"
    src.mkdir()
    out.mkdir()
    (src / "notes.txt").write_text("x")
    with pytest.raises(ValueError, match="No CSV files"):
        downsampling['process_directory_parallel'](str(src), str(out), "uniform", 5, 20, 100, 1)


def test_downsampling_jobs_widget_offers_only_valid_values():
    from streamlit.testing.v1 import AppTest

    script = f"""
import sys
sys.path[:0] = [{os.path.realpath(CAFE_APP)!r}, {os.path.realpath(os.path.join(CAFE_APP, 'tools'))!r}]
import streamlit as st
st.session_state['downsampling_configured'] = True
exec(compile(open({os.path.join(ADVANCED, 'Downsampling.py')!r}).read(), 'Downsampling.py', 'exec'))
"""
    at = AppTest.from_string(script, default_timeout=60).run()
    assert not at.exception
    jobs = [s for s in at.select_slider if 'number of jobs' in s.label]
    assert len(jobs) == 1
    assert list(jobs[0].options) == ["-1", "1", "2", "3", "4"]
    assert jobs[0].value == -1


# ── S2: Selective Clustering "Apply" resets Data Processing state too ─────────

def test_selective_apply_clears_data_processing_flags(tmp_path, monkeypatch):
    import session_store
    from streamlit.testing.v1 import AppTest

    _assert_in_worktree(session_store)
    monkeypatch.setattr(session_store, "store_base", lambda: tmp_path)
    script = f"""
import sys
sys.path[:0] = [{os.path.realpath(CAFE_APP)!r}, {os.path.realpath(os.path.join(CAFE_APP, 'tools'))!r}]
import matplotlib
matplotlib.use("Agg")
exec(compile(open({os.path.join(ADVANCED, 'Selective_Clustering.py')!r}).read(), 'Selective_Clustering.py', 'exec'))
"""
    at = AppTest.from_string(script, default_timeout=60)
    at.session_state['qc_done'] = True
    at.session_state['_umap_pending'] = True
    at.session_state['pca_done'] = True
    at.run()
    assert not at.exception
    at.button[0].click().run()
    assert not at.exception
    for key in ('qc_done', '_umap_pending', 'pca_done'):
        assert key not in at.session_state
    assert at.session_state['sel_pca_done'] is False
    assert at.session_state['adata'] is None
