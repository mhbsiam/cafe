"""Page-level (AppTest) checks for Data_Processing.py: new uploads start at QC, and use_rep "None" means raw X."""
import io
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from streamlit.testing.v1 import AppTest

CAFE_APP = Path(__file__).resolve().parents[1] / "src" / "cafe_app"
PAGE = CAFE_APP / "tools" / "Data_Processing.py"
sys.path.insert(0, str(CAFE_APP))  # for `theme`; conftest already adds tools/
matplotlib.use("Agg")  # AppTest runs the page off the main thread; the macOS GUI backend refuses that

import session_store  # noqa: E402
import theme  # noqa: E402
import utils  # noqa: E402

_STALE_FLAGS = [
    "qc_done", "pca_done", "pca_selected", "batch_correction_done",
    "umap_computed", "leiden_computed",
]


@pytest.fixture(autouse=True)
def _temp_store(monkeypatch, tmp_path):
    for mod in (session_store, theme, utils):
        assert Path(mod.__file__).resolve().is_relative_to(CAFE_APP), mod.__file__
    base = tmp_path / "sessions"
    base.mkdir(mode=0o700)
    monkeypatch.setattr(session_store, "store_base", lambda: base)


def _csv(name, text):
    f = io.BytesIO(text.encode())
    f.name = name
    return f


def test_new_upload_resets_stale_flags_and_shows_load_notice_once():
    rows = "".join(f"{i}.0,{i + 1}.0,\n" for i in range(20))
    files = [
        _csv("S1_Young.csv", "CD3,CD4,Empty\n" + rows),  # all-empty column -> load notice
        _csv("S2_Aged.csv", "CD3,CD4\n" + rows.replace(",\n", "\n")),
    ]
    review = utils.inspect_uploaded_csv_files(files)
    assert not review["bad_names"] and not review["mismatches"]

    at = AppTest.from_file(str(PAGE), default_timeout=60)
    for key in _STALE_FLAGS:
        at.session_state[key] = True
    at.session_state["_umap_pending"] = [("success", "stale")]
    at.session_state["_load_notice"] = "stale notice from an earlier load"
    at.session_state["qc_nmads_S1"] = 7.0
    at.session_state["_umap_fig_leiden"] = b"stale-png"
    at.session_state["_upload_files"] = files
    at.session_state["_upload_review"] = review
    at.run()
    assert not at.exception
    at.button(key="dp_upload_proceed").click().run()

    assert not at.exception
    for key in _STALE_FLAGS:
        assert at.session_state[key] is False, key
    assert "_umap_pending" not in at.session_state
    assert "_umap_fig_leiden" not in at.session_state
    assert "qc_nmads_S1" not in at.session_state
    md = " ".join(m.value for m in at.markdown)
    assert "Quality Control" in md and "UMAP + Leiden Clustering" not in md
    warnings = [w.value for w in at.warning]
    assert any("empty column" in w for w in warnings), warnings
    assert not any("stale notice" in w for w in warnings)
    assert "_load_notice" not in at.session_state

    at.run()  # shown once only
    assert not any("empty column" in w.value for w in at.warning)


def test_duplicate_samples_block_proceed_with_own_card():
    rows = "CD3,CD4\n1.0,2.0\n"
    files = [_csv("S1_Young.csv", rows), _csv("S1_Aged.csv", rows)]
    at = AppTest.from_file(str(PAGE), default_timeout=60)
    at.session_state["_upload_files"] = files
    at.session_state["_upload_review"] = utils.inspect_uploaded_csv_files(files)
    at.run()
    assert not at.exception
    assert at.button(key="dp_upload_proceed").disabled
    html = " ".join(m.value for m in at.markdown)
    assert "Duplicate samples" in html
    assert "Filename convention" not in html


def _click(at, label):
    next(b for b in at.button if b.label == label).click().run()


def _clustering_ready_page(n=600, timeout=180):
    """The page at the UMAP + Leiden step, with four separable populations loaded."""
    rng = np.random.default_rng(0)
    centers = rng.normal(0, 6, (4, 8))
    adata = sc.AnnData(
        centers[rng.integers(0, 4, n)] + rng.normal(0, 1, (n, 8)),
        obs=pd.DataFrame(
            {"SampleID": np.repeat(["S1", "S2", "S3"], n // 3), "Group": np.repeat(["A", "B", "A"], n // 3)},
            index=[str(i) for i in range(n)],
        ),
    )
    adata.var_names = [f"M{i}" for i in range(8)]
    at = AppTest.from_file(str(PAGE), default_timeout=timeout)
    at.session_state["adata"] = adata
    at.session_state["_dp_started"] = True
    for key in ("qc_done", "pca_done", "batch_correction_done"):
        at.session_state[key] = True
    return at.run()


def test_stability_sweep_suggests_a_resolution_and_feeds_the_report():
    at = _clustering_ready_page()
    _click(at, "Compute graph and UMAP")
    next(t for t in at.text_input if t.label.startswith("Resolutions")).set_value("0.1, 0.5")
    next(n for n in at.number_input if n.label == "Repeats").set_value(3)
    _click(at, "Run stability check")
    assert not at.exception
    stab = at.session_state["_params_stability"]
    assert list(stab["table"]["resolution"]) == [0.1, 0.5] and stab["subsamples"] == 3 and stab["seeds"] == 0
    assert next(s for s in at.slider if s.label == "Leiden Resolution").value == stab["suggested"]

    _click(at, "Run Leiden clustering")
    assert not at.exception
    result = at.session_state["adata"]
    scores = at.session_state["_cluster_stability"]
    assert set(scores["cluster"]) == set(result.obs["leiden"].astype(str))
    assert result.obs["leiden_stability"].between(0, 1).all()
    at.run()  # the results page packages the download
    import zipfile
    names = zipfile.ZipFile(io.BytesIO(at.session_state["_last_zip_buffer"])).namelist()
    assert any(n.startswith("cluster_stability_") for n in names)
    assert any(n.startswith("resolution_stability_") for n in names)
    report = zipfile.ZipFile(io.BytesIO(at.session_state["_last_zip_buffer"])).read(
        next(n for n in names if n.startswith("analysis_parameters"))
    ).decode()
    assert "RESOLUTION STABILITY" in report and "Clustering stability across Leiden resolutions" in report


def test_change_resolution_keeps_the_graph_and_drops_old_scores():
    at = _clustering_ready_page()
    _click(at, "Compute graph and UMAP")
    umap_before = at.session_state["adata"].obsm["X_umap"].copy()
    _click(at, "Run Leiden clustering")
    assert "_cluster_stability" in at.session_state
    _click(at, "Change resolution")
    assert not at.exception
    assert "_cluster_stability" not in at.session_state and not at.session_state["leiden_computed"]
    next(c for c in at.checkbox if c.label.startswith("Score each cluster")).uncheck()
    next(s for s in at.slider if s.label == "Leiden Resolution").set_value(1.0)
    _click(at, "Run Leiden clustering")
    assert not at.exception
    result = at.session_state["adata"]
    assert "leiden_stability" not in result.obs  # no stale scores from the first resolution
    assert at.session_state["_params_umap"]["resolution"] == 1.0
    assert (result.obsm["X_umap"] == umap_before).all()

    _click(at, "Change resolution")
    _click(at, "Change graph settings")
    _click(at, "Compute graph and UMAP")
    assert not at.exception
    assert "leiden" not in at.session_state["adata"].obs  # clusters from the old graph are not kept


def test_none_representation_uses_raw_x_without_pca():
    rng = np.random.default_rng(0)
    n, m = 300, 60  # > 50 markers: scanpy would use/compute X_pca if use_rep were None
    adata = sc.AnnData(
        rng.normal(size=(n, m)).astype("float32"),
        obs=pd.DataFrame(
            {"SampleID": np.repeat(["S1", "S2", "S3"], n // 3), "Group": np.repeat(["A", "B", "A"], n // 3)},
            index=[str(i) for i in range(n)],
        ),
    )
    adata.var_names = [f"M{i}" for i in range(m)]

    at = AppTest.from_file(str(PAGE), default_timeout=180)
    at.session_state["adata"] = adata
    at.session_state["_dp_started"] = True
    for key in ("qc_done", "pca_done", "batch_correction_done"):
        at.session_state[key] = True
    at.run()
    assert not at.exception
    assert at.selectbox[0].value == "None"
    assert not any(s.label == "Leiden flavor" for s in at.selectbox)
    _click(at, "Compute graph and UMAP")
    assert not at.exception
    assert "leiden" not in at.session_state["adata"].obs  # graph and UMAP only; Leiden is its own part
    _click(at, "Run Leiden clustering")

    assert not at.exception
    result = at.session_state["adata"]
    assert "leiden" in result.obs
    assert "X_pca" not in result.obsm
    assert result.uns["neighbors"]["params"].get("use_rep") == "X"
    params = at.session_state["_params_umap"]
    assert params["use_rep"] == "raw expression (X)"
    assert params["leiden_flavor"] == "igraph"
