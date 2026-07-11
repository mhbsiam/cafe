"""Tests for Data_Processing.py bug fixes.

Covers: B1 (selected_flavor), B2 (selected_dim), B11 (duplicate columns),
B12 (PCA variance threshold).

These tests verify the logic of the fixes without requiring a full Streamlit
runtime. They test the underlying scanpy/anndata calls directly.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import pytest

from utils import (
    safe_sort_clusters,
    compute_batch_emd,
    summarize_batch_effect,
    recommend_nmads,
)
from scipy.stats import median_abs_deviation

# Check if leidenalg is available (it's an optional dependency)
try:
    import leidenalg  # noqa: F401
    HAS_LEIDENALG = True
except ImportError:
    HAS_LEIDENALG = False


class TestB1LeidenFlavor:
    """B1: selected_flavor must be passed to sc.tl.leiden, not hardcoded."""

    def test_leiden_igraph_flavor(self, small_dense_adata):
        """Verify leiden with flavor='igraph' produces clusters."""
        adata = small_dense_adata.copy()
        sc.pp.neighbors(adata, n_neighbors=15, use_rep=None, random_state=50)
        sc.tl.leiden(adata, resolution=0.5, random_state=50, flavor="igraph",
                     n_iterations=2, directed=False)
        assert 'leiden' in adata.obs.columns
        assert adata.obs['leiden'].nunique() > 0

    def test_leiden_leidenalg_flavor(self, small_dense_adata):
        """Verify leiden with flavor='leidenalg' works (the previously ignored option)."""
        if not HAS_LEIDENALG:
            pytest.skip("leidenalg not installed in this environment")
        adata = small_dense_adata.copy()
        sc.pp.neighbors(adata, n_neighbors=15, use_rep=None, random_state=50)
        # This is the call that was previously hardcoded to "igraph"
        sc.tl.leiden(adata, resolution=0.5, random_state=50, flavor="leidenalg",
                     n_iterations=2, directed=False)
        assert 'leiden' in adata.obs.columns
        assert adata.obs['leiden'].nunique() > 0

    def test_flavors_produce_results(self, small_dense_adata):
        """Both flavors should produce valid clustering results."""
        flavors = ['igraph']
        if HAS_LEIDENALG:
            flavors.append('leidenalg')
        for flavor in flavors:
            adata = small_dense_adata.copy()
            sc.pp.neighbors(adata, n_neighbors=15, use_rep=None, random_state=50)
            sc.tl.leiden(adata, resolution=0.5, random_state=50, flavor=flavor,
                         n_iterations=2, directed=False)
            assert adata.obs['leiden'].dtype == 'category'


class TestB2SelectedDim:
    """B2: selected_dim must be used as use_rep, not just pca_choice boolean."""

    def test_use_rep_none(self, small_dense_adata):
        """When selected_dim='None', use_rep should be None."""
        adata = small_dense_adata.copy()
        use_rep = None  # This is what the fix produces for 'None'
        sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep, random_state=50)
        assert 'neighbors' in adata.uns

    def test_use_rep_x_pca(self, small_dense_adata):
        """When selected_dim='X_pca', use_rep should be 'X_pca'."""
        adata = small_dense_adata.copy()
        sc.tl.pca(adata, n_comps=3, svd_solver='full')
        use_rep = 'X_pca'
        sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep, random_state=50)
        assert 'neighbors' in adata.uns

    def test_use_rep_x_umap(self, small_dense_adata):
        """When selected_dim='X_umap', use_rep should be 'X_umap'."""
        adata = small_dense_adata.copy()
        sc.tl.pca(adata, n_comps=3, svd_solver='full')
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca', random_state=50)
        sc.tl.umap(adata, random_state=50)
        use_rep = 'X_umap'
        sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep, random_state=50)
        assert 'neighbors' in adata.uns

    def test_invalid_use_rep_raises(self, small_dense_adata):
        """If selected_dim references a non-existent representation, should error."""
        adata = small_dense_adata.copy()
        use_rep = 'X_umap'
        # X_umap doesn't exist yet
        assert use_rep not in adata.obsm
        # Actually exercise the error path: neighbors on a missing representation must raise.
        with pytest.raises((KeyError, ValueError)):
            sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep, random_state=50)


class TestB11DuplicateColumns:
    """B11: CSV with existing SampleID/Group columns must not produce duplicates."""

    def test_no_duplicate_columns(self):
        """Verify that the B11 fix logic drops existing columns before appending."""
        import pyarrow as pa
        import pyarrow.csv as pv
        import io as io_module

        # Create a CSV with existing SampleID and Group columns
        csv_content = "CD1,CD2,SampleID,Group\n1.0,2.0,existing_S1,existing_G1\n3.0,4.0,existing_S2,existing_G2\n"
        table = pv.read_csv(io_module.BytesIO(csv_content.encode()))

        # Simulate the B11 fix: drop existing SampleID/Group before appending
        if 'SampleID' in table.column_names:
            table = table.drop(['SampleID'])
        if 'Group' in table.column_names:
            table = table.drop(['Group'])

        # Now append filename-derived columns
        table = table.append_column('SampleID', pa.array(['S1', 'S1']))
        table = table.append_column('Group', pa.array(['G1', 'G1']))

        # Verify no duplicates
        col_names = table.column_names
        assert col_names.count('SampleID') == 1
        assert col_names.count('Group') == 1
        # Verify the values are from the filename, not the original CSV
        assert table.column('SampleID').to_pylist() == ['S1', 'S1']
        assert table.column('Group').to_pylist() == ['G1', 'G1']


class TestB12PCAThreshold:
    """B12: PCA variance threshold edge case when no component reaches threshold."""

    def test_threshold_reached(self, small_dense_adata):
        """Normal case: threshold is reached by some component."""
        adata = small_dense_adata.copy()
        sc.tl.pca(adata, n_comps=3, svd_solver='full')
        variance_threshold = 50  # 50% should be reachable
        pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()
        threshold_mask = pca_variance_ratio >= variance_threshold / 100
        assert threshold_mask.any()
        num_components = threshold_mask.argmax() + 1
        assert num_components >= 1

    def test_threshold_not_reached(self, small_dense_adata):
        """Edge case: threshold too high, no component reaches it."""
        adata = small_dense_adata.copy()
        sc.tl.pca(adata, n_comps=3, svd_solver='full')
        variance_threshold = 99.99  # Almost certainly unreachable with random data
        pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()
        threshold_mask = pca_variance_ratio >= variance_threshold / 100

        # B12 fix: when no component reaches threshold, use all components
        if not threshold_mask.any():
            num_components = len(pca_variance_ratio)
        else:
            num_components = threshold_mask.argmax() + 1

        # Should use all components, not just 1 (which was the old buggy behavior)
        assert num_components == len(pca_variance_ratio)
        assert num_components > 1  # Old code would have returned 1


class TestBatchEMD:
    """Batch effect diagnostics: compute_batch_emd + summarize_batch_effect."""

    @staticmethod
    def _make_adata(shift_marker=None, shift=0.0, confounded=False, seed=0):
        """Two batches x two groups; optionally shift one marker in batch '2'.

        When ``confounded`` is True, Group is nested within Batch (batch 1 is all
        Healthy, batch 2 is all Disease) so the confounding guard should fire.
        """
        rng = np.random.default_rng(seed)
        n_per = 200
        n_markers = 5
        marker_names = [f"CD{i}" for i in range(n_markers)]

        # 4 samples, 2 per batch.
        batch = np.array(["1"] * (2 * n_per) + ["2"] * (2 * n_per))
        sample = np.array(
            ["S1"] * n_per + ["S2"] * n_per + ["S3"] * n_per + ["S4"] * n_per
        )
        if confounded:
            group = np.where(batch == "1", "Healthy", "Disease")
        else:
            # Group crossed with batch: each batch spans both groups.
            group = np.array(
                (["Healthy"] * n_per + ["Disease"] * n_per) * 2
            )

        X = rng.normal(0.0, 1.0, size=(4 * n_per, n_markers))
        if shift_marker is not None:
            j = marker_names.index(shift_marker)
            X[batch == "2", j] += shift

        obs = pd.DataFrame(
            {"SampleID": sample, "Group": group, "Batch": batch},
            index=[str(i) for i in range(4 * n_per)],
        )
        adata = sc.AnnData(X, obs=obs)
        adata.var_names = marker_names
        return adata

    def test_shifted_marker_ranks_top(self):
        """A marker shifted in one batch should dominate the EMD ranking."""
        adata = self._make_adata(shift_marker="CD2", shift=5.0)
        emd = compute_batch_emd(adata, batch_key="Batch")
        assert list(emd.columns) == ["marker", "emd_max", "emd_mean", "worst_batch"]
        assert emd.iloc[0]["marker"] == "CD2"
        top = emd.set_index("marker")["emd_max"]
        # Shifted marker's EMD should dwarf every other marker.
        others = top.drop("CD2").max()
        assert top["CD2"] > 5 * others
        # worst_batch is the batch farthest from the pooled reference. With a
        # symmetric two-batch shift both batches are ~equidistant from the pool,
        # so the winner is a floating-point tie that flips across BLAS backends
        # (Accelerate vs OpenBLAS) and numpy/scipy versions. Assert only that it
        # names a real batch; the dominance checks above carry the actual intent.
        assert emd.iloc[0]["worst_batch"] in {"1", "2"}

    def test_well_mixed_batches(self):
        """No shift -> tiny EMDs, success verdict, nothing flagged."""
        adata = self._make_adata()
        emd = compute_batch_emd(adata, batch_key="Batch")
        assert (emd["emd_max"] < 0.5).all()
        verdict, kind, n_flagged, confounded = summarize_batch_effect(emd, adata)
        assert kind == "success"
        assert n_flagged == 0
        assert confounded == []

    def test_shifted_marker_is_flagged(self):
        """A real shift should be flagged and produce a warning verdict."""
        adata = self._make_adata(shift_marker="CD2", shift=5.0)
        emd = compute_batch_emd(adata, batch_key="Batch")
        verdict, kind, n_flagged, confounded = summarize_batch_effect(emd, adata)
        assert kind == "warning"
        assert n_flagged >= 1
        assert confounded == []

    def test_confounded_batch_reported(self):
        """Batch nested within Group -> confounding guard fires."""
        adata = self._make_adata(shift_marker="CD2", shift=5.0, confounded=True)
        emd = compute_batch_emd(adata, batch_key="Batch")
        verdict, kind, n_flagged, confounded = summarize_batch_effect(emd, adata)
        assert confounded == ["Group"]
        assert kind == "warning"
        assert "confounded" in verdict.lower()

    def test_single_batch_returns_empty(self):
        """Fewer than two batches -> empty frame and an info verdict."""
        adata = self._make_adata()
        adata.obs["Batch"] = "1"
        emd = compute_batch_emd(adata, batch_key="Batch")
        assert emd.empty
        verdict, kind, n_flagged, confounded = summarize_batch_effect(emd, adata)
        assert kind == "info"
        assert n_flagged == 0


class TestRecommendNmads:
    """Data-driven per-sample MAD recommendation (target-% heuristic)."""

    GRID = np.arange(2.0, 8.0001, 0.5)

    @staticmethod
    def _max_tail_pct(vals, nm):
        """Larger tail-removal percentage, which is the basis of the per-tail heuristic."""
        vals = np.asarray(vals, dtype=float)
        med = np.median(vals)
        mad = median_abs_deviation(vals)
        low = np.count_nonzero(vals < med - nm * mad) / vals.shape[0] * 100.0
        high = np.count_nonzero(vals > med + nm * mad) / vals.shape[0] * 100.0
        return max(low, high)

    def test_normal_sample_honors_per_tail_cap(self):
        """A clean, near-normal sample recommends the tightest per-tail-capped n."""
        rng = np.random.default_rng(0)
        vals = rng.normal(1000.0, 100.0, size=20000)
        rec = recommend_nmads(vals, cap_pct=5.0)
        assert rec is not None
        assert rec <= 5.0
        # Each tail at the recommendation must be at/under the cap.
        assert self._max_tail_pct(vals, rec) <= 5.0 + 1e-9
        # ...and it must be the *smallest* such n on the grid (tightest window).
        smaller = self.GRID[self.GRID < rec]
        assert all(self._max_tail_pct(vals, nm) > 5.0 for nm in smaller)

    def test_heavy_tailed_falls_back_to_max(self):
        """When even n=8 exceeds the cap on a tail, fall back to the most lenient value."""
        rng = np.random.default_rng(1)
        vals = rng.standard_cauchy(size=20000)  # very heavy tails
        # Sanity: even the widest window still removes more than a tiny cap per tail.
        assert self._max_tail_pct(vals, 8.0) > 0.5
        rec = recommend_nmads(vals, cap_pct=0.5)
        assert rec == 8.0

    def test_degenerate_returns_none(self):
        """All-equal values (MAD == 0) yield no recommendation."""
        assert recommend_nmads(np.full(500, 42.0)) is None
        assert recommend_nmads(np.array([])) is None

    def test_raising_cap_never_increases_recommendation(self):
        """A more permissive cap must never require a *larger* (more lenient) n."""
        rng = np.random.default_rng(2)
        vals = rng.normal(500.0, 80.0, size=20000)
        prev = None
        for cap in [1.0, 2.0, 5.0, 10.0, 15.0]:
            rec = recommend_nmads(vals, cap_pct=cap)
            assert rec is not None
            if prev is not None:
                assert rec <= prev
            prev = rec


class TestSeedPerSampleNmads:
    """The per-sample n_MADs auto-seeding helper in Data_Processing.py.

    That module runs Streamlit at import time, so we extract just this function's
    AST node and execute it against a stub ``st`` without a Streamlit runtime.
    Verifies: first entry seeds from recommendations; a degenerate rec (None or
    the NaN pandas coerces it to) falls back to 5.0; changing the cap re-seeds
    everything; and a manual override survives while the cap is unchanged.
    """

    @staticmethod
    def _load():
        import ast
        import types
        from pathlib import Path

        src_path = (
            Path(__file__).resolve().parents[1]
            / "src" / "cafe_app" / "tools" / "Data_Processing.py"
        )
        tree = ast.parse(src_path.read_text())
        stub_st = types.SimpleNamespace(session_state={})
        stub_st.session_state = dict(stub_st.session_state)  # a plain dict
        ns = {"pd": pd, "st": stub_st}
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "_seed_per_sample_nmads":
                exec(compile(ast.Module([node], []), "<extract>", "exec"), ns)
                return ns["_seed_per_sample_nmads"], stub_st
        raise AssertionError("_seed_per_sample_nmads not found")

    def test_first_entry_seeds_and_handles_degenerate(self):
        seed, stub_st = self._load()
        seed({"S1": 3.0, "S2": float("nan"), "S3": None, "S4": 6.5}, cap=5.0)
        assert stub_st.session_state["qc_nmads_S1"] == 3.0
        assert stub_st.session_state["qc_nmads_S2"] == 5.0  # NaN -> fallback
        assert stub_st.session_state["qc_nmads_S3"] == 5.0  # None -> fallback
        assert stub_st.session_state["qc_nmads_S4"] == 6.5

    def test_cap_change_reseeds_but_same_cap_preserves_override(self):
        seed, stub_st = self._load()
        seed({"S1": 3.0}, cap=5.0)
        # User overrides S1 in the fine-tune expander.
        stub_st.session_state["qc_nmads_S1"] = 7.0
        # Same cap, new recommendation -> override is preserved (not clobbered).
        seed({"S1": 3.5}, cap=5.0)
        assert stub_st.session_state["qc_nmads_S1"] == 7.0
        # Cap changes -> everything is re-seeded to the fresh recommendation.
        seed({"S1": 2.5}, cap=8.0)
        assert stub_st.session_state["qc_nmads_S1"] == 2.5


class TestReviewGate:
    """The step-advancement gate in Data_Processing.py.

    Guards the core fix: a finished step must NOT auto-advance because its results are
    stashed and shown, and only clicking Continue sets the step's done flag. The
    module runs Streamlit at import time, so we extract the helper's AST node and
    exec it against a stub ``st`` whose ``stop``/``rerun`` raise sentinels.
    """

    class _Stop(Exception):
        pass

    class _Rerun(Exception):
        pass

    def _load(self, button_returns):
        import ast
        import types
        from pathlib import Path

        stop_exc, rerun_exc = self._Stop, self._Rerun
        st = types.SimpleNamespace(session_state={}, rendered=[])
        st.button = lambda *a, **k: button_returns
        st.divider = lambda *a, **k: None
        st.stop = lambda: (_ for _ in ()).throw(stop_exc())
        st.rerun = lambda: (_ for _ in ()).throw(rerun_exc())
        for name in ["success", "info", "warning", "write", "caption",
                     "subheader", "image", "dataframe"]:
            st.__dict__[name] = (
                lambda n: (lambda p, **k: st.rendered.append((n, p)))
            )(name)

        src_path = (
            Path(__file__).resolve().parents[1]
            / "src" / "cafe_app" / "tools" / "Data_Processing.py"
        )
        tree = ast.parse(src_path.read_text())
        ns = {"st": st}
        wanted = {"_render_review_items", "_review_gate"}
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in wanted:
                exec(compile(ast.Module([node], []), "<x>", "exec"), ns)
        return ns["_review_gate"], st

    def test_no_pending_passes_through(self):
        gate, st = self._load(button_returns=False)
        gate("_qc_pending", "qc_done", "PCA")  # must not raise
        assert "qc_done" not in st.session_state

    def test_pending_without_click_renders_and_blocks(self):
        gate, st = self._load(button_returns=False)
        st.session_state["_qc_pending"] = [("success", "done"), ("write", "table")]
        with pytest.raises(self._Stop):
            gate("_qc_pending", "qc_done", "PCA")
        # Critical: it did NOT advance, and results stay pending.
        assert st.session_state.get("qc_done") is not True
        assert "_qc_pending" in st.session_state
        assert [k for k, _ in st.rendered] == ["success", "write"]

    def test_click_advances_and_clears_pending(self):
        gate, st = self._load(button_returns=True)
        st.session_state["_qc_pending"] = [("success", "done")]
        with pytest.raises(self._Rerun):
            gate("_qc_pending", "qc_done", "PCA")
        assert st.session_state["qc_done"] is True
        assert "_qc_pending" not in st.session_state


class TestPcaResultScatterItems:
    """PCA result plots must match the annotation columns that are actually present."""

    @staticmethod
    def _load():
        import ast
        from pathlib import Path

        src_path = (
            Path(__file__).resolve().parents[1]
            / "src" / "cafe_app" / "tools" / "Data_Processing.py"
        )
        tree = ast.parse(src_path.read_text())
        calls = []
        ns = {
            "GROUP_PALETTE": ["group-color"],
            "_pca_scatter": lambda adata, color, palette=None: calls.append(color) or color,
        }
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "_pca_result_scatter_items":
                exec(compile(ast.Module([node], []), "<extract>", "exec"), ns)
                return ns["_pca_result_scatter_items"], calls
        raise AssertionError("_pca_result_scatter_items not found")

    def test_missing_batch_ignores_stale_combat_plots(self, small_dense_adata):
        build_items, calls = self._load()

        items = build_items(
            small_dense_adata,
            before_batch=b"stale-batch",
            before_group=b"stale-group",
        )

        assert calls == ["Group"]
        assert items[0] == ("write", "PCA results by group:")
        assert all("ComBat" not in str(item) for item in items)

    def test_combat_comparison_requires_batch_and_group(self, small_dense_adata):
        build_items, calls = self._load()
        small_dense_adata.obs["Batch"] = "1"

        items = build_items(
            small_dense_adata,
            before_batch=b"before-batch",
            before_group=b"before-group",
        )

        assert calls == ["Batch", "Group"]
        assert any("colored by Batch" in str(item) for item in items)
        assert any("colored by Group" in str(item) for item in items)


class TestBuildOutputsZip:
    """The lazily-called output ZIP builder in Data_Processing.py.

    Guards the refactor that moved ZIP packaging out of the clustering step (so
    it can run under a spinner on the Export click), the ZIP must still contain
    exactly the same set of outputs. Extracted via AST with its two helper deps
    stubbed, so no Streamlit runtime is needed.
    """

    @staticmethod
    def _load():
        import ast
        import io
        import os
        import random
        import re
        import tempfile
        import zipfile
        from pathlib import Path

        src_path = (
            Path(__file__).resolve().parents[1]
            / "src" / "cafe_app" / "tools" / "Data_Processing.py"
        )
        tree = ast.parse(src_path.read_text())
        ns = {
            "io": io, "os": os, "random": random, "tempfile": tempfile,
            "zipfile": zipfile, "pd": pd, "np": np,
            "build_parameters_report": lambda a: "PARAMS REPORT",
            "_safe_filename": lambda s: (
                re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_") or "sample"
            ),
        }
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "build_outputs_zip":
                exec(compile(ast.Module([node], []), "<x>", "exec"), ns)
                return ns["build_outputs_zip"]
        raise AssertionError("build_outputs_zip not found")

    @staticmethod
    def _clustered_adata():
        rng = np.random.default_rng(0)
        n = 300
        X = rng.normal(size=(n, 4)).astype("float32")
        obs = pd.DataFrame(
            {
                "SampleID": np.repeat(["S1", "S2", "S3"], n // 3),
                "Group": np.repeat(["A", "B", "A"], n // 3),
                "leiden": pd.Categorical(rng.integers(0, 4, n).astype(str)),
            },
            index=[str(i) for i in range(n)],
        )
        adata = sc.AnnData(X, obs=obs)
        adata.var_names = [f"CD{i}" for i in range(4)]
        adata.obsm["X_umap"] = rng.normal(size=(n, 2))
        return adata

    def test_zip_contains_expected_outputs(self):
        import io
        import zipfile

        build_outputs_zip = self._load()
        zip_bytes, name = build_outputs_zip(self._clustered_adata(), resolution=0.5)
        assert name.startswith("analysis_outputs_") and name.endswith(".zip")

        names = zipfile.ZipFile(io.BytesIO(zip_bytes)).namelist()
        assert any(x.startswith("adata_0.5_") and x.endswith(".h5ad") for x in names)
        assert any(x.startswith("cluster_sample_counts_") for x in names)
        assert any(x.startswith("cluster_sample_frequencies_") for x in names)
        assert any(x.startswith("median_expression_") for x in names)
        assert any(x.startswith("analysis_parameters_") and x.endswith(".txt") for x in names)
        per_sample = sorted(x for x in names if x.startswith("per_sample_clustered/"))
        assert per_sample == [f"per_sample_clustered/{s}.csv" for s in ("S1", "S2", "S3")]

    def test_per_sample_csv_has_roundtrip_columns(self):
        import io
        import zipfile

        build_outputs_zip = self._load()
        zip_bytes, _ = build_outputs_zip(self._clustered_adata(), resolution=0.5)
        zf = zipfile.ZipFile(io.BytesIO(zip_bytes))
        sub = pd.read_csv(io.BytesIO(zf.read("per_sample_clustered/S1.csv")))
        assert {"Leiden_Cluster", "Group", "UMAP_1", "UMAP_2"}.issubset(sub.columns)

    def test_umap_pngs_are_bundled_when_provided(self):
        import io
        import zipfile

        build_outputs_zip = self._load()
        umap_pngs = {"umap_leiden_clusters": b"\x89PNG-leiden", "umap_Group": b"\x89PNG-group"}
        zip_bytes, _ = build_outputs_zip(
            self._clustered_adata(), resolution=0.5, umap_pngs=umap_pngs
        )
        names = zipfile.ZipFile(io.BytesIO(zip_bytes)).namelist()
        assert any(x.startswith("umap_plots/umap_leiden_clusters_") for x in names)
        assert any(x.startswith("umap_plots/umap_Group_") for x in names)

    def test_no_umap_plots_folder_without_pngs(self):
        import io
        import zipfile

        build_outputs_zip = self._load()
        zip_bytes, _ = build_outputs_zip(self._clustered_adata(), resolution=0.5)
        names = zipfile.ZipFile(io.BytesIO(zip_bytes)).namelist()
        assert not any(x.startswith("umap_plots/") for x in names)
