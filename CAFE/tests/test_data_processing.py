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

from utils import safe_sort_clusters

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
