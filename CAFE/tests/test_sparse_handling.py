"""Tests for sparse/dense matrix handling (B7/B8).

Verifies that safe_flatten and safe_expression produce identical results
for sparse and dense representations of the same data, and that they
work correctly with AnnData slicing.
"""
import numpy as np
import pytest
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

from utils import safe_flatten, safe_expression


class TestAnnDataSparseDense:
    """Test that AnnData operations work identically for sparse and dense."""

    @pytest.fixture
    def dense_adata(self):
        np.random.seed(42)
        X = np.random.randn(50, 3) * 100 + 50
        obs = sc.AnnData(X).obs
        obs['leiden'] = ['0'] * 25 + ['1'] * 25
        adata = sc.AnnData(X)
        adata.obs['leiden'] = ['0'] * 25 + ['1'] * 25
        adata.var_names = ['CD1', 'CD2', 'CD3']
        return adata

    @pytest.fixture
    def sparse_adata(self, dense_adata):
        adata = dense_adata.copy()
        adata.X = csr_matrix(adata.X)
        return adata

    def test_flatten_same_result(self, dense_adata, sparse_adata):
        """safe_flatten must give same result for sparse and dense."""
        for marker in dense_adata.var_names:
            dense_flat = safe_flatten(dense_adata[:, marker].X)
            sparse_flat = safe_flatten(sparse_adata[:, marker].X)
            np.testing.assert_array_equal(dense_flat, sparse_flat)

    def test_expression_same_result(self, dense_adata, sparse_adata):
        """safe_expression must give same result for sparse and dense."""
        for marker in dense_adata.var_names:
            for method in ['Mean', 'Median']:
                dense_val = safe_expression(dense_adata[:, marker].X, method)
                sparse_val = safe_expression(sparse_adata[:, marker].X, method)
                assert dense_val == pytest.approx(sparse_val)

    def test_cluster_subset_expression(self, dense_adata, sparse_adata):
        """Expression for a cluster subset must match between sparse/dense."""
        for cluster in ['0', '1']:
            dense_subset = dense_adata[dense_adata.obs['leiden'] == cluster]
            sparse_subset = sparse_adata[sparse_adata.obs['leiden'] == cluster]

            for marker in dense_adata.var_names:
                dense_val = safe_expression(dense_subset[:, marker].X, 'Mean')
                sparse_val = safe_expression(sparse_subset[:, marker].X, 'Mean')
                assert dense_val == pytest.approx(sparse_val)

    def test_flatten_returns_ndarray(self, sparse_adata):
        """Critical: must return np.ndarray, not np.matrix."""
        for marker in sparse_adata.var_names:
            result = safe_flatten(sparse_adata[:, marker].X)
            assert isinstance(result, np.ndarray)
            assert not isinstance(result, np.matrix)

    def test_expression_returns_float(self, sparse_adata):
        """Must return a Python float, not np.matrix(1x1)."""
        for marker in sparse_adata.var_names:
            result = safe_expression(sparse_adata[:, marker].X, 'Mean')
            assert isinstance(result, float)

    def test_dataframe_from_sparse(self, sparse_adata):
        """Creating a DataFrame from sparse X must work (B8 fix pattern)."""
        markers = list(sparse_adata.var_names)
        X = sparse_adata[:, markers].X
        import pandas as pd
        df = pd.DataFrame(
            X.toarray() if issparse(X) else np.asarray(X),
            columns=markers
        )
        assert df.shape == (50, 3)
        assert list(df.columns) == markers
