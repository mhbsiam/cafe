"""Pytest fixtures for CAFE tests.

Provides small AnnData objects (both sparse and dense) and helper utilities
for testing the bug fixes from the code review.
"""
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import pytest

# Make the Tools directory importable
TOOLS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'cafe_app', 'tools')
sys.path.insert(0, TOOLS_DIR)


@pytest.fixture
def small_dense_adata():
    """A small dense AnnData with leiden clusters, SampleID, and Group."""
    np.random.seed(42)
    n_cells = 100
    n_markers = 5
    marker_names = [f"CD{i}" for i in range(n_markers)]

    X = np.random.randn(n_cells, n_markers) * 100 + 50
    obs = pd.DataFrame({
        'leiden': ['0'] * 50 + ['1'] * 30 + ['2'] * 20,
        'SampleID': ['S1'] * 50 + ['S2'] * 30 + ['S3'] * 20,
        'Group': ['Healthy'] * 50 + ['Disease'] * 30 + ['Healthy'] * 20,
    }, index=[str(i) for i in range(n_cells)])

    adata = sc.AnnData(X, obs=obs)
    adata.var_names = marker_names
    return adata


@pytest.fixture
def small_sparse_adata(small_dense_adata):
    """Same data as small_dense_adata but with a sparse X matrix."""
    adata = small_dense_adata.copy()
    adata.X = csr_matrix(adata.X)
    return adata


@pytest.fixture
def adata_with_cell_type(small_dense_adata):
    """AnnData with a 'cell_type' column (for testing get_cluster_col)."""
    adata = small_dense_adata.copy()
    adata.obs['cell_type'] = ['T cell'] * 50 + ['B cell'] * 30 + ['NK cell'] * 20
    return adata


@pytest.fixture
def adata_non_integer_clusters(small_dense_adata):
    """AnnData with non-integer cluster labels (for testing safe_sort_clusters)."""
    adata = small_dense_adata.copy()
    adata.obs['leiden'] = ['T_cells'] * 50 + ['B_cells'] * 30 + ['NK_cells'] * 20
    return adata
