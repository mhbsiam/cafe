"""Tests for CAFE shared utilities (utils.py).

Covers: safe_flatten (B8), safe_expression (B7), safe_sort_clusters (B9),
get_cluster_col, validate_path (S1).
"""
import os
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from utils import (
    safe_flatten,
    safe_expression,
    safe_sort_clusters,
    get_cluster_col,
    validate_path,
    RANDOM_STATE,
    MARKER_HIGH_THRESHOLD,
    MARKER_LOW_THRESHOLD,
)


class TestSafeFlatten:
    """B8: safe_flatten must work on both sparse and dense matrices."""

    def test_dense_1d(self):
        X = np.array([1.0, 2.0, 3.0])
        result = safe_flatten(X)
        assert isinstance(result, np.ndarray)
        assert result.shape == (3,)
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])

    def test_dense_2d_single_column(self):
        X = np.array([[1.0], [2.0], [3.0]])
        result = safe_flatten(X)
        assert isinstance(result, np.ndarray)
        assert result.shape == (3,)

    def test_sparse_1d(self):
        X = csr_matrix(np.array([[1.0], [2.0], [3.0]]))
        result = safe_flatten(X)
        assert isinstance(result, np.ndarray)
        assert result.shape == (3,)
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])

    def test_sparse_with_explicit_zeros(self):
        """Sparse matrices with explicit zeros must densify correctly."""
        dense = np.array([[0.0], [5.0], [0.0]])
        sparse = csr_matrix(dense)
        result = safe_flatten(sparse)
        np.testing.assert_array_equal(result, [0.0, 5.0, 0.0])

    def test_returns_ndarray_not_matrix(self):
        """Critical: .X.flatten() on sparse returns np.matrix, not ndarray."""
        X = csr_matrix(np.array([[1.0], [2.0]]))
        result = safe_flatten(X)
        assert not isinstance(result, np.matrix)


class TestSafeExpression:
    """B7: safe_expression must return scalars for both sparse and dense."""

    def test_mean_dense(self):
        X = np.array([1.0, 2.0, 3.0])
        result = safe_expression(X, 'Mean')
        assert isinstance(result, float)
        assert result == pytest.approx(2.0)

    def test_median_dense(self):
        X = np.array([1.0, 2.0, 3.0, 4.0])
        result = safe_expression(X, 'Median')
        assert isinstance(result, float)
        assert result == pytest.approx(2.5)

    def test_mean_sparse(self):
        """B7: data.X.mean() on sparse returns np.matrix (1x1), not scalar."""
        X = csr_matrix(np.array([[1.0], [2.0], [3.0]]))
        result = safe_expression(X, 'Mean')
        assert isinstance(result, float)
        assert result == pytest.approx(2.0)

    def test_median_sparse(self):
        """B7: np.median on sparse operates on raw sparse data (wrong)."""
        X = csr_matrix(np.array([[1.0], [2.0], [3.0], [4.0]]))
        result = safe_expression(X, 'Median')
        assert isinstance(result, float)
        assert result == pytest.approx(2.5)

    def test_sparse_and_dense_match(self):
        """Same data in sparse and dense form must give same result."""
        dense = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])
        sparse = csr_matrix(dense)
        assert safe_expression(dense, 'Mean') == pytest.approx(safe_expression(sparse, 'Mean'))
        assert safe_expression(dense, 'Median') == pytest.approx(safe_expression(sparse, 'Median'))


class TestSafeSortClusters:
    """B9: safe_sort_clusters must handle non-integer cluster labels."""

    def test_integer_labels(self):
        labels = ['2', '0', '1', '3']
        result = safe_sort_clusters(labels)
        assert result == ['0', '1', '2', '3']

    def test_string_labels(self):
        labels = ['T_cells', 'B_cells', 'NK_cells']
        result = safe_sort_clusters(labels)
        assert result == ['B_cells', 'NK_cells', 'T_cells']

    def test_mixed_labels(self):
        labels = ['10', '2', '1']
        result = safe_sort_clusters(labels)
        assert result == ['1', '2', '10']

    def test_non_integer_does_not_crash(self):
        """B9: sorted(..., key=int) crashes on non-integer labels."""
        labels = ['T_cells', 'B_cells']
        # Should not raise
        result = safe_sort_clusters(labels)
        assert len(result) == 2


class TestGetClusterCol:
    def test_returns_cell_type_when_present(self, adata_with_cell_type):
        assert get_cluster_col(adata_with_cell_type) == 'cell_type'

    def test_returns_leiden_when_no_cell_type(self, small_dense_adata):
        assert get_cluster_col(small_dense_adata) == 'leiden'


class TestValidatePath:
    """S1: validate_path must reject paths outside the allowed base."""

    def test_valid_path(self, tmp_path):
        result = validate_path(str(tmp_path), base=str(tmp_path))
        assert result == str(tmp_path)

    def test_rejects_traversal(self):
        with pytest.raises(ValueError, match="outside the allowed directory"):
            validate_path("/etc/passwd", base="/tmp/safe_dir")

    def test_rejects_dotdot(self, tmp_path):
        base = str(tmp_path / "safe")
        os.makedirs(base, exist_ok=True)
        with pytest.raises(ValueError):
            validate_path(f"{base}/../../../etc/passwd", base=base)


class TestConstants:
    def test_random_state(self):
        assert RANDOM_STATE == 50

    def test_marker_thresholds(self):
        assert MARKER_HIGH_THRESHOLD == 1000
        assert MARKER_LOW_THRESHOLD == 0
