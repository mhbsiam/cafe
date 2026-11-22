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


class TestCombatCovariates:
    """ComBat must pass Group as a covariate so an unbalanced design keeps the group shift."""

    @staticmethod
    def _make_adata(design, shift=2.0, seed=0):
        import pandas as pd
        import scanpy as sc
        rng = np.random.default_rng(seed)
        X, obs = [], []
        for i, (batch, group) in enumerate(design):
            x = rng.normal(0, 1, (200, 5)) + (shift if group == 'Aged' else 0) + (3.0 if batch == 'b2' else 0)
            X.append(x)
            obs += [(f'S{i}', group, batch)] * 200
        obs = pd.DataFrame(obs, columns=['SampleID', 'Group', 'Batch'])
        obs.index = obs.index.astype(str)
        return sc.AnnData(np.vstack(X), obs=obs)

    @staticmethod
    def _group_shift(adata):
        aged = (adata.obs['Group'] == 'Aged').values
        return adata.X[aged].mean() - adata.X[~aged].mean()

    def test_unbalanced_design_preserves_group_shift(self):
        from utils import run_combat_correction
        design = [('b1', 'Young')] * 3 + [('b1', 'Aged'), ('b2', 'Young')] + [('b2', 'Aged')] * 3
        adata = self._make_adata(design)
        obs_cols = list(adata.obs.columns)
        ok, msg = run_combat_correction(adata)
        assert ok, msg
        assert 'covariate' in msg
        assert self._group_shift(adata) == pytest.approx(2.0, rel=0.15)
        assert list(adata.obs.columns) == obs_cols

    def test_fully_confounded_design_refused(self):
        from utils import run_combat_correction
        design = [('b1', 'Young')] * 4 + [('b2', 'Aged')] * 4
        adata = self._make_adata(design)
        X_before = adata.X.copy()
        ok, _ = run_combat_correction(adata)
        assert not ok
        np.testing.assert_array_equal(adata.X, X_before)

    def test_single_level_group_runs_without_covariates(self):
        from utils import run_combat_correction
        adata = self._make_adata([('b1', 'Young')] * 2 + [('b2', 'Young')] * 2)
        obs_cols = list(adata.obs.columns)
        for covariates in (('Group',), ()):
            ok, msg = run_combat_correction(adata, covariates=covariates)
            assert ok, msg
            assert 'covariate' not in msg
            assert list(adata.obs.columns) == obs_cols

    def test_uploaded_loaders_read_named_buffers(self, monkeypatch):
        import io
        import utils
        monkeypatch.setattr(utils.st, 'session_state', {})

        def uploads():
            out = []
            for name, text in (('S1_Aged.csv', b"CD1,CD2\n1,2\n"), ('S2_Young.csv', b"CD2,CD1\n3,4\n")):
                f = io.BytesIO(text)
                f.name = name
                out.append(f)
            return out

        adata = utils.process_uploaded_csv_files(uploads())
        assert adata.n_obs == 2 and list(adata.var_names) == ['CD1', 'CD2']
        df = utils.process_uploaded_csv_to_df(uploads())
        assert df['CD1'].tolist() == [1.0, 4.0]

    def test_singular_design_refused_and_temp_columns_removed(self):
        from utils import run_combat_correction
        design = [('b1', 'Young')] * 3 + [('b1', 'Aged'), ('b2', 'Young')] + [('b2', 'Aged')] * 3
        adata = self._make_adata(design)
        adata.obs['Group2'] = adata.obs['Group']  # duplicate covariate -> singular design
        obs_cols = list(adata.obs.columns)
        X_before = adata.X.copy()
        ok, msg = run_combat_correction(adata, covariates=('Group', 'Group2'))
        assert not ok
        assert 'not applied' in msg
        np.testing.assert_array_equal(adata.X, X_before)
        assert list(adata.obs.columns) == obs_cols


class TestCsvLoading:
    """Shared CSV loader: empty columns, column order/type mismatches, missing rows, .CSV names, duplicates."""

    @staticmethod
    def _table(text):
        import io
        import pyarrow.csv as pv
        return pv.read_csv(io.BytesIO(text.encode()))

    def test_trailing_comma_keeps_rows_and_drops_empty_column(self):
        from utils import _combine_tables
        df, skipped, notices = _combine_tables([('S1_Aged.csv', self._table("CD1,CD2,\n1,2,\n3,4,\n"))])
        assert len(df) == 2
        assert list(df.columns) == ['CD1', 'CD2', 'SampleID', 'Group']
        assert not skipped
        assert any('empty column' in n for n in notices)

    def test_reordered_columns_load(self):
        from utils import _combine_tables
        df, skipped, _ = _combine_tables([
            ('S1_Aged.csv', self._table("CD1,CD2\n1,2\n")),
            ('S2_Young.csv', self._table("CD2,CD1\n5.5,6\n")),
        ])
        assert not skipped
        assert df['CD1'].tolist() == [1.0, 6.0]
        assert df['CD2'].tolist() == [2.0, 5.5]
        assert df['CD1'].dtype == np.float64

    def test_string_marker_file_is_reported_not_raised(self):
        from utils import _combine_tables
        df, skipped, _ = _combine_tables([
            ('S1_Aged.csv', self._table("CD1,CD2\n1,2\n")),
            ('S2_Young.csv', self._table("CD1,CD2\n1,oops\n")),
        ])
        assert df['SampleID'].tolist() == ['S1']
        assert len(skipped) == 1 and 'S2_Young.csv' in skipped[0]

    def test_text_column_in_every_file_is_dropped_with_notice(self, monkeypatch):
        import io
        import utils
        monkeypatch.setattr(utils.st, 'session_state', {})
        uploads = []
        for name, text in (('S1_Aged.csv', b"CD1,Well,CD2\n1,A1,2\n"), ('S2_Young.csv', b"Well,CD2,CD1\nB2,3,4\n")):
            f = io.BytesIO(text)
            f.name = name
            uploads.append(f)
        adata = utils.process_uploaded_csv_files(uploads)
        assert adata.n_obs == 2
        assert list(adata.var_names) == ['CD1', 'CD2']
        assert 'Well' in utils.st.session_state['_load_notice']

    def test_skip_messages_reach_load_notice(self, monkeypatch):
        import utils
        monkeypatch.setattr(utils.st, 'session_state', {'_load_notice': 'earlier'})
        utils._report_load(['Skipping file X.csv'], [], save_notice=True)
        assert utils.st.session_state['_load_notice'] == 'Skipping file X.csv'
        utils._report_load([], [], save_notice=True)
        assert '_load_notice' not in utils.st.session_state

    def test_repeated_load_does_not_grow_notice(self, monkeypatch):
        import io
        import utils
        monkeypatch.setattr(utils.st, 'session_state', {})

        def load():
            f = io.BytesIO(b"CD1,CD2\n1,2\n3,\n")
            f.name = 'S1_Aged.csv'
            utils.process_uploaded_csv_files([f])

        load()
        once = utils.st.session_state['_load_notice']
        load()
        assert utils.st.session_state['_load_notice'] == once

    def test_other_loaders_leave_data_processing_notice_alone(self, tmp_path, monkeypatch):
        import io
        import utils
        monkeypatch.setattr(utils.st, 'session_state', {'_load_notice': 'unseen'})
        (tmp_path / 'S1_Aged.csv').write_text("CD1,CD2\n1,2\n")
        assert utils.process_directory_csv_files(str(tmp_path)) is not None
        f = io.BytesIO(b"CD1,CD2\n1,2\n3,\n")
        f.name = 'S1_Aged.csv'
        assert utils.process_uploaded_csv_to_df([f]) is not None  # warns, but must not overwrite either
        assert utils.st.session_state['_load_notice'] == 'unseen'

    def test_mismatched_columns_file_is_reported(self):
        from utils import _combine_tables
        df, skipped, _ = _combine_tables([
            ('S1_Aged.csv', self._table("CD1,CD2\n1,2\n")),
            ('S2_Young.csv', self._table("CD1,CD3\n1,2\n")),
        ])
        assert df['SampleID'].tolist() == ['S1']
        assert len(skipped) == 1 and 'S2_Young.csv' in skipped[0]

    def test_all_rows_missing_returns_none(self):
        from utils import _combine_tables
        df, _, notices = _combine_tables([('S1_Aged.csv', self._table("CD1,CD2\n1,\n,2\n"))])
        assert df is None
        assert any('No rows left' in n for n in notices)

    def test_directory_loader_uppercase_csv_keeps_load_notice(self, tmp_path, monkeypatch):
        import utils
        monkeypatch.setattr(utils.st, 'session_state', {'_load_notice': 'earlier'})
        (tmp_path / 'S1_Aged.CSV').write_text("CD1,CD2\n1,2\n3,\n")
        adata = utils.process_directory_csv_files(str(tmp_path))
        assert adata.n_obs == 1
        assert adata.obs['SampleID'].tolist() == ['S1']
        assert adata.obs['Group'].tolist() == ['Aged']
        assert utils.st.session_state['_load_notice'] == 'earlier'  # only Data Processing's loader owns the key

    def test_inspect_flags_duplicate_sample_ids_and_file_names(self):
        import io
        from utils import inspect_uploaded_csv_files

        def upload(name):
            f = io.BytesIO(b"CD1,CD2\n1,2\n")
            f.name = name
            return f

        review = inspect_uploaded_csv_files([upload('S1_Aged.csv'), upload('S1_Young.csv'), upload('S2_Aged.CSV'),
                                             upload('S2_Aged.CSV')])
        assert len(review['duplicates']) == 2
        assert all(d in review['bad_names'] for d in review['duplicates'])
        assert review['files'][2]['Group'] == 'Aged'
        ok = inspect_uploaded_csv_files([upload('S1_Aged.csv'), upload('S2_Young.csv')])
        assert not ok['bad_names'] and not ok['duplicates']
