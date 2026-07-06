"""Tests for the per-sample differential abundance test (tab34 fix).

The old tab34 ran chi2_contingency on a per-*cell* contingency table, which is
pseudoreplicated: with many cells the p-value collapses to ~0 for any dataset.
These tests exercise the replacement (compute_sample_proportions +
differential_abundance_test in utils.py), which uses the sample as the unit of
analysis.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import pytest

from utils import compute_sample_proportions, differential_abundance_test


def _build_adata(cells_per_sample_cluster):
    """Build an AnnData from a {(sample, group, cluster): n_cells} dict."""
    rows = []
    for (sample, group, cluster), n in cells_per_sample_cluster.items():
        rows.extend([(sample, group, cluster)] * n)
    obs = pd.DataFrame(rows, columns=['SampleID', 'Group', 'leiden'])
    obs.index = obs.index.astype(str)
    X = np.zeros((len(obs), 2))
    adata = sc.AnnData(X, obs=obs)
    adata.var_names = ['m0', 'm1']
    return adata


def _da_dataset(scale=1):
    """5 samples/group, 2 groups, 3 clusters (100 cells/sample before scaling).

    Cluster 'A' is strongly enriched in Disease (10% vs 50%), cluster 'B' is its
    complement (enriched in Healthy), and cluster 'C' is genuinely balanced —
    its per-sample proportions overlap between groups, so it must NOT come out
    significant.  (With only two clusters a "balanced" one is impossible, since
    proportions sum to 1.)  ``scale`` multiplies all cell counts without
    changing proportions — used for the pseudoreplication regression check.
    """
    c_healthy = [30, 28, 32, 29, 31]
    c_disease = [31, 29, 33, 28, 30]
    data = {}
    for i in range(5):
        s = f"H{i}"
        c = c_healthy[i]
        data[(s, 'Healthy', 'A')] = 10 * scale
        data[(s, 'Healthy', 'C')] = c * scale
        data[(s, 'Healthy', 'B')] = (100 - 10 - c) * scale
    for i in range(5):
        s = f"D{i}"
        c = c_disease[i]
        data[(s, 'Disease', 'A')] = 50 * scale
        data[(s, 'Disease', 'C')] = c * scale
        data[(s, 'Disease', 'B')] = (100 - 50 - c) * scale
    return _build_adata(data)


class TestComputeSampleProportions:
    def test_rows_sum_to_one(self):
        adata = _da_dataset()
        prop_df, sample_group = compute_sample_proportions(adata, 'leiden')
        assert np.allclose(prop_df.sum(axis=1), 1.0)
        assert set(sample_group.unique()) == {'Healthy', 'Disease'}
        assert prop_df.shape == (10, 3)

    def test_zero_cell_samples_get_zero_not_dropped(self):
        # S2 has no cluster '1' cells at all -> proportion 0, still present.
        adata = _build_adata({
            ('S1', 'G', '0'): 10, ('S1', 'G', '1'): 10,
            ('S2', 'G', '0'): 20,  # no cluster '1'
        })
        prop_df, _ = compute_sample_proportions(adata, 'leiden')
        assert prop_df.loc['S2', '1'] == 0.0
        assert prop_df.loc['S2', '0'] == 1.0


class TestDifferentialAbundanceTest:
    def test_one_row_per_cluster_valid_pvalues(self):
        adata = _da_dataset()
        prop_df, sample_group = compute_sample_proportions(adata, 'leiden')
        res = differential_abundance_test(prop_df, sample_group)
        assert set(res['cluster']) == {'A', 'B', 'C'}
        pvals = res['p_value'].dropna()
        assert ((pvals >= 0) & (pvals <= 1)).all()
        assert ((res['p_adj'].dropna() >= 0) & (res['p_adj'].dropna() <= 1)).all()

    def test_associated_cluster_significant_balanced_not(self):
        adata = _da_dataset()
        prop_df, sample_group = compute_sample_proportions(adata, 'leiden')
        res = differential_abundance_test(prop_df, sample_group).set_index('cluster')
        assert res.loc['A', 'p_adj'] < 0.05
        assert bool(res.loc['A', 'significant']) is True
        assert 'enriched in Disease' == res.loc['A', 'direction']
        # Balanced cluster C must NOT be flagged (proves it no longer flags everything).
        assert not bool(res.loc['C', 'significant'])

    def test_pvalues_invariant_to_cell_count_scaling(self):
        """The core pseudoreplication regression check.

        Scaling all cell counts 1000x leaves proportions unchanged, so a
        sample-level test must return identical p-values — unlike the old
        per-cell chi-square whose p-value would collapse toward 0.
        """
        small = _da_dataset(scale=1)
        big = _da_dataset(scale=1000)

        res_small = differential_abundance_test(
            *compute_sample_proportions(small, 'leiden')).set_index('cluster')
        res_big = differential_abundance_test(
            *compute_sample_proportions(big, 'leiden')).set_index('cluster')

        for cluster in ['A', 'B', 'C']:
            assert np.isclose(res_small.loc[cluster, 'p_value'],
                              res_big.loc[cluster, 'p_value'])

    def test_group_with_single_sample_yields_nan_not_crash(self):
        adata = _build_adata({
            ('H0', 'Healthy', 'A'): 10, ('H0', 'Healthy', 'B'): 90,
            ('H1', 'Healthy', 'A'): 12, ('H1', 'Healthy', 'B'): 88,
            ('D0', 'Disease', 'A'): 50, ('D0', 'Disease', 'B'): 50,  # only 1 Disease sample
        })
        prop_df, sample_group = compute_sample_proportions(adata, 'leiden')
        res = differential_abundance_test(prop_df, sample_group)
        # Disease has <2 samples -> no cluster testable -> all p-values NaN, no crash.
        assert res['p_value'].isna().all()
        assert not res['significant'].any()
