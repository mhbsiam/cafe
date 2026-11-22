"""Shared-module fixes: embedding explorer background statistics cache."""
import numpy as np
import anndata as ad


def _adata(X):
    a = ad.AnnData(X.astype("float32"))
    a.var_names = [f"CD{i}" for i in range(X.shape[1])]
    return a


def test_background_stats_recomputed_for_same_shape_different_values(monkeypatch):
    import embedding_explorer

    monkeypatch.setattr(embedding_explorer.st, "session_state", {})
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 3))
    mean1, median1 = embedding_explorer._background_stats(_adata(X))
    mean2, median2 = embedding_explorer._background_stats(_adata(X + 5.0))
    assert not np.allclose(mean1, mean2)
    assert not np.allclose(median1, median2)
    np.testing.assert_allclose(mean2, (X + 5.0).mean(axis=0), rtol=1e-5)
    # Same data again hits the cache and returns the same numbers.
    mean3, _ = embedding_explorer._background_stats(_adata(X + 5.0))
    np.testing.assert_allclose(mean3, mean2)
