"""Guards run_umap_parallel's direct umap-learn call against dependency upgrades."""
import numpy as np
import scanpy as sc

from utils import run_umap_parallel


def test_run_umap_parallel_fills_x_umap():
    adata = sc.AnnData(np.random.default_rng(0).normal(size=(500, 8)))
    sc.pp.neighbors(adata, use_rep="X", random_state=50)
    run_umap_parallel(adata, min_dist=0.1)
    assert adata.obsm["X_umap"].shape == (500, 2)
    assert np.isfinite(adata.obsm["X_umap"]).all()
    assert {"a", "b"} <= adata.uns["umap"]["params"].keys()
