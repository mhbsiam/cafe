"""Leiden stability helpers: ARI sweep, resolution suggestion, per-cluster scores."""
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from utils import _leiden_labels, _neighbor_graph, cluster_stability, leiden_stability, suggest_resolution


@pytest.fixture(scope="module")
def blobs():
    rng = np.random.default_rng(0)
    centers = rng.normal(0, 6, (4, 8))
    adata = sc.AnnData(centers[rng.integers(0, 4, 1200)] + rng.normal(0, 1, (1200, 8)))
    sc.pp.neighbors(adata, n_neighbors=15, method="umap", use_rep="X", random_state=50)
    return adata


@pytest.mark.parametrize("resolution", [0.2, 1.0])
def test_direct_leiden_matches_scanpy(blobs, resolution):
    sc.tl.leiden(blobs, resolution=resolution, random_state=50, flavor="igraph", n_iterations=2, directed=False)
    direct = _leiden_labels(_neighbor_graph(blobs), resolution, 50)
    assert (blobs.obs["leiden"].astype(int).to_numpy() == direct).all()


def test_leiden_stability_columns_and_progress(blobs):
    calls = []
    table, n = leiden_stability(blobs, [0.1, 0.5], seeds=3, subsamples=2, progress=lambda d, t: calls.append((d, t)))
    assert n == 1200 and list(table["resolution"]) == [0.1, 0.5]
    assert {"clusters", "seed_ari", "seed_sd", "sub_ari", "sub_sd"} <= set(table.columns)
    assert table["sub_ari"].iloc[0] > 0.95  # four well-separated blobs at low resolution
    assert calls[-1] == (calls[-1][1], calls[-1][1]) and len(calls) == calls[-1][1]


def test_leiden_stability_samples_above_max_cells(blobs):
    table, n = leiden_stability(blobs, [0.1], seeds=0, subsamples=2, max_cells=500)
    assert n == 500 and "seed_ari" not in table and table["sub_ari"].between(0, 1).all()


# Measured on a 100,000-cell sample of the FR-FCM-Z3WR dataset, resolutions 0.1 to 1.0, with four different
# random states (reference seeds and subsample draws). A single reference run had flipped between 0.1 and 0.2.
_REAL = [
    ([.73, .68, .56, .53, .55, .53, .52, .52, .53, .52], [.19, .13, .06, .05, .06, .06, .04, .04, .04, .04]),
    ([.76, .75, .62, .55, .55, .55, .55, .56, .53, .53], [.19, .11, .10, .07, .07, .07, .07, .06, .07, .06]),
    ([.82, .78, .59, .59, .57, .55, .54, .53, .53, .53], [.21, .14, .08, .07, .05, .06, .06, .05, .05, .05]),
    ([.76, .73, .62, .57, .54, .56, .56, .54, .54, .53], [.21, .16, .11, .06, .06, .06, .06, .05, .05, .04]),
]


@pytest.mark.parametrize("sub, sub_sd", _REAL)
def test_suggestion_is_the_same_for_every_random_state(sub, sub_sd):
    table = pd.DataFrame({"resolution": [round(0.1 * i, 1) for i in range(1, 11)], "sub_ari": sub, "sub_sd": sub_sd})
    suggested, lowest, best = suggest_resolution(table)
    assert (suggested, lowest) == (0.2, 0.1) and best == max(sub)


def test_suggestion_averages_both_measures():
    table = pd.DataFrame({"resolution": [0.1, 0.2, 0.3], "sub_ari": [0.9, 0.5, 0.5], "sub_sd": [0.02] * 3,
                          "seed_ari": [0.5, 0.9, 0.5], "seed_sd": [0.02] * 3})
    assert suggest_resolution(table, seeds=10) == (0.2, 0.1, 0.7)


def test_suggestion_ignores_a_far_resolution_that_ties():
    table = pd.DataFrame({"resolution": [0.1, 0.2, 0.3, 0.9], "sub_ari": [0.9, 0.9, 0.5, 0.9], "sub_sd": [0.05] * 4})
    assert suggest_resolution(table)[0] == 0.2


def test_cluster_stability_scores_every_cluster(blobs):
    sc.tl.leiden(blobs, resolution=0.2, random_state=50, flavor="igraph", n_iterations=2, directed=False)
    out = cluster_stability(blobs, 0.2, n_seeds=3)
    assert set(out["cluster"]) == set(blobs.obs["leiden"].astype(str))
    assert out["cells"].sum() == blobs.n_obs and out["stability"].between(0, 1).all()
    assert out["stability"].max() > 0.95
