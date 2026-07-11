"""Tests for the probabilistic GMM cell-type annotation core (celltype_gmm.py)."""
import os
import sys

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy.sparse import csr_matrix

TOOLS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'cafe_app', 'tools')
sys.path.insert(0, TOOLS_DIR)

from celltype_gmm import (  # noqa: E402
    fit_marker_gmm,
    compute_positive_posteriors,
    score_cell_types,
    assign_cell_types,
    suggest_thresholds,
)


def _bimodal(neg_mean, pos_mean, n_neg, n_pos, sd=1.0, seed=0):
    """Two well-separated Gaussian humps concatenated (neg then pos)."""
    rng = np.random.default_rng(seed)
    neg = rng.normal(neg_mean, sd, n_neg)
    pos = rng.normal(pos_mean, sd, n_pos)
    return np.concatenate([neg, pos])


@pytest.fixture
def bimodal_adata():
    """AnnData with two samples; each has clearly bimodal CD8 and CD4 markers.

    Sample S1's positive population sits at a *higher* raw intensity than S2's
    (intensity drift) so per-sample fitting is exercised.  Layout per sample:
    first half are CD8+ CD4-, second half are CD8- CD4+.
    """
    rng = np.random.default_rng(1)
    per_sample = 400
    half = per_sample // 2

    def sample_block(cd8_neg, cd8_pos, cd4_neg, cd4_pos, seed):
        # first `half` cells: CD8 high, CD4 low; next `half`: CD8 low, CD4 high
        cd8 = np.concatenate([
            rng.normal(cd8_pos, 1.0, half),
            rng.normal(cd8_neg, 1.0, half),
        ])
        cd4 = np.concatenate([
            rng.normal(cd4_neg, 1.0, half),
            rng.normal(cd4_pos, 1.0, half),
        ])
        return np.column_stack([cd8, cd4])

    s1 = sample_block(cd8_neg=2, cd8_pos=12, cd4_neg=2, cd4_pos=12, seed=1)
    s2 = sample_block(cd8_neg=5, cd8_pos=18, cd4_neg=5, cd4_pos=18, seed=2)
    X = np.vstack([s1, s2])

    obs = pd.DataFrame({
        "SampleID": ["S1"] * per_sample + ["S2"] * per_sample,
        "Group": ["A"] * per_sample + ["B"] * per_sample,
        "leiden": (["0"] * half + ["1"] * half) * 2,
    }, index=[str(i) for i in range(2 * per_sample)])

    adata = sc.AnnData(X, obs=obs)
    adata.var_names = ["CD8", "CD4"]
    return adata


# ── fit_marker_gmm ────────────────────────────────────────────────────────────

def test_fit_marker_gmm_separates_two_humps():
    values = _bimodal(neg_mean=0, pos_mean=10, n_neg=300, n_pos=300)
    proba, info = fit_marker_gmm(values, random_state=0)

    assert info["informative"] is True
    # low hump -> low posterior, high hump -> high posterior
    assert proba[:300].mean() < 0.1
    assert proba[300:].mean() > 0.9
    # cutoff lands between the two means
    assert 0 < info["cutoff"] < 10


def test_fit_marker_gmm_unimodal_is_uninformative():
    rng = np.random.default_rng(0)
    values = rng.normal(5.0, 1.0, 500)  # single hump
    proba, info = fit_marker_gmm(values, random_state=0)
    assert info["informative"] is False


def test_fit_marker_gmm_zero_variance_returns_half():
    values = np.full(200, 3.0)
    proba, info = fit_marker_gmm(values)
    assert np.allclose(proba, 0.5)
    assert info["informative"] is False


def test_fit_marker_gmm_too_few_cells_returns_half():
    proba, info = fit_marker_gmm(np.array([1.0, 2.0, 3.0]))
    assert np.allclose(proba, 0.5)
    assert info["informative"] is False


def test_fit_marker_gmm_accepts_sparse_column():
    values = _bimodal(0, 10, 200, 200)
    sparse_col = csr_matrix(values.reshape(-1, 1))
    proba, info = fit_marker_gmm(sparse_col, random_state=0)
    assert proba.shape[0] == 400
    assert info["informative"] is True


# ── compute_positive_posteriors (per-sample independence) ─────────────────────

def test_per_sample_posteriors_gate_each_sample(bimodal_adata):
    posteriors, diag = compute_positive_posteriors(bimodal_adata, random_state=0)

    assert posteriors.shape == (bimodal_adata.n_obs, 2)
    assert list(posteriors.columns) == ["CD8", "CD4"]
    # diagnostics: one row per (sample, marker) = 2 x 2
    assert len(diag) == 4
    assert diag["informative"].all()

    obs = bimodal_adata.obs
    # In both samples the CD8-high block (first half of each sample) should read high.
    for sample in ["S1", "S2"]:
        block = obs["SampleID"] == sample
        cd8 = posteriors.loc[block.values, "CD8"].to_numpy()
        n = cd8.shape[0]
        assert cd8[: n // 2].mean() > 0.9   # CD8+ block
        assert cd8[n // 2:].mean() < 0.1    # CD8- block


def test_missing_sampleid_falls_back_to_global(bimodal_adata):
    ad = bimodal_adata.copy()
    del ad.obs["SampleID"]
    posteriors, diag = compute_positive_posteriors(ad, random_state=0)
    assert posteriors.shape[0] == ad.n_obs
    # one global fit per marker
    assert set(diag["SampleID"]) == {"__all__"}


def test_pooled_fallback_rescues_minority_sample():
    """A marker that's bimodal overall but nearly all-negative in one sample
    should still be graded in that sample via the pooled gate (not 0.5)."""
    rng = np.random.default_rng(3)
    # S1: balanced neg/pos. S2: almost all negative (only ~3% positive), which is too
    # few positives to fit a per-sample bimodal, but pooled is clearly bimodal.
    s1 = _bimodal(neg_mean=0, pos_mean=10, n_neg=250, n_pos=250, seed=3)
    s2 = np.concatenate([
        rng.normal(0, 1.0, 485),
        rng.normal(10, 1.0, 15),
    ])
    X = np.concatenate([s1, s2]).reshape(-1, 1)
    obs = pd.DataFrame(
        {"SampleID": ["S1"] * 500 + ["S2"] * 500},
        index=[str(i) for i in range(1000)],
    )
    ad = sc.AnnData(X, obs=obs)
    ad.var_names = ["CD8"]

    posteriors, diag = compute_positive_posteriors(ad, random_state=0)

    s2_diag = diag[diag["SampleID"] == "S2"].iloc[0]
    # per-sample fit is not usable, so it must not be neutralised
    assert s2_diag["source"] in ("per-sample", "pooled")
    assert s2_diag["source"] != "neutral"

    # the 15 genuine S2 positives (last 15 rows) should read clearly positive
    s2_pos = posteriors["CD8"].to_numpy()[-15:]
    assert s2_pos.mean() > 0.8


def test_uninformative_marker_is_neutralised():
    rng = np.random.default_rng(0)
    n = 300
    X = np.column_stack([
        _bimodal(0, 10, n // 2, n // 2),   # informative
        rng.normal(5.0, 1.0, n),           # unimodal -> neutralised
    ])
    obs = pd.DataFrame({"SampleID": ["S1"] * n}, index=[str(i) for i in range(n)])
    ad = sc.AnnData(X, obs=obs)
    ad.var_names = ["good", "flat"]

    posteriors, diag = compute_positive_posteriors(ad, random_state=0)
    assert np.allclose(posteriors["flat"].to_numpy(), 0.5)


# ── score_cell_types ──────────────────────────────────────────────────────────

def test_score_cell_types_geometric_mean():
    posteriors = pd.DataFrame({"CD8": [0.9, 0.1], "CD3": [0.4, 0.8]})
    definitions = {"CD8 T": {"CD8": "+", "CD3": "+"}}
    scores = score_cell_types(posteriors, definitions)
    expected = np.sqrt(np.array([0.9 * 0.4, 0.1 * 0.8]))
    assert np.allclose(scores["CD8 T"].to_numpy(), expected)


def test_score_cell_types_negative_marker_inverts():
    posteriors = pd.DataFrame({"CD8": [0.9], "CD4": [0.2]})
    # CD8 T = CD8+ CD4-  -> geomean(0.9, 1-0.2)
    scores = score_cell_types(posteriors, {"CD8 T": {"CD8": "+", "CD4": "-"}})
    assert np.isclose(scores["CD8 T"].iloc[0], np.sqrt(0.9 * 0.8))


def test_score_cell_types_skips_types_without_markers():
    posteriors = pd.DataFrame({"CD8": [0.9]})
    scores = score_cell_types(posteriors, {"Junk": {"CD8": ""}, "Missing": {"XYZ": "+"}})
    assert scores.shape[1] == 0


# ── assign_cell_types ─────────────────────────────────────────────────────────

def test_assign_high_confidence():
    scores = pd.DataFrame({"CD8 T": [0.95], "CD4 T": [0.05]})
    calls = assign_cell_types(scores, min_score=0.5, high_margin=0.5, low_margin=0.15)
    assert calls["cell_type"].iloc[0] == "CD8 T"
    assert calls["confidence_level"].iloc[0] == "High"


def test_assign_ambiguous_tie():
    # two types score nearly equally -> small margin -> Ambiguous
    scores = pd.DataFrame({"CD8 T": [0.8], "CD4 T": [0.78]})
    calls = assign_cell_types(scores, min_score=0.5, high_margin=0.5, low_margin=0.15)
    assert calls["confidence_level"].iloc[0] == "Ambiguous"


def test_assign_low_confidence_middle_margin():
    scores = pd.DataFrame({"CD8 T": [0.7], "CD4 T": [0.3]})
    calls = assign_cell_types(scores, min_score=0.5, high_margin=0.5, low_margin=0.15)
    # margin = 0.7 - 0.3 = 0.4  -> between low(0.15) and high(0.5)
    assert calls["confidence_level"].iloc[0] == "Low"
    assert calls["cell_type"].iloc[0] == "CD8 T"


def test_assign_unassigned_below_min_score():
    scores = pd.DataFrame({"CD8 T": [0.2], "CD4 T": [0.1]})
    calls = assign_cell_types(scores, min_score=0.5, high_margin=0.5, low_margin=0.15)
    assert calls["cell_type"].iloc[0] == "Unassigned"
    assert calls["confidence_level"].iloc[0] == "Ambiguous"


def test_assign_empty_scores_all_unassigned():
    scores = pd.DataFrame(index=[0, 1, 2])
    calls = assign_cell_types(scores)
    assert (calls["cell_type"] == "Unassigned").all()
    assert (calls["confidence_level"] == "Ambiguous").all()


def test_assign_single_type_high_when_matched():
    scores = pd.DataFrame({"CD8 T": [0.9, 0.2]})
    calls = assign_cell_types(scores, min_score=0.5, high_margin=0.5, low_margin=0.15)
    assert calls["cell_type"].iloc[0] == "CD8 T"
    assert calls["confidence_level"].iloc[0] == "High"
    assert calls["cell_type"].iloc[1] == "Unassigned"


# ── suggest_thresholds ────────────────────────────────────────────────────────

def test_suggest_thresholds_returns_keys_and_ranges():
    rng = np.random.default_rng(0)
    n = 2000
    # a "matched" mode near 0.9 and an "unmatched" mode near 0.3 for two types
    best_a = np.concatenate([rng.normal(0.9, 0.03, n), rng.normal(0.3, 0.05, n)])
    other = np.concatenate([rng.normal(0.2, 0.03, n), rng.normal(0.3, 0.05, n)])
    scores = pd.DataFrame({"A": np.clip(best_a, 0, 1), "B": np.clip(other, 0, 1)})
    th = suggest_thresholds(scores)
    assert set(th) == {"min_score", "high_margin", "low_margin"}
    assert 0.5 <= th["min_score"] <= 0.65
    assert 0.05 <= th["low_margin"] <= 0.30
    assert th["high_margin"] >= th["low_margin"] + 0.10


def test_suggest_thresholds_single_type_zeroes_margins():
    scores = pd.DataFrame({"T cell": np.linspace(0.2, 0.95, 500)})
    th = suggest_thresholds(scores)
    assert th["high_margin"] == 0.0
    assert th["low_margin"] == 0.0
    assert th["min_score"] >= 0.5


def test_suggest_thresholds_empty_scores_defaults():
    th = suggest_thresholds(pd.DataFrame(index=range(10)))
    assert th == {"min_score": 0.5, "high_margin": 0.5, "low_margin": 0.15}


def test_suggest_thresholds_never_below_coin_flip():
    # all cells match well -> min_score must not drop below 0.5
    scores = pd.DataFrame({"A": np.full(500, 0.95), "B": np.full(500, 0.1)})
    th = suggest_thresholds(scores)
    assert th["min_score"] >= 0.5


# ── end-to-end ────────────────────────────────────────────────────────────────

def test_end_to_end_pipeline(bimodal_adata):
    posteriors, _ = compute_positive_posteriors(bimodal_adata, random_state=0)
    definitions = {
        "CD8 T": {"CD8": "+", "CD4": "-"},
        "CD4 T": {"CD8": "-", "CD4": "+"},
    }
    scores = score_cell_types(posteriors, definitions)
    calls = assign_cell_types(scores)

    assert len(calls) == bimodal_adata.n_obs
    obs = bimodal_adata.obs
    # leiden "0" blocks are CD8+CD4-, "1" blocks are CD8-CD4+
    is0 = obs["leiden"].to_numpy() == "0"
    top_for_0 = calls.loc[is0, "cell_type"].value_counts().index[0]
    top_for_1 = calls.loc[~is0, "cell_type"].value_counts().index[0]
    assert top_for_0 == "CD8 T"
    assert top_for_1 == "CD4 T"
