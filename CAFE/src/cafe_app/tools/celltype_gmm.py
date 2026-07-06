"""Probabilistic marker-based cell-type annotation via per-marker Gaussian mixtures.

FlowJo-scaled expression values place each marker's negative and positive
populations as two roughly-Gaussian humps.  Fitting a 2-component Gaussian
Mixture Model per marker recovers the pos/neg gate and yields, for every cell,
the posterior probability of belonging to the positive component.  That
posterior doubles as a per-marker confidence.

Gating strategy (per marker):
    1. Fit a **pooled** GMM across all cells — stabilises minority-positive
       markers (e.g. CD8a within all leukocytes) whose positive shoulder is too
       small to fit reliably in a single sample.
    2. For each sample, prefer a **per-sample** fit (absorbs intensity drift)
       when it is genuinely bimodal, otherwise fall back to the pooled gate,
       and only neutralise the marker (posterior 0.5) when it is unimodal even
       pooled.  This graceful degradation stops a hard-to-fit marker from being
       silently discarded.

A fit is considered bimodal ("informative") only when a 2-component fit beats a
1-component fit by BIC, *and* both components carry a non-trivial weight, *and*
their means are separated (Cohen's d over the average component SD).  The
average-SD denominator avoids the failure of a pooled-variance metric, which is
inflated by a broad positive shoulder and under-reports real separation.

These functions are deliberately Streamlit-free so they can be unit tested.

Pipeline:
    posteriors, diag = compute_positive_posteriors(adata)
    scores           = score_cell_types(posteriors, definitions)
    calls            = assign_cell_types(scores)
"""
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.mixture import GaussianMixture

_MIN_SEPARATION = 0.8   # Cohen's d (average-SD) between the two component means
_MIN_WEIGHT = 0.02      # each component must hold >= 2% of cells to count
_MIN_CELLS = 10


def _as_1d(values):
    """Return a dense 1-D float array from a dense or sparse column."""
    if issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=float).ravel()


def _cutoff_from_proba(x, proba):
    """First value where the positive posterior crosses 0.5 (else NaN)."""
    order = np.argsort(x)
    xs, ps = x[order], proba[order]
    crossings = np.where(np.diff(np.sign(ps - 0.5)) != 0)[0]
    return float(xs[crossings[0]]) if crossings.size else float("nan")


def fit_marker_gmm(values, random_state=0, max_fit=50000):
    """Fit a 2-component GMM to one marker; return posteriors, info, and model.

    Parameters
    ----------
    values : array-like or sparse column
        Expression values for a single marker (one entry per cell).
    random_state : int
        Seed for GMM initialisation and fit subsampling (reproducibility).
    max_fit : int
        Cap on the number of cells used to *fit* the mixture.  All cells are
        still scored via ``predict_proba``.

    Returns
    -------
    pos_posterior : np.ndarray, shape (n,)
        Probability each cell belongs to the positive (higher-mean) component.
    info : dict
        ``cutoff``, ``means`` (neg, pos), ``weights`` (neg, pos),
        ``separation`` (Cohen's d), ``informative`` (bool) and ``model`` — a
        reusable ``{"gmm", "pos"}`` bundle (``None`` for a degenerate fit) so a
        pooled fit can be applied to individual samples.
    """
    x = _as_1d(values)
    n = x.size
    degenerate = {
        "cutoff": float("nan"),
        "means": (float("nan"), float("nan")),
        "weights": (float("nan"), float("nan")),
        "separation": 0.0,
        "informative": False,
        "model": None,
    }
    if n < _MIN_CELLS or not np.isfinite(x).all() or np.ptp(x) < 1e-9 or x.std() < 1e-9:
        return np.full(n, 0.5), dict(degenerate)

    xcol = x.reshape(-1, 1)
    if n > max_fit:
        rng = np.random.default_rng(random_state)
        fit_x = xcol[rng.choice(n, size=max_fit, replace=False)]
    else:
        fit_x = xcol

    gmm2 = GaussianMixture(n_components=2, n_init=3, random_state=random_state).fit(fit_x)
    gmm1 = GaussianMixture(n_components=1, random_state=random_state).fit(fit_x)

    means = gmm2.means_.ravel()
    sds = np.sqrt(gmm2.covariances_.ravel())
    weights = gmm2.weights_.ravel()
    pos = int(np.argmax(means))          # positive = higher-mean component
    neg = 1 - pos

    proba = gmm2.predict_proba(xcol)[:, pos]

    # Cohen's d over the *average* component SD: not dominated by a broad
    # positive shoulder the way sqrt(mean(variance)) is.
    separation = float(abs(means[pos] - means[neg]) / (0.5 * (sds[pos] + sds[neg]) + 1e-9))
    bimodal = gmm2.bic(fit_x) < gmm1.bic(fit_x)
    informative = bool(
        bimodal and weights.min() >= _MIN_WEIGHT and separation >= _MIN_SEPARATION
    )

    info = {
        "cutoff": _cutoff_from_proba(x, proba),
        "means": (float(means[neg]), float(means[pos])),
        "weights": (float(weights[neg]), float(weights[pos])),
        "separation": separation,
        "informative": informative,
        "model": {"gmm": gmm2, "pos": pos},
    }
    return proba, info


def _score_with_model(values, model):
    """Positive-component posteriors for ``values`` under a fitted GMM bundle."""
    x = _as_1d(values).reshape(-1, 1)
    return model["gmm"].predict_proba(x)[:, model["pos"]]


def compute_positive_posteriors(adata, sample_key="SampleID", random_state=0,
                                max_fit=50000, progress=None):
    """Per-sample, per-marker positive posteriors for every cell.

    For each marker a pooled GMM is fit first; each sample then uses its own fit
    when bimodal, the pooled fit as a fallback, and 0.5 only when the marker is
    unimodal even pooled.

    Returns
    -------
    posteriors : pd.DataFrame, shape (n_cells, n_markers)
        Indexed like ``adata.obs_names``, columns are ``adata.var_names``.
    diagnostics : pd.DataFrame
        One row per (SampleID, marker) with ``cutoff`` (of the gate actually
        used), per-sample ``separation`` / ``informative``, and ``source`` in
        {``per-sample``, ``pooled``, ``neutral``}.
    """
    markers = list(adata.var_names)
    obs_names = adata.obs_names.astype(str)
    posteriors = pd.DataFrame(0.5, index=obs_names, columns=markers, dtype=float)

    X = adata.X
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    if sample_key in adata.obs.columns:
        sample_values = adata.obs[sample_key].astype(str).to_numpy()
        samples = list(pd.unique(sample_values))
    else:
        sample_values = np.array(["__all__"] * adata.n_obs)
        samples = ["__all__"]
    single_sample = len(samples) == 1

    diagnostics = []
    total = max(len(markers) * (len(samples) + 1), 1)
    step = 0
    for m_idx, marker in enumerate(markers):
        col_full = X[:, m_idx]

        # Pooled fit — stabiliser for minority-positive markers and fallback gate.
        _, pooled_info = fit_marker_gmm(col_full, random_state=random_state, max_fit=max_fit)
        step += 1
        if progress is not None:
            progress(step / total)

        for sample in samples:
            row_pos = np.where(sample_values == sample)[0]
            vals = col_full[row_pos]

            if single_sample:
                s_info = pooled_info
                s_proba = (
                    _score_with_model(vals, pooled_info["model"])
                    if pooled_info["model"] is not None
                    else np.full(vals.size, 0.5)
                )
            else:
                s_proba, s_info = fit_marker_gmm(
                    vals, random_state=random_state, max_fit=max_fit
                )

            if s_info["informative"]:
                proba, source, used_cutoff = s_proba, "per-sample", s_info["cutoff"]
            elif pooled_info["informative"] and pooled_info["model"] is not None:
                proba = _score_with_model(vals, pooled_info["model"])
                source, used_cutoff = "pooled", pooled_info["cutoff"]
            else:
                proba = np.full(vals.size, 0.5)
                source, used_cutoff = "neutral", float("nan")

            posteriors.iloc[row_pos, m_idx] = proba
            diagnostics.append({
                "SampleID": sample,
                "marker": marker,
                "cutoff": used_cutoff,
                "separation": s_info["separation"],
                "informative": s_info["informative"],
                "source": source,
            })
            step += 1
            if progress is not None:
                progress(step / total)

    return posteriors, pd.DataFrame(diagnostics)


def score_cell_types(posteriors, definitions):
    """Geometric-mean match score for each cell type.

    Parameters
    ----------
    posteriors : pd.DataFrame
        Output of :func:`compute_positive_posteriors`.
    definitions : dict[str, dict[str, str]]
        ``{cell_type: {marker: "+"|"-"}}``.  Markers absent from a type's dict
        (or with a blank state) are ignored for that type.

    Returns
    -------
    scores : pd.DataFrame, shape (n_cells, n_types)
        Geometric mean of per-marker match probabilities.  The geometric mean
        normalises across types defined by different numbers of markers (unlike
        a raw product).  Types with no usable markers are skipped; if none
        remain the frame is empty.
    """
    eps = 1e-9
    scores = {}
    for cell_type, marker_states in definitions.items():
        used = [
            (m, s) for m, s in marker_states.items()
            if s in ("+", "-") and m in posteriors.columns
        ]
        if not used:
            continue
        log_terms = []
        for marker, state in used:
            p = posteriors[marker].to_numpy()
            match = p if state == "+" else (1.0 - p)
            log_terms.append(np.log(np.clip(match, eps, 1.0)))
        scores[cell_type] = np.exp(np.mean(np.vstack(log_terms), axis=0))
    return pd.DataFrame(scores, index=posteriors.index)


def assign_cell_types(scores, min_score=0.5, high_margin=0.5, low_margin=0.15):
    """Assign each cell to its best-scoring type with a confidence level.

    Parameters
    ----------
    scores : pd.DataFrame
        Output of :func:`score_cell_types`.
    min_score : float
        Cells whose best raw score is below this match no type well and become
        ``"Unassigned"`` / ``Ambiguous``.
    high_margin, low_margin : float
        Thresholds on the margin between the best and runner-up type
        (measured on scores normalised to pseudo-probabilities).

    Returns
    -------
    calls : pd.DataFrame
        Columns ``cell_type``, ``best_score``, ``confidence_score`` (margin),
        ``confidence_level`` in {``High``, ``Low``, ``Ambiguous``}.
    """
    index = scores.index
    n = len(index)
    if scores.shape[1] == 0:
        return pd.DataFrame({
            "cell_type": pd.Series(["Unassigned"] * n, index=index),
            "best_score": pd.Series(np.zeros(n), index=index),
            "confidence_score": pd.Series(np.zeros(n), index=index),
            "confidence_level": pd.Series(["Ambiguous"] * n, index=index),
        })

    types = np.asarray(scores.columns)
    S = scores.to_numpy()

    best_idx = np.argmax(S, axis=1)
    best_score = S[np.arange(n), best_idx]

    # Normalise to pseudo-probabilities for the best-vs-runner-up margin.
    row_sum = S.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    P_sorted = np.sort(S / row_sum, axis=1)
    p_best = P_sorted[:, -1]
    p_second = P_sorted[:, -2] if S.shape[1] > 1 else np.zeros(n)
    margin = p_best - p_second

    cell_type = types[best_idx].astype(object)
    level = np.empty(n, dtype=object)

    matched = best_score >= min_score
    cell_type[~matched] = "Unassigned"
    level[~matched] = "Ambiguous"
    level[matched & (margin < low_margin)] = "Ambiguous"
    level[matched & (margin >= high_margin)] = "High"
    level[matched & (margin >= low_margin) & (margin < high_margin)] = "Low"

    return pd.DataFrame({
        "cell_type": pd.Series(cell_type, index=index),
        "best_score": pd.Series(best_score, index=index),
        "confidence_score": pd.Series(margin, index=index),
        "confidence_level": pd.Series(level, index=index),
    })


def _otsu_threshold(values, bins=256):
    """Otsu's method: the value that best splits a distribution into two modes.

    Returns the split point, or ``None`` when the data is degenerate (constant /
    too few points / no between-class separation).
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 8 or np.ptp(v) < 1e-9:
        return None
    hist, edges = np.histogram(v, bins=bins)
    hist = hist.astype(float)
    total = hist.sum()
    if total == 0:
        return None
    centers = (edges[:-1] + edges[1:]) / 2.0
    w0 = np.cumsum(hist)
    w1 = total - w0
    cum = np.cumsum(hist * centers)
    m_total = cum[-1]
    valid = (w0 > 0) & (w1 > 0)
    between = np.zeros_like(w0, dtype=float)
    mu0 = cum[valid] / w0[valid]
    mu1 = (m_total - cum[valid]) / w1[valid]
    between[valid] = w0[valid] * w1[valid] * (mu0 - mu1) ** 2
    if between.max() <= 0:
        return None
    return float(centers[int(np.argmax(between))])


def suggest_thresholds(scores):
    """Data-driven confidence thresholds from the score / margin distributions.

    ``min_score`` is placed at the valley of the per-cell best-score distribution
    (never below the coin-flip level 0.5, since a fully-neutral cell scores 0.5).
    The margin thresholds are placed at valleys of the best-vs-runner-up margin
    distribution.  Everything is clamped to sane ranges and falls back to fixed
    defaults when a distribution is unimodal.  With a single defined type there
    is no runner-up, so the margins are set to 0 (every matched cell is High).

    Returns
    -------
    dict with keys ``min_score``, ``high_margin``, ``low_margin``.
    """
    defaults = {"min_score": 0.5, "high_margin": 0.5, "low_margin": 0.15}
    if scores.shape[1] == 0:
        return dict(defaults)

    S = scores.to_numpy()
    best = S.max(axis=1)

    ms = _otsu_threshold(best)
    # Only raise min_score above the 0.5 floor if a real "unmatched" mode exists.
    min_score = float(np.clip(ms, 0.5, 0.65)) if ms is not None else 0.5

    if S.shape[1] == 1:
        return {"min_score": min_score, "high_margin": 0.0, "low_margin": 0.0}

    row_sum = S.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    P_sorted = np.sort(S / row_sum, axis=1)
    margin = P_sorted[:, -1] - P_sorted[:, -2]

    lm = _otsu_threshold(margin)
    low_margin = float(np.clip(lm, 0.05, 0.30)) if lm is not None else 0.15

    upper = margin[margin >= low_margin]
    hm = _otsu_threshold(upper) if upper.size > 20 else None
    if hm is None:
        hm = float(np.quantile(margin, 0.60))
    high_margin = float(np.clip(hm, low_margin + 0.10, 0.85))

    return {"min_score": min_score, "high_margin": high_margin, "low_margin": low_margin}
