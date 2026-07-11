"""Interactive UMAP explorer: a live Plotly Scattergl companion to the static figures; additive, never writes files or mutates adata."""
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from scipy.sparse import issparse
import streamlit as st

from utils import (
    safe_flatten,
    get_cluster_col,
    safe_sort_clusters,
    RANDOM_STATE,
    MARKER_HIGH_THRESHOLD,
)

# Categorical obs columns worth colouring by, in the order researchers reach
# for them. Filtered against what's actually present at render time.
_CATEGORICAL_CANDIDATES = ["leiden", "cell_type", "Group", "SampleID"]

# Above this, the explorer samples down to this size (static figures still use every cell).
_DISPLAY_CAP = 100_000

# Continuous marker-expression scale: light warm grey -> signature red.
_MARKER_SCALE = [
    [0.0, "#e9e4e0"],
    [0.5, "#ff9e8f"],
    [1.0, "#d6203a"],
]


def _category_colors(categories, cmap_name):
    """Map category labels -> hex colours via a matplotlib colormap (mirrors Visualization.get_palette)."""
    cmap = plt.get_cmap(cmap_name)
    listed = getattr(cmap, "colors", None)
    n = len(categories)
    colors = {}
    for i, cat in enumerate(categories):
        if listed is not None and len(listed) >= n:
            colors[cat] = to_hex(listed[i])
        else:
            colors[cat] = to_hex(cmap(i / max(n - 1, 1)))
    return colors


def _base_layout(height):
    return dict(
        height=height,
        margin=dict(l=8, r=8, t=8, b=8),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False, scaleanchor="x", scaleratio=1),
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        legend=dict(
            title="",
            itemsizing="constant",
            bgcolor="rgba(255,255,255,0.75)",
            borderwidth=0,
            font=dict(size=11),
        ),
        # Constant uirevision preserves pan/zoom across recolour reruns.
        uirevision="cafe-umap-explorer",
        dragmode="lasso",
    )


def _build_categorical_figure(coords, labels, obs, cluster_col, cmap_name, dot_size, opacity, height, row_ids):
    categories = safe_sort_clusters(labels.unique())
    color_map = _category_colors(categories, cmap_name)

    group = obs["Group"].astype(str).values if "Group" in obs.columns else None
    sample = obs["SampleID"].astype(str).values if "SampleID" in obs.columns else None
    cluster = obs[cluster_col].astype(str).values if cluster_col in obs.columns else None

    fig = go.Figure()
    for cat in categories:
        mask = (labels.astype(str).values == str(cat))
        idx = np.where(mask)[0]  # positions within the (possibly sampled) view
        if idx.size == 0:
            continue
        # customdata[0] is the ORIGINAL adata row index (so selection is correct
        # even when the view is downsampled); the rest are hover fields.
        cd = np.column_stack([
            row_ids[idx].astype(object),
            (cluster[idx] if cluster is not None else np.array([""] * idx.size)),
            (group[idx] if group is not None else np.array([""] * idx.size)),
            (sample[idx] if sample is not None else np.array([""] * idx.size)),
        ])
        fig.add_trace(go.Scattergl(
            x=coords[idx, 0],
            y=coords[idx, 1],
            mode="markers",
            name=str(cat),
            marker=dict(size=dot_size, color=color_map[cat], opacity=opacity,
                        line=dict(width=0)),
            customdata=cd,
            hovertemplate=(
                "cluster %{customdata[1]}<br>"
                "group %{customdata[2]}<br>"
                "sample %{customdata[3]}"
                "<extra></extra>"
            ),
        ))
    fig.update_layout(**_base_layout(height))
    return fig


def _build_continuous_figure(coords, values, obs, marker, cluster_col, dot_size, opacity, height, row_ids):
    n = coords.shape[0]
    group = obs["Group"].astype(str).values if "Group" in obs.columns else np.array([""] * n)
    cluster = obs[cluster_col].astype(str).values if cluster_col in obs.columns else np.array([""] * n)
    # customdata[0] is the ORIGINAL adata row index (see categorical builder).
    cd = np.column_stack([row_ids.astype(object), cluster, group, values])

    fig = go.Figure(go.Scattergl(
        x=coords[:, 0],
        y=coords[:, 1],
        mode="markers",
        marker=dict(
            size=dot_size,
            color=values,
            colorscale=_MARKER_SCALE,
            opacity=opacity,
            line=dict(width=0),
            colorbar=dict(title=marker, thickness=12, len=0.6, x=1.0),
        ),
        customdata=cd,
        hovertemplate=(
            f"{marker} " + "%{customdata[3]:.1f}<br>"
            "cluster %{customdata[1]}<br>"
            "group %{customdata[2]}"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(**_base_layout(height))
    return fig


def _render_selection_summary(adata, event, cluster_col):
    """Turn a lasso/box selection into a live composition read-out."""
    points = []
    if event and isinstance(event, dict):
        points = event.get("selection", {}).get("points", []) or []

    if not points:
        st.caption(
            "Lasso or box-select cells on the plot to see how many you picked "
            "and what clusters and groups they come from."
        )
        return

    rows = []
    for p in points:
        cd = p.get("customdata")
        if cd:
            rows.append(int(cd[0]))
    if not rows:
        return

    sel = adata.obs.iloc[rows]
    total = len(sel)

    st.markdown(f"**{total:,}** cells selected  ·  {total / adata.n_obs:.1%} of all cells")

    cols = st.columns(2)
    with cols[0]:
        st.caption(f"By {cluster_col}")
        by_cluster = sel[cluster_col].astype(str).value_counts()
        by_cluster = by_cluster.reindex(safe_sort_clusters(by_cluster.index))
        st.bar_chart(by_cluster, height=180, color="#0f5070")
    if "Group" in sel.columns:
        with cols[1]:
            st.caption("By group")
            st.bar_chart(sel["Group"].astype(str).value_counts(), height=180, color="#3d8ba8")

    _render_marker_profile(adata, rows)


def _mean_expression(X):
    """Mean expression per column (marker), sparse- and dense-safe."""
    if issparse(X):
        return np.asarray(X.mean(axis=0)).ravel()
    return np.asarray(X).mean(axis=0)


def _densify(X):
    return np.asarray(X.toarray()) if issparse(X) else np.asarray(X)


def _background_stats(adata):
    """All-cell mean & median per marker, cached per session (median computed column-by-column to avoid densifying)."""
    key = (adata.n_obs, adata.n_vars, tuple(map(str, adata.var_names)))
    cached = st.session_state.get("_cafe_bg_stats")
    if cached and cached["key"] == key:
        return cached["mean"], cached["median"]

    mean = _mean_expression(adata.X)
    median = np.empty(adata.n_vars)
    for i in range(adata.n_vars):
        median[i] = float(np.median(safe_flatten(adata.X[:, i])))
    st.session_state["_cafe_bg_stats"] = {"key": key, "mean": mean, "median": median}
    return mean, median


def _render_marker_profile(adata, rows):
    """Per-marker expression profile of the selected cells vs the all-cell background, with positive-fraction annotations."""
    markers = list(adata.var_names)
    if not markers:
        return

    metric = st.radio(
        "Rank markers by",
        ["Median", "Mean", "Enrichment (Δ vs all)"],
        horizontal=True,
        key="umap_explorer_marker_metric",
    )

    sel = _densify(adata.X[rows])
    sel_mean = sel.mean(axis=0)
    sel_median = np.median(sel, axis=0)
    pct_pos = (sel > MARKER_HIGH_THRESHOLD).mean(axis=0) * 100.0
    bg_mean, bg_median = _background_stats(adata)

    if metric == "Mean":
        sel_vals, ref_vals, xtitle, diverging = sel_mean, bg_mean, "Mean expression", False
    elif metric.startswith("Enrichment"):
        # Difference of medians vs background. A ratio/log2 "fold-change" is
        # unsafe here: raw spectral-flow values are frequently negative.
        sel_vals = sel_median - bg_median
        ref_vals, xtitle, diverging = None, "Δ median vs all cells", True
    else:  # Median (default)
        sel_vals, ref_vals, xtitle, diverging = sel_median, bg_median, "Median expression", False

    order = np.argsort(sel_vals)  # ascending -> largest ends on top of the h-bar
    labels = [markers[i] for i in order]
    sel_ord = sel_vals[order]
    pct_ord = pct_pos[order]
    text = [f"{p:.0f}%+" for p in pct_ord]

    top3 = ", ".join(markers[i] for i in np.argsort(sel_vals)[::-1][:3])
    thr = int(MARKER_HIGH_THRESHOLD)
    st.caption(
        f"Marker profile of the {len(rows):,} selected cells · top: **{top3}** · "
        f"“%+” = share of selected cells above the positive threshold ({thr:,})."
    )

    fig = go.Figure()
    if diverging:
        colors = ["#ff4b4b" if v >= 0 else "#4b7bff" for v in sel_ord]
        fig.add_trace(go.Bar(
            y=labels, x=sel_ord, orientation="h", marker=dict(color=colors),
            text=text, textposition="outside", cliponaxis=False,
            customdata=pct_ord,
            hovertemplate="%{y}: Δ %{x:.1f}<br>%{customdata:.0f}% positive<extra></extra>",
        ))
    else:
        ref_ord = ref_vals[order]
        fig.add_trace(go.Bar(
            y=labels, x=ref_ord, orientation="h", name="All cells",
            marker=dict(color="#d8d4d0"),
            hovertemplate="%{y}: %{x:.1f}<extra>All cells</extra>",
        ))
        fig.add_trace(go.Bar(
            y=labels, x=sel_ord, orientation="h", name="Selected",
            marker=dict(color="#ff4b4b"), opacity=0.9,
            text=text, textposition="outside", cliponaxis=False,
            customdata=pct_ord,
            hovertemplate="%{y}: %{x:.1f}<br>%{customdata:.0f}% positive<extra>Selected</extra>",
        ))

    fig.update_layout(
        barmode="overlay",
        bargap=0.25,
        showlegend=not diverging,
        height=max(240, 24 * len(markers) + 90),
        margin=dict(l=8, r=8, t=8, b=8),
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        xaxis=dict(title=xtitle, zeroline=True, zerolinecolor="#c9c4bf",
                   gridcolor="#f0ece8"),
        yaxis=dict(automargin=True),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="right", x=1.0,
                    font=dict(size=11)),
    )
    st.plotly_chart(fig, key="umap_explorer_marker_profile", width="stretch",
                    config={"displaylogo": False})


def render_embedding_explorer(adata, key_prefix="umap_explorer"):
    """Render the interactive embedding explorer; no-ops with a hint if no UMAP yet. Never mutates adata or writes files."""
    if "X_umap" not in adata.obsm:
        st.info("Compute a UMAP embedding first. The interactive explorer will appear here.")
        return

    coords_all = np.asarray(adata.obsm["X_umap"])[:, :2]
    cluster_col = get_cluster_col(adata)
    n_total = adata.n_obs

    # Downsample the view past the cap (seeded); row_ids map back to true cells via customdata.
    if n_total > _DISPLAY_CAP:
        rng = np.random.default_rng(RANDOM_STATE)
        row_ids = np.sort(rng.choice(n_total, _DISPLAY_CAP, replace=False))
        st.caption(
            f":material/bolt: Showing a random **{_DISPLAY_CAP:,}**-cell sample of "
            f"{n_total:,} for a responsive view. The static plots below use every cell."
        )
    else:
        row_ids = np.arange(n_total)

    coords = coords_all[row_ids]
    obs = adata.obs.iloc[row_ids]

    cat_options = [c for c in _CATEGORICAL_CANDIDATES if c in obs.columns]
    marker_options = [f"marker: {m}" for m in adata.var_names]
    options = cat_options + marker_options

    # Default to the cluster identity researchers see in the static plots.
    default_idx = cat_options.index(cluster_col) if cluster_col in cat_options else 0

    ctrl = st.columns([2, 1, 1])
    with ctrl[0]:
        color_by = st.selectbox("Color by", options, index=default_idx,
                                key=f"{key_prefix}_colorby")
    with ctrl[1]:
        dot_size = st.slider("Point size", 1, 12, 4, key=f"{key_prefix}_size")
    with ctrl[2]:
        opacity = st.slider("Opacity", 0.1, 1.0, 0.75, 0.05, key=f"{key_prefix}_opacity")

    height = 620
    if color_by.startswith("marker: "):
        marker = color_by[len("marker: "):]
        values = safe_flatten(adata[:, marker].X)[row_ids]
        fig = _build_continuous_figure(coords, values, obs, marker, cluster_col,
                                       dot_size, opacity, height, row_ids)
    else:
        labels = obs[color_by]
        fig = _build_categorical_figure(coords, labels, obs, cluster_col,
                                        _cmap_for(color_by), dot_size, opacity, height, row_ids)

    event = st.plotly_chart(
        fig,
        key=f"{key_prefix}_chart",
        on_select="rerun",
        selection_mode=("box", "lasso"),
        config={"displaylogo": False, "scrollZoom": True},
        width="stretch",
    )

    _render_selection_summary(adata, event, cluster_col)


def _cmap_for(color_by):
    """Choose a categorical colormap that echoes the static UMAP defaults."""
    # 'Group' is low-cardinality (usually 2); a punchy qualitative map reads
    # better there. Clusters/samples can be many, so use the 40-colour set.
    if color_by == "Group":
        return "Set1"
    return "tab20c"
