"""Interactive flow plot: one marker against another on an arcsinh (biexponential-like) scale, in FlowJo's six styles; never mutates adata."""
import io
import re

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from PIL import Image
from scipy.ndimage import gaussian_filter

from _plot_export import plotly_image_bytes
from embedding_explorer import _DISPLAY_CAP, _render_selection_summary
from utils import RANDOM_STATE, get_cluster_col, safe_flatten, safe_sort_clusters

STYLES = ["Pseudocolor", "Dot plot", "Contour", "Density", "Zebra", "Heatmap statistic"]
_POINT_STYLES = {"Pseudocolor", "Dot plot", "Heatmap statistic"}  # every cell is a point, so lasso works
_BINS = 200
_OUTLIER_FRACTION = 0.02  # lowest-density 2% of cells are drawn as dots outside contours / density
_PSEUDOCOLOR = [[0.0, "#1f2bff"], [0.3, "#00b7ff"], [0.55, "#18d64b"], [0.8, "#ffd400"], [1.0, "#ff1a1a"]]
# Colour scales for the coloured styles. Light-starting scales (Blues, Greys…) are left out: their low end vanishes on white.
_COLOR_SCALES = {"FlowJo classic": _PSEUDOCOLOR, "Viridis": "Viridis", "Plasma": "Plasma", "Inferno": "Inferno",
                 "Magma": "Magma", "Cividis": "Cividis", "Turbo": "Turbo", "Jet": "Jet", "Portland": "Portland"}
# White below the outlier level (smooth outer edge, unlike NaN masking), then light grey to black.
_GREYS = [[0.0, "#ffffff"], [0.02, "#ffffff"], [0.0201, "#e4e4e4"], [1.0, "#000000"]]
_SUPERSCRIPT = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def asinh_transform(values, cofactor):
    """Arcsinh scale used by spectral cytometry: linear around 0, logarithmic for large |values|."""
    return np.arcsinh(np.asarray(values, dtype=float) / cofactor)


def decade_ticks(lo, hi, cofactor):
    """Tick positions (in arcsinh space) and FlowJo-style labels (−10³, 0, 10³ …) for raw values in [lo, hi]."""
    raw = [-(10.0 ** e) for e in range(7, 1, -1)] + [0.0] + [10.0 ** e for e in range(2, 8)]
    # Decades inside the linear zone (|v| < 2 × cofactor) would crowd the 0 tick, so skip them.
    raw = [v for v in raw if lo <= v <= hi and (v == 0 or abs(v) >= 2 * cofactor)]
    labels = ["0" if v == 0 else ("−" if v < 0 else "") + "10" + str(int(np.log10(abs(v)))).translate(_SUPERSCRIPT)
              for v in raw]
    return asinh_transform(raw, cofactor), labels


def density_rank(x, y, x_range, y_range):
    """Smoothed 2D density on a grid, expressed as 'fraction of cells in sparser places' (0–1).

    Equal steps of this rank hold equal shares of cells, which gives FlowJo's equal-probability contours.
    Returns (rank per cell, rank grid with shape [y, x], x bin centres, y bin centres).
    """
    hist, xe, ye = np.histogram2d(x, y, bins=_BINS, range=[x_range, y_range])
    hist = gaussian_filter(hist, sigma=3)  # FlowJo-smooth contours; less smoothing gives ragged edges
    ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, _BINS - 1)
    iy = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, _BINS - 1)
    cell_density = hist[ix, iy]
    sorted_density = np.sort(cell_density)
    n = max(len(sorted_density), 1)
    cell_rank = np.searchsorted(sorted_density, cell_density, side="right") / n
    grid_rank = np.searchsorted(sorted_density, hist, side="right") / n
    centres = lambda e: (e[:-1] + e[1:]) / 2
    return cell_rank, grid_rank.T, centres(xe), centres(ye)


def _axis_range(t):
    lo, hi = np.percentile(t, [0.05, 99.95])
    pad = 0.05 * (hi - lo or 1.0)
    return float(lo - pad), float(hi + pad)


def _axis(title, rng, cofactor):
    tickvals, ticktext = decade_ticks(np.sinh(rng[0]) * cofactor, np.sinh(rng[1]) * cofactor, cofactor)
    return dict(title=title, range=list(rng), tickvals=tickvals, ticktext=ticktext, ticks="outside",
                showgrid=False, zeroline=False, showline=True, linecolor="#0b1a20", mirror=True)


def _cell_filter(adata, cluster_col, key_prefix):
    """'Show cells from' control: returns a boolean mask over adata rows."""
    columns = [c for c in (cluster_col, "Group", "SampleID") if c in adata.obs.columns]
    by = st.selectbox("Show cells from", ["All cells"] + columns, key=f"{key_prefix}_filter_by")
    if by == "All cells":
        return np.ones(adata.n_obs, dtype=bool)
    col = adata.obs[by].astype(str)
    chosen = st.multiselect(by, safe_sort_clusters(col.unique()), key=f"{key_prefix}_filter_values",
                            placeholder=f"Choose {by} values")
    return col.isin(chosen).to_numpy() if chosen else np.ones(adata.n_obs, dtype=bool)


def quadrant_shares(x_raw, y_raw, x_thr=None, y_thr=None):
    """Percent of cells on each side of the threshold(s): keys like 'left', 'top-right'; a value on the line counts as high."""
    sides = []
    if x_thr is not None:
        hi = x_raw >= x_thr
        sides.append({"left": ~hi, "right": hi})
    if y_thr is not None:
        hi = y_raw >= y_thr
        sides.append({"bottom": ~hi, "top": hi})
    if len(sides) == 2:
        parts = {f"{v}-{h}": sides[0][h] & sides[1][v] for h in sides[0] for v in sides[1]}
    else:
        parts = sides[0] if sides else {}
    n = max(len(x_raw), 1)
    return {k: 100.0 * m.sum() / n for k, m in parts.items()}


def _add_threshold_lines(fig, x_raw, y_raw, x_thr, y_thr, cofactor):
    """Draw the vertical/horizontal threshold lines and label each region with its share of cells."""
    line = dict(color="#000000", width=1.5)
    if x_thr is not None:
        fig.add_vline(x=float(asinh_transform(x_thr, cofactor)), line=line)
    if y_thr is not None:
        fig.add_hline(y=float(asinh_transform(y_thr, cofactor)), line=line)
    corner = {"left": 0.02, "right": 0.98, "bottom": 0.02, "top": 0.98}
    for name, pct in quadrant_shares(x_raw, y_raw, x_thr, y_thr).items():
        v, h = (name.split("-") if "-" in name else
                (("top", name) if name in ("left", "right") else (name, "left")))
        fig.add_annotation(x=corner[h], y=corner[v], xref="x domain", yref="y domain", text=f"{pct:.1f}%",
                           showarrow=False, xanchor=h, yanchor=v,
                           font=dict(size=13, color="#000000"), bgcolor="rgba(255,255,255,0.8)")


EXPORT_FORMATS = {"PNG (300 DPI)": "png", "PDF (vector)": "pdf"}


def export_flow_plot(fig, fmt, inches=6.0):
    """Render fig as a square `inches` figure: 'png' at 300 DPI (pixels plus DPI tag) or 'pdf' as vector; None if export is unavailable."""
    layout = go.Layout(fig.layout)
    # Room for tick labels and titles; no selection/drag state in a static file; black text for publication.
    layout.update(margin=dict(l=80, r=30, t=30, b=70, pad=0), dragmode=False, uirevision=None,
                  font=dict(color="#000000"))
    if fmt == "png":
        png = plotly_image_bytes(go.Figure(fig.data, layout), "png", width=600, height=600, scale=300 * inches / 600)
        if png is None:
            return None
        buf = io.BytesIO()
        Image.open(io.BytesIO(png)).save(buf, format="PNG", dpi=(300, 300))
        return buf.getvalue()
    # PDF: WebGL points would be embedded as a ~190 DPI bitmap, so redraw them as vector markers.
    # The invisible selection layer (contour styles) is dropped: it is most of the cells and draws nothing.
    traces = []
    for t in fig.data:
        if t.type == "scattergl":
            if t.marker.opacity == 0:
                continue
            t = go.Scatter({k: v for k, v in t.to_plotly_json().items() if k != "type"})
        traces.append(t)
    # Kaleido sizes PDF pages at 0.75 pt per px: 600 px × scale → inches × 72 pt, same layout as the PNG.
    return plotly_image_bytes(go.Figure(traces, layout), "pdf", width=600, height=600, scale=inches * 72 / 0.75 / 600)


def render_flow_explorer(adata, key_prefix="flow_explorer"):
    """Render the two-marker flow plot with lasso selection (point styles); never mutates adata or writes files."""
    markers = list(map(str, adata.var_names))
    if len(markers) < 2:
        st.info("A flow plot needs at least two markers.")
        return
    cluster_col = get_cluster_col(adata)

    row1 = st.columns(3)
    x_marker = row1[0].selectbox("X axis", markers, index=0, key=f"{key_prefix}_x")
    y_marker = row1[1].selectbox("Y axis", markers, index=1, key=f"{key_prefix}_y")
    style = row1[2].selectbox("Plot style", STYLES, key=f"{key_prefix}_style")

    row2 = st.columns(4)
    with row2[0]:
        mask = _cell_filter(adata, cluster_col, key_prefix)
    cofactor = row2[1].number_input(
        "Arcsinh cofactor", min_value=1, max_value=100_000, value=150, step=50, key=f"{key_prefix}_cofactor",
        help="Width of the linear zone around 0. Larger values squeeze dim and negative events together.",
    )
    color_marker, dot_size, colorscale = None, 2, None
    if style in ("Pseudocolor", "Heatmap statistic"):
        names = list(_COLOR_SCALES)
        # Each style remembers its own choice; Pseudocolor starts FlowJo-like, Heatmap statistic on Viridis.
        scale_name = row2[2].selectbox("Color scale", names, index=names.index(
            "FlowJo classic" if style == "Pseudocolor" else "Viridis"), key=f"{key_prefix}_scale_{style}")
        colorscale = _COLOR_SCALES[scale_name]
    if style == "Heatmap statistic":
        color_marker = row2[3].selectbox("Color by marker", markers, index=min(2, len(markers) - 1),
                                         key=f"{key_prefix}_color")
    elif style in _POINT_STYLES:
        dot_size = row2[3].slider("Point size", 1, 8, 2, key=f"{key_prefix}_size")

    # Threshold lines (FlowJo quadrant style), in raw units; ± steppers move them.
    # Each line is one control: a switch with its stepper directly below (keeps them aligned).
    row3 = st.columns(2)
    with row3[0]:
        show_v = st.toggle(f"Vertical line · {x_marker}", key=f"{key_prefix}_vline")
        x_thr = st.number_input(f"{x_marker} threshold", value=1000, step=100, disabled=not show_v,
                                label_visibility="collapsed", key=f"{key_prefix}_xthr")
    with row3[1]:
        show_h = st.toggle(f"Horizontal line · {y_marker}", key=f"{key_prefix}_hline")
        y_thr = st.number_input(f"{y_marker} threshold", value=1000, step=100, disabled=not show_h,
                                label_visibility="collapsed", key=f"{key_prefix}_ythr")

    rows = np.flatnonzero(mask)
    if rows.size < 3:
        st.info("Too few cells match this filter to draw a flow plot.")
        return

    x_raw = safe_flatten(adata[:, x_marker].X)[rows]
    y_raw = safe_flatten(adata[:, y_marker].X)[rows]
    x, y = asinh_transform(x_raw, cofactor), asinh_transform(y_raw, cofactor)
    x_range, y_range = _axis_range(x), _axis_range(y)
    # Density uses every filtered cell; only the drawn points are capped.
    cell_rank, grid_rank, xc, yc = density_rank(x, y, x_range, y_range)

    # Every style carries the cells as points so lasso/box selection works over contours and density too.
    shown = np.arange(rows.size)
    if shown.size > _DISPLAY_CAP:
        shown = np.sort(np.random.default_rng(RANDOM_STATE).choice(shown, _DISPLAY_CAP, replace=False))
        st.caption(f":material/bolt: Drawing a random **{_DISPLAY_CAP:,}**-cell sample of {rows.size:,}. "
                   "Density and contours use every cell.")

    obs = adata.obs.iloc[rows[shown]]
    customdata = np.column_stack([
        rows[shown].astype(object),
        obs[cluster_col].astype(str).to_numpy() if cluster_col in obs.columns else np.full(shown.size, ""),
        obs["Group"].astype(str).to_numpy() if "Group" in obs.columns else np.full(shown.size, ""),
        x_raw[shown], y_raw[shown],
    ])
    hover = (f"{x_marker} " + "%{customdata[3]:,.0f}<br>" + f"{y_marker} " + "%{customdata[4]:,.0f}<br>"
             "cluster %{customdata[1]} · group %{customdata[2]}<extra></extra>")

    fig = go.Figure()
    if style in _POINT_STYLES:
        if style == "Pseudocolor":
            marker = dict(color=cell_rank[shown], colorscale=colorscale, cmin=0, cmax=1)
        elif style == "Dot plot":
            marker = dict(color="#000000")
        else:
            c_raw = safe_flatten(adata[:, color_marker].X)[rows][shown]
            # Label the colour bar in the same decades as the axes, so the legend can be read.
            tickvals, ticktext = decade_ticks(c_raw.min(), c_raw.max(), cofactor)
            marker = dict(color=asinh_transform(c_raw, cofactor), colorscale=colorscale,
                          colorbar=dict(title=color_marker, thickness=12, len=0.6,
                                        tickvals=tickvals, ticktext=ticktext))
        marker.update(size=dot_size, line=dict(width=0))
        fig.add_trace(go.Scattergl(x=x[shown], y=y[shown], mode="markers", marker=marker,
                                   customdata=customdata, hovertemplate=hover, showlegend=False))
    else:
        levels = dict(start=_OUTLIER_FRACTION, end=0.92, size=0.1)  # each band holds ~10% of cells
        if style == "Density":
            fig.add_trace(go.Heatmap(x=xc, y=yc, z=grid_rank, colorscale=_GREYS, zmin=0, zmax=1,
                                     zsmooth="best", showscale=False, hoverinfo="skip"))
        elif style == "Contour":
            fig.add_trace(go.Contour(x=xc, y=yc, z=grid_rank, contours=dict(coloring="lines", **levels),
                                     line=dict(color="#000000", width=1), showscale=False, hoverinfo="skip"))
        else:  # Zebra: filled equal-probability bands with black edges
            fig.add_trace(go.Contour(x=xc, y=yc, z=grid_rank, contours=dict(coloring="fill", **levels),
                                     colorscale=_GREYS, zmin=0, zmax=1, line=dict(color="#000000", width=0.6),
                                     showscale=False, hoverinfo="skip"))
        # Two disjoint point layers: visible outlier dots outside the contours, and invisible
        # selectable cells inside them (so a lasso counts every cell exactly once).
        outlier = cell_rank[shown] <= _OUTLIER_FRACTION
        for part, opacity in ((~outlier, 0.0), (outlier, 1.0)):
            fig.add_trace(go.Scattergl(
                x=x[shown][part], y=y[shown][part], mode="markers",
                marker=dict(size=2, color="#000000", opacity=opacity, line=dict(width=0)),
                selected=dict(marker=dict(opacity=opacity)), unselected=dict(marker=dict(opacity=opacity * 0.3)),
                customdata=customdata[part], hovertemplate=hover if opacity else None,
                hoverinfo=None if opacity else "skip", showlegend=False))

    fig.update_layout(
        # pad=0: Streamlit's theme pads the frame 8px off the data area, leaving a gap at the threshold lines.
        height=620, margin=dict(l=8, r=8, t=8, b=8, pad=0), plot_bgcolor="#ffffff", paper_bgcolor="#ffffff",
        xaxis=_axis(x_marker, x_range, cofactor),
        yaxis=_axis(y_marker, y_range, cofactor),
        # Keep the view while the style changes; reset it when the axes change.
        uirevision=f"flow-{x_marker}-{y_marker}-{cofactor}",
        dragmode="lasso",
    )

    _add_threshold_lines(fig, x_raw, y_raw, x_thr if show_v else None, y_thr if show_h else None, cofactor)

    batch_corrected = "Batch" in adata.obs.columns and adata.obs["Batch"].nunique() > 1
    st.caption(
        f"{rows.size:,} cells · arcsinh scale, cofactor {cofactor:,} · "
        + ("values are **ComBat batch-corrected**, not the original CSV values."
           if batch_corrected else "values are the scaled CSV exports.")
    )

    event = st.plotly_chart(fig, key=f"{key_prefix}_chart", on_select="rerun", selection_mode=("box", "lasso"),
                            config={"displaylogo": False, "scrollZoom": True}, width="stretch")
    # Rendered on request: a vector PDF of 100k points takes ~25 s, a PNG a few seconds.
    exp = st.columns([2, 1, 3], vertical_alignment="bottom")
    label = exp[0].radio("Download format", list(EXPORT_FORMATS), horizontal=True, key=f"{key_prefix}_export_fmt")
    if exp[1].button("Prepare file", key=f"{key_prefix}_export"):
        fmt = EXPORT_FORMATS[label]
        with st.spinner(f"Rendering {label}…"):
            data = export_flow_plot(fig, fmt)
        if data is not None:
            name = re.sub(r"[^A-Za-z0-9._-]+", "_", f"flowplot_{x_marker}_vs_{y_marker}_{style}")
            exp[2].download_button(f"Download {label}, 6 × 6 in", data, file_name=f"{name}.{fmt}",
                                   mime="image/png" if fmt == "png" else "application/pdf",
                                   on_click="ignore", key=f"{key_prefix}_export_download")
    # Same read-out as the UMAP tab: cell count, cluster and group make-up, marker profile vs all cells.
    _render_selection_summary(adata, event, cluster_col)
