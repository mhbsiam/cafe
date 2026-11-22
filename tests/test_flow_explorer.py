import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "CAFE", "src", "cafe_app", "tools"))

from flow_explorer import asinh_transform, decade_ticks, density_rank, export_flow_plot, quadrant_shares  # noqa: E402


def test_decade_ticks_label_negative_zero_and_positive_decades():
    vals, labels = decade_ticks(-5_000, 200_000, cofactor=150)
    assert labels == ["−10³", "0", "10³", "10⁴", "10⁵"]
    assert np.allclose(vals, asinh_transform([-1e3, 0, 1e3, 1e4, 1e5], 150))


def test_density_rank_puts_dense_cluster_above_scattered_cells():
    rng = np.random.default_rng(0)
    dense = rng.normal(0, 0.05, size=(5_000, 2))
    sparse = rng.uniform(-3, 3, size=(200, 2))
    xy = np.vstack([dense, sparse])
    cell_rank, grid, xc, yc = density_rank(xy[:, 0], xy[:, 1], (-3, 3), (-3, 3))
    assert grid.shape == (len(yc), len(xc))
    assert 0 < cell_rank.min() and cell_rank.max() == 1.0
    assert np.median(cell_rank[:5_000]) > np.median(cell_rank[5_000:])


def test_quadrant_shares_split_cells_by_both_thresholds():
    x = np.array([10, 10, 5000, 5000, 5000, 1000])
    y = np.array([10, 5000, 10, 5000, 5000, 1000])  # last cell sits on both lines: counts as high
    q = quadrant_shares(x, y, 1000, 1000)
    assert q == {"bottom-left": 100 / 6, "top-left": 100 / 6, "bottom-right": 100 / 6, "top-right": 300 / 6}
    assert quadrant_shares(x, y, x_thr=1000) == {"left": 200 / 6, "right": 400 / 6}
    assert quadrant_shares(x, y) == {}


def test_exports_are_6_inches_square():
    import io
    import re

    import plotly.graph_objects as go
    from PIL import Image

    fig = go.Figure([go.Scattergl(x=[0, 1], y=[0, 1], mode="markers"),
                     go.Scattergl(x=[0.5], y=[0.5], mode="markers", marker=dict(opacity=0))])
    img = Image.open(io.BytesIO(export_flow_plot(fig, "png")))
    assert img.size == (1800, 1800)
    assert all(round(d) == 300 for d in img.info["dpi"])

    pdf = export_flow_plot(fig, "pdf")
    assert pdf.startswith(b"%PDF")
    box = [float(v) for v in re.search(rb"/MediaBox\s*\[([^\]]*)\]", pdf).group(1).split()]
    assert box[2:] == [432.0, 432.0]  # 6 in × 72 pt
    assert b"/Subtype /Image" not in pdf  # points are vector, not an embedded WebGL bitmap
