"""Page-level (AppTest) checks for the stage-3 Visualization.py fixes: dataset-shape guards and always-failing options."""
import json
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from streamlit.testing.v1 import AppTest

CAFE_APP = Path(__file__).resolve().parents[1] / "src" / "cafe_app"
PAGE = Path(os.environ.get("CAFE_VIZ_PAGE", CAFE_APP / "tools" / "Visualization.py"))  # override to run against another copy
sys.path.insert(0, str(CAFE_APP))  # for `theme`; conftest already adds tools/
matplotlib.use("Agg")  # AppTest runs the page off the main thread; the macOS GUI backend refuses that

import _plot_export  # noqa: E402
import session_store  # noqa: E402
import theme  # noqa: E402
import utils  # noqa: E402


@pytest.fixture(autouse=True)
def _temp_store(monkeypatch, tmp_path):
    for mod in (_plot_export, session_store, theme, utils):
        assert Path(mod.__file__).resolve().is_relative_to(CAFE_APP), mod.__file__
    base = tmp_path / "sessions"
    base.mkdir(mode=0o700)
    monkeypatch.setattr(session_store, "store_base", lambda: base)


def _page(page, open_tabs):
    """Run Visualization.py with only the named lazy tabs open, plotly export stubbed, and per-sample UMAP colours recorded."""
    import runpy

    import _plot_export
    import scanpy as sc
    import streamlit as st
    from matplotlib.colors import to_hex

    real_tabs = st.tabs

    def tabs(labels, **kwargs):
        containers = real_tabs(labels, **kwargs)
        for label, container in zip(labels, containers):
            container.open = label in open_tabs
        return containers

    st.tabs = tabs
    real_export = _plot_export.plotly_image_bytes
    _plot_export.plotly_image_bytes = lambda *a, **k: None  # no Chrome needed

    real_umap = sc.pl.umap

    def umap(adata, *args, **kwargs):
        out = real_umap(adata, *args, **kwargs)
        col = kwargs.get("color")
        if col in ("leiden", "cell_type") and kwargs.get("ax") is not None:
            faces = kwargs["ax"].collections[0].get_facecolors()
            seen = {}
            for cluster, face in zip(adata.obs[col], faces):
                seen.setdefault(str(cluster), set()).add(to_hex(face))
            st.session_state.setdefault(f"_umap_colors_{col}", []).append(seen)
        return out

    sc.pl.umap = umap
    real_download = st.download_button

    def download_button(*args, **kwargs):
        label = kwargs.get("label", args[0] if args else None)
        st.session_state.setdefault("_downloads", []).append((label, kwargs.get("file_name")))
        return real_download(*args, **kwargs)

    st.download_button = download_button
    real_pyplot = st.pyplot

    def pyplot(fig=None, *args, **kwargs):
        import matplotlib.pyplot as plt

        fig = fig or plt.gcf()
        texts = set()
        for ax in getattr(fig, "axes", []):  # a seaborn ClusterGrid has none
            texts.update(t.get_text() for t in [*ax.get_xticklabels(), *ax.get_yticklabels()])
            if ax.get_legend():
                texts.update(t.get_text() for t in ax.get_legend().get_texts())
        st.session_state.setdefault("_pyplot_text", []).append(texts)
        st.session_state.setdefault("_pyplot_seq", []).append([  # same texts, in drawing order
            t.get_text() for ax in getattr(fig, "axes", [])
            for t in [*ax.get_xticklabels(), *ax.get_yticklabels(), *(ax.get_legend().get_texts() if ax.get_legend() else [])]])
        return real_pyplot(fig, *args, **kwargs)

    st.pyplot = pyplot
    try:
        runpy.run_path(page, run_name="__main__")
    finally:
        st.tabs = real_tabs
        sc.pl.umap = real_umap
        st.download_button = real_download
        st.pyplot = real_pyplot
        _plot_export.plotly_image_bytes = real_export


def _adata(n_groups=2, tiny_cluster=False, negative_marker=False, sample_missing_cluster=False):
    rng = np.random.default_rng(0)
    n = 240
    X = rng.normal(loc=500, scale=100, size=(n, 4)).astype("float32")
    if negative_marker:
        X[:, 3] = rng.normal(loc=-5, scale=1, size=n)
    leiden = np.tile(["0", "1", "2"], n // 3)
    samples = np.repeat(["S1", "S2", "S3", "S4"], n // 4)
    if sample_missing_cluster:
        leiden[(samples == "S1") & (leiden == "1")] = "0"
    if tiny_cluster:
        leiden[0] = "9"
    obs = pd.DataFrame(
        {
            "SampleID": samples,
            "Group": np.resize([f"G{i}" for i in range(n_groups)], n),
            "leiden": pd.Categorical(leiden),
            "total_signal": rng.normal(size=n),
        },
        index=[str(i) for i in range(n)],
    )
    adata = sc.AnnData(X, obs=obs)
    adata.var_names = ["CD3", "CD4", "CD8", "NEG"]
    adata.obsm["X_umap"] = rng.normal(size=(n, 2))
    return adata


def _run(adata, open_tabs=()):
    at = AppTest.from_function(_page, args=(str(PAGE), tuple(open_tabs)), default_timeout=120)
    at.session_state["adata"] = adata
    at.run()
    assert not at.exception, at.exception
    return at


def _click(at, label):
    next(b for b in at.button if b.label == label).click().run()
    assert not at.exception, at.exception
    return at


def _warnings(at):
    return " ".join(w.value for w in at.warning)


def test_missing_group_column_shows_friendly_error():
    adata = _adata()
    del adata.obs["Group"]
    at = _run(adata)
    assert any("Group" in e.value for e in at.error)


def test_one_group_sankey_and_bar_by_group_warn_instead_of_raising():
    at = _run(_adata(n_groups=1), ["Sankey: Cluster/Group", "Bar (by Group)"])
    _click(at, "Generate Sankey Diagram")
    assert "at least two groups" in _warnings(at)
    _click(at, "Plot Bar Plot")
    assert "at least two groups" in _warnings(at)


def test_six_groups_cell_percentage_and_count():
    adata = _adata(n_groups=6)
    adata.obs["SampleID"] = "S_" + adata.obs["Group"]  # the tab takes its groups from the sample -> group map
    at = _run(adata, ["Cell % & Count"])
    _click(at, "Plot Cell Percentage and Count")


def test_one_cell_cluster_skips_marker_ranking():
    at = _run(_adata(tiny_cluster=True), ["Dot Plot", "Compare 3 Clusters", "Fold Change"])
    _click(at, "Apply")  # Overview
    assert "fewer than 2 cells" in _warnings(at)
    for label in ("Plot Dot Plot and Dendrogram", "Submit", "Run Analysis"):
        _click(at, label)
        assert "group(s) 9 have fewer than 2 cells" in _warnings(at), label


def test_marker_expression_single_marker_single_column():
    at = _run(_adata(), ["Marker Expression"])
    at.multiselect(key="selected_markers_no_split").set_value(["CD3"])
    at.slider(key="num_of_col_no_split").set_value(1)
    _click(at, "Plot Marker Expression UMAPs")
    assert len(at.get("download_button")) == 1


def test_tab40bc_colormap_plots():
    at = _run(_adata(), ["Clusters"])
    at.selectbox(key="umap_clusters_cmap").set_value("tab40bc")
    _click(at, "Plot UMAPs")
    assert len(at.get("download_button")) >= 3


def test_group_by_options_exclude_numeric_columns_and_histogram_skips_flat_groups():
    adata = _adata()
    adata.obs["SampleID"] = adata.obs["SampleID"].astype(str)
    adata.obs.iloc[0, adata.obs.columns.get_loc("SampleID")] = "Solo"  # one-cell group -> no KDE
    at = _run(adata, ["Histogram", "Violin + Bar"])
    for key in ("kde_groupby_selector", "violin_groupby_selector"):
        options = at.selectbox(key=key).options
        assert "total_signal" not in options and "Group" in options, options
    at.selectbox(key="kde_groupby_selector").set_value("SampleID")
    at.slider(key="kde_num_cols").set_value(1)
    _click(at, "Plot Histogram")


def test_selected_clusters_empty_selection_warns():
    at = _run(_adata(), ["Selected Clusters"])
    at.multiselect(key="selected_clusters_unique").set_value([])
    _click(at, "Plot UMAP")
    assert "at least one cluster" in _warnings(at)


def test_fold_change_keeps_marker_with_negative_mean():
    at = _run(_adata(negative_marker=True), ["Fold Change"])
    _click(at, "Run Analysis")
    fold_chart = next(
        c for c in at.get("plotly_chart") if "Fold Change of Markers" in c.proto.spec
    )
    assert "NEG" in json.loads(fold_chart.proto.spec)["data"][0]["x"]


def test_per_sample_umaps_use_one_cluster_colour_map():
    at = _run(_adata(sample_missing_cluster=True), ["Per Sample"])
    _click(at, "Plot UMAP of each sample")
    per_plot = at.session_state["_umap_colors_leiden"]
    assert len(per_plot) == 4
    colour = {}
    for seen in per_plot:
        for cluster, hexes in seen.items():
            assert len(hexes) == 1, (cluster, hexes)
            colour.setdefault(cluster, set()).update(hexes)
    assert all(len(h) == 1 for h in colour.values()), colour
    assert any("1" not in seen for seen in per_plot)  # a sample really lacks cluster 1


def _adata12(categorical, cell_type=False):
    """12 clusters ('0'..'11', so text order != natural order); sample S1 lacks cluster 11."""
    rng = np.random.default_rng(1)
    n = 240
    leiden = np.tile([str(i) for i in range(12)], n // 12)
    samples = np.repeat(["S1", "S2", "S3", "S4"], n // 4)
    leiden[(samples == "S1") & (leiden == "11")] = "0"
    obs = pd.DataFrame(
        {
            "SampleID": samples,
            "Group": np.resize(["G0", "G1"], n),
            "leiden": pd.Categorical(leiden, categories=[str(i) for i in range(12)]) if categorical else leiden,
            "total_signal": rng.normal(size=n),
        },
        index=[str(i) for i in range(n)],
    )
    if cell_type:
        obs["cell_type"] = "CT" + obs["leiden"].astype(str)
    adata = sc.AnnData(rng.normal(loc=500, scale=100, size=(n, 4)).astype("float32"), obs=obs)
    adata.var_names = ["CD3", "CD4", "CD8", "NEG"]
    adata.obsm["X_umap"] = rng.normal(size=(n, 2))
    return adata


def _one_colour_per_category(plots):
    colour = {}
    for seen in plots:
        for cat, hexes in seen.items():
            colour.setdefault(cat, set()).update(hexes)
    assert all(len(h) == 1 for h in colour.values()), colour
    return {cat: next(iter(h)) for cat, h in colour.items()}


@pytest.mark.parametrize("categorical", [False, True])
def test_twelve_clusters_keep_one_colour_across_cluster_and_sample_umaps(categorical):
    at = _run(_adata12(categorical), ["Clusters", "Per Sample"])
    _click(at, "Plot UMAPs")
    _click(at, "Plot UMAP of each sample")
    plots = at.session_state["_umap_colors_leiden"]
    assert len(plots) == 1 + 2 + 4  # combined, per group, per sample
    assert any("11" not in seen for seen in plots[3:])  # a sample really lacks cluster 11
    colour = _one_colour_per_category(plots)
    assert len(set(colour.values())) == 12


def test_per_sample_default_colours_follow_scanpy_order_without_private_api():
    assert "sc.pl._utils" not in PAGE.read_text()
    adata = _adata12(categorical=False)
    at = _run(adata, ["Per Sample"])
    _click(at, "Plot UMAP of each sample")
    colour = _one_colour_per_category(at.session_state["_umap_colors_leiden"])
    expected = {str(i): matplotlib.colors.to_hex(c) for i, c in enumerate(sc.pl.palettes.default_20[:12])}
    assert colour == expected
    assert "leiden_colors" not in at.session_state["adata"].uns  # session data not modified


def test_clusters_tab_per_group_plots_use_cell_type_colours():
    at = _run(_adata12(categorical=False, cell_type=True), ["Clusters"])
    _click(at, "Plot UMAPs")
    plots = at.session_state["_umap_colors_cell_type"]
    assert len(plots) == 1 + 2  # combined and both groups coloured by cell type
    assert len(set(_one_colour_per_category(plots).values())) == 12


def test_dot_plot_font_size_does_not_leak():
    before = matplotlib.rcParams["font.size"]
    at = _run(_adata(), ["Dot Plot"])
    at.slider(key="font_size").set_value(31)
    _click(at, "Plot Dot Plot and Dendrogram")
    assert matplotlib.rcParams["font.size"] == before


def test_mime_helper():
    import ast

    tree = ast.parse(PAGE.read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_mime")
    ns = {}
    exec(compile(ast.Module([node], []), str(PAGE), "exec"), ns)
    mime = ns["_mime"]
    assert [mime(f) for f in ("PDF", "svg", "PNG", "jpg", "jpeg")] == [
        "application/pdf", "image/svg+xml", "image/png", "image/jpeg", "image/jpeg"]
    assert 'image/{' not in PAGE.read_text()


def test_profile_bar_group_by_excludes_numeric_columns():
    at = _run(_adata(), ["Profile (Bar)"])
    options = at.selectbox(key="bar_groupby_selector_50").options
    assert "total_signal" not in options and "Group" in options, options
    _click(at, "Plot Bar Plot")


def test_per_sample_keeps_existing_categorical_order():
    adata = _adata12(categorical=True)
    order = [str(i) for i in reversed(range(12))]
    adata.obs["leiden"] = adata.obs["leiden"].cat.reorder_categories(order)
    at = _run(adata, ["Per Sample"])
    _click(at, "Plot UMAP of each sample")
    colour = _one_colour_per_category(at.session_state["_umap_colors_leiden"])
    assert colour == {c: matplotlib.colors.to_hex(p) for c, p in zip(order, sc.pl.palettes.default_20)}


@pytest.mark.parametrize("stored_colours,categorical", [(False, False), (False, True), (True, True)])
def test_per_sample_colours_by_cell_type(stored_colours, categorical):
    adata = _adata12(categorical=True, cell_type=True)
    if categorical:
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    cats = sorted(set(adata.obs["cell_type"]))
    if stored_colours:
        adata.uns["cell_type_colors"] = [matplotlib.colors.to_hex(c) for c in sc.pl.palettes.default_102[30:42]]
    at = _run(adata, ["Per Sample"])
    _click(at, "Plot UMAP of each sample")
    assert "_umap_colors_leiden" not in at.session_state
    plots = at.session_state["_umap_colors_cell_type"]
    assert len(plots) == 4 and any("CT11" not in seen for seen in plots)  # S1 lacks CT11
    colour = _one_colour_per_category(plots)
    assert len(set(colour.values())) == 12
    if stored_colours:
        assert colour == dict(zip(cats, adata.uns["cell_type_colors"]))


@pytest.mark.parametrize("cell_type", [False, True])
def test_clusters_tab_download_names_follow_colour_column(cell_type):
    at = _run(_adata12(categorical=True, cell_type=cell_type), ["Clusters"])
    _click(at, "Plot UMAPs")
    label, file_name = at.session_state["_downloads"][0]
    if cell_type:
        assert label == "Download UMAP Cell type (PNG)" and file_name.startswith("UMAP_cell_type_"), (label, file_name)
    else:
        assert label == "Download UMAP Leiden (PNG)" and file_name.startswith("UMAP_leiden_"), (label, file_name)


def test_page_keeps_categorical_cluster_column():
    adata = _adata12(categorical=True)
    cats = list(adata.obs["leiden"].cat.categories)
    at = _run(adata, ["Heatmap"])
    _click(at, "Generate Heatmap")
    col = at.session_state["adata"].obs["leiden"]
    assert isinstance(col.dtype, pd.CategoricalDtype) and list(col.cat.categories) == cats


def test_page_converts_integer_cluster_column_to_strings():
    adata = _adata12(categorical=False)
    adata.obs["leiden"] = adata.obs["leiden"].astype(int)
    at = _run(adata)
    col = at.session_state["adata"].obs["leiden"]
    assert pd.api.types.is_string_dtype(col) and set(col) == {str(i) for i in range(12)}


def test_unused_cluster_category_does_not_block_marker_ranking():
    adata = _adata12(categorical=True)
    adata.obs["leiden"] = adata.obs["leiden"].cat.add_categories(["99"])
    at = _run(adata, ["Dot Plot", "Fold Change"])
    _click(at, "Plot Dot Plot and Dendrogram")
    _click(at, "Run Analysis")
    assert "fewer than 2 cells" not in _warnings(at)


def test_compare_three_clusters_has_own_test_selector():
    at = _run(_adata(), ["Compare 3 Clusters"])
    radio = at.radio(key="test_19")
    assert list(radio.options) == ["Parametric (T-test)", "Non-parametric (Wilcoxon)"]
    assert radio.value == "Parametric (T-test)"
    radio.set_value("Non-parametric (Wilcoxon)")
    _click(at, "Submit")
    assert at.session_state["adata"].uns["rank_genes_groups"]["params"]["method"] == "wilcoxon"


def _one_group_per_sample(categorical):
    adata = _adata(n_groups=3)
    adata.obs["SampleID"] = "S_" + adata.obs["Group"].astype(str)
    if categorical:
        adata.obs["SampleID"] = adata.obs["SampleID"].astype("category")
        adata.obs["Group"] = adata.obs["Group"].astype("category")
    return adata


@pytest.mark.parametrize("categorical", [False, True])
def test_cell_percentage_warns_for_sample_in_two_groups(categorical):
    adata = _one_group_per_sample(categorical)
    adata.obs.iloc[0, adata.obs.columns.get_loc("Group")] = "G1"  # cell 0 is S_G0's: S_G0 now spans G0 and G1
    at = _run(adata, ["Cell % & Count"])
    _click(at, "Plot Cell Percentage and Count")
    text = _warnings(at)
    assert "S_G0" in text and "more than one Group" in text and "S_G1" not in text and "S_G2" not in text
    assert len(at.get("download_button")) == 3


@pytest.mark.parametrize("categorical", [False, True])
def test_cell_percentage_one_group_per_sample_does_not_warn(categorical):
    at = _run(_one_group_per_sample(categorical), ["Cell % & Count"])
    _click(at, "Plot Cell Percentage and Count")
    assert "more than one Group" not in _warnings(at)


@pytest.mark.parametrize("tab,button", [
    ("Per Sample", "Plot UMAP of each sample"),
    ("Cluster Heatmaps", "Plot Combined Heatmaps"),
    ("Composition Bars", "Plot Combined Bar Charts"),
    ("Cell % & Count", "Plot Cell Percentage and Count"),
])
def test_sample_group_tabs_keep_every_group_with_categorical_columns(tab, button):
    at = _run(_one_group_per_sample(categorical=True), [tab])
    _click(at, button)
    groups = {"G0", "G1", "G2"}
    if tab == "Per Sample":  # only uses the samples, never their group
        assert {label for label, _ in at.session_state["_downloads"][:3]} == {
            f"Download UMAP for S_{g} (PNG)" for g in groups}
    elif tab == "Composition Bars":
        specs = [json.loads(c.proto.spec) for c in at.get("plotly_chart") if "per Group" in c.proto.spec]
        assert len(specs) == 2 and all(set(s["data"][0]["x"]) == groups for s in specs), specs
    else:
        texts = at.session_state["_pyplot_text"]
        assert texts and all(groups <= t for t in texts), texts


@pytest.mark.parametrize("categorical,run_overview", [(True, False), (False, True), (False, False)])
def test_cluster_frequency_boxplot_with_categorical_group(categorical, run_overview):
    at = _run(_one_group_per_sample(categorical), ["Cluster Frequency"])
    if run_overview:
        _click(at, "Apply")  # Overview's scanpy calls turn the string Group into a categorical
        assert isinstance(at.session_state["adata"].obs["Group"].dtype, pd.CategoricalDtype)
    for group2 in ("G0", "G2"):  # same group twice (the default), then two of the three groups
        for label, value in (("Select Group 1:", "G0"), ("Select Group 2:", group2)):
            next(s for s in at.selectbox if s.label == label).set_value(value)
        _click(at, "Submit")
    downloads = [label for label, _ in at.session_state["_downloads"]]
    assert sum(label.startswith("Download Cell Count Plot") for label in downloads) == 2, downloads


def test_overview_ranks_despite_unused_leiden_category():
    adata = _adata12(categorical=True, cell_type=True)  # the page's own cleanup only touches cell_type here
    adata.obs["leiden"] = adata.obs["leiden"].cat.add_categories(["99"])
    at = _click(_run(adata), "Apply")  # Overview
    assert "fewer than 2 cells" not in _warnings(at)
    assert at.session_state["adata"].uns["rank_genes_groups"]["params"]["groupby"] == "leiden"


def test_group_correlation_skips_single_cluster_group():
    adata = _adata(n_groups=2)
    adata.obs.loc[adata.obs["Group"] == "G0", "leiden"] = "0"  # G0 now has one cluster, G1 keeps three
    at = _click(_run(adata, ["Group Correlation"]), "Run Correlation and Dendrogram Analysis for Each Group")
    text = _warnings(at)
    assert "Group G0 skipped" in text and "found 1 cluster(s)" in text and "G1" not in text, text
    labels = [label for label, _ in at.session_state["_downloads"]]
    assert "Download Correlation_Group_G1 (PNG)" in labels and not any("Group_G0" in l for l in labels), labels


GROUP_NAMES = ["G10", "G2", "G1", "G21", "G3", "G5"]  # 6 strings: set order matches either expected order by chance ~1/720


@pytest.mark.parametrize("categorical", [False, True])
@pytest.mark.parametrize("tab,button", [
    ("Cluster Heatmaps", "Plot Combined Heatmaps"),
    ("Composition Bars", "Plot Combined Bar Charts"),
    ("Cell % & Count", "Plot Cell Percentage and Count"),
])
def test_group_tabs_use_categorical_or_natural_group_order(tab, button, categorical):
    adata = _adata()
    adata.obs["Group"] = np.resize(GROUP_NAMES, adata.n_obs)
    adata.obs["SampleID"] = "S_" + adata.obs["Group"]
    if categorical:
        custom = ["G21", "G1", "G99", "G10", "G3", "G5", "G2"]  # G99 has no cells
        adata.obs["Group"] = pd.Categorical(adata.obs["Group"], categories=custom)
        expected = [g for g in custom if g != "G99"]
    else:
        expected = ["G1", "G2", "G3", "G5", "G10", "G21"]
    at = _click(_run(adata, [tab]), button)
    if tab == "Composition Bars":
        specs = [json.loads(c.proto.spec) for c in at.get("plotly_chart") if "per Group" in c.proto.spec]
        assert specs and all(s["data"][0]["x"] == expected for s in specs), specs
    else:
        seqs = [[t for t in seq if t in GROUP_NAMES] for seq in at.session_state["_pyplot_seq"]]
        seqs = [s for s in seqs if s]
        assert seqs and all(s == expected for s in seqs), seqs


def test_differential_abundance_ignores_unused_group_category():
    adata = _adata()
    adata.obs["Group"] = pd.Categorical(
        adata.obs["SampleID"].map({"S1": "G0", "S2": "G0", "S3": "G1", "S4": "G1"}), categories=["G0", "G1", "G9"])
    at = _click(_run(adata, ["Differential Abundance"]), "Run Differential Abundance Test")
    assert "G9" not in _warnings(at) and "G9" not in " ".join(m.value for m in at.markdown)
    assert len(at.dataframe) >= 1


@pytest.mark.parametrize("tab,button", [
    ("Cluster Heatmaps", "Plot Combined Heatmaps"),
    ("Composition Bars", "Plot Combined Bar Charts"),
    ("Cell % & Count", "Plot Cell Percentage and Count"),
])
def test_group_tabs_categorical_group_with_unused_category_no_pandas4_warning(tab, button):
    import warnings

    adata = _one_group_per_sample(categorical=True)
    adata.obs["Group"] = adata.obs["Group"].cat.add_categories(["G9"])  # G9 has no cells
    at = _run(adata, [tab])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _click(at, button)
    assert not [w for w in caught if issubclass(w.category, pd.errors.Pandas4Warning)], [str(w.message) for w in caught]


def test_group_correlation_heading_only_for_plotted_groups():
    adata = _adata(n_groups=2)
    adata.obs.loc[adata.obs["Group"] == "G0", "leiden"] = "0"  # G0 is skipped, G1 is plotted
    at = _click(_run(adata, ["Group Correlation"]), "Run Correlation and Dendrogram Analysis for Each Group")
    headings = [m.value for m in at.markdown if m.value.startswith("### Correlation Plot for Group")]
    assert headings == ["### Correlation Plot for Group: G1"], headings
    assert "Group G0 skipped" in _warnings(at)
