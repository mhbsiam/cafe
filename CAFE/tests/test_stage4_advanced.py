"""Stage 4 fixes: resolution parsing, blank signature names, PCA-by-group plot."""
import ast
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

ADVANCED = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'cafe_app', 'advanced')


def _load_functions(page, names):
    """Exec only the named top-level functions of a Streamlit page (the page body needs a running app)."""
    path = os.path.join(ADVANCED, page)
    tree = ast.parse(open(path).read())
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    ns = {'__file__': path, 'pd': pd}
    exec(compile(tree, path, 'exec'), ns)
    return ns


def test_parse_resolutions():
    parse = _load_functions('Cluster_Evaluation.py', {'parse_resolutions'})['parse_resolutions']
    assert parse("0.1, 0.5,1,") == [0.1, 0.5, 1.0]
    with pytest.raises(ValueError, match="'0' must be a positive"):
        parse("0,5;0.7")  # decimal-comma typo splits into '0' and '5;0.7'
    with pytest.raises(ValueError, match="'5;0.7' is not a number"):
        parse("0.1,5;0.7")
    with pytest.raises(ValueError, match="'abc'"):
        parse("0.5, abc")
    with pytest.raises(ValueError, match="at least one"):
        parse("")
    with pytest.raises(ValueError, match="at least one"):
        parse(" , ")
    with pytest.raises(ValueError, match="'-0.5' must be a positive"):
        parse("0.5,-0.5")
    with pytest.raises(ValueError, match="positive"):
        parse("0")
    with pytest.raises(ValueError, match="positive"):
        parse("nan")


@pytest.mark.parametrize("blank", [None, np.nan, pd.NA, "", "  "])
def test_signature_blank_name_is_never_a_cell_type(blank):
    build = _load_functions('Cluster_Annotation.py', {'_signature_definitions'})['_signature_definitions']
    edited = pd.DataFrame({
        "Cell type": [" T cell ", blank, blank, "B cell"],
        "CD3": ["+", "+", None, ""],
        "CD19": [np.nan, "", None, "+"],
    }, dtype=object)
    definitions, unnamed, duplicates = build(edited, ["CD3", "CD19"])
    assert definitions == {"T cell": {"CD3": "+"}, "B cell": {"CD19": "+"}}
    assert unnamed == [2]  # row 3 (no name, no markers) is skipped silently
    assert duplicates == []


def test_signature_duplicate_names_are_reported_not_overwritten():
    build = _load_functions('Cluster_Annotation.py', {'_signature_definitions'})['_signature_definitions']
    edited = pd.DataFrame({
        "Cell type": ["T cells", "B cell", "t cells "],
        "CD3": ["+", "", "-"],
        "CD19": ["", "+", "+"],
    }, dtype=object)
    definitions, unnamed, duplicates = build(edited, ["CD3", "CD19"])
    assert duplicates == ["T cells"]
    assert unnamed == []
    assert definitions["T cells"] == {"CD3": "+"}  # the later row does not silently win


def test_selective_clustering_pca_by_group_renders_a_figure():
    import anndata as ad
    import scanpy as sc

    src = open(os.path.join(ADVANCED, 'Selective_Clustering.py')).read()
    call = "sc.pl.pca(adata, color='Group', show=False, return_fig=True)"
    assert call in src and "show=True" not in src

    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.normal(size=(40, 5)).astype(np.float32))
    adata.obs['Group'] = pd.Categorical(['a', 'b'] * 20)
    sc.pp.pca(adata, n_comps=3)
    plt.close('all')
    fig = eval(call)
    assert isinstance(fig, matplotlib.figure.Figure) and fig.axes
    plt.close(fig)
    assert not plt.get_fignums()
