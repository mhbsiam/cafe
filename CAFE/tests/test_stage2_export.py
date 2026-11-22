"""Tests for the plotly image export helpers in _plot_export.py."""
import os

import plotly.graph_objects as go
import pytest
from streamlit.testing.v1 import AppTest

import _plot_export
from _plot_export import plotly_image_bytes

TOOLS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "cafe_app", "tools")


def _small_fig():
    return go.Figure(go.Bar(x=["a", "b"], y=[1, 2]))


def test_module_loaded_from_this_checkout():
    assert os.path.realpath(_plot_export.__file__).startswith(os.path.realpath(TOOLS_DIR))


def test_plotly_image_bytes_returns_png_bytes():
    data = plotly_image_bytes(_small_fig(), "png")
    if data is None:
        pytest.skip("Chrome/Chromium not installed")
    assert data[:8] == b"\x89PNG\r\n\x1a\n"


def _chrome_missing_script():
    import plotly.graph_objects as go
    import streamlit as st
    from kaleido.errors import ChromeNotFoundError
    from _plot_export import plotly_image_bytes

    def broken_write_image(*args, **kwargs):
        raise ChromeNotFoundError("no chrome")

    results = []
    for _ in range(2):
        fig = go.Figure()
        fig.write_image = broken_write_image
        results.append(plotly_image_bytes(fig, "png"))
    st.write(repr(results))


def _plotly_runtime_error_script():
    import plotly.graph_objects as go
    import streamlit as st
    from plotly.io._kaleido import PLOTLY_GET_CHROME_ERROR_MSG
    from _plot_export import plotly_image_bytes

    def broken_write_image(*args, **kwargs):
        raise RuntimeError(PLOTLY_GET_CHROME_ERROR_MSG)

    fig = go.Figure()
    fig.write_image = broken_write_image
    st.write(repr([plotly_image_bytes(fig, "png"), plotly_image_bytes(fig, "pdf")]))


@pytest.mark.parametrize("script", [_chrome_missing_script, _plotly_runtime_error_script])
def test_missing_chrome_returns_none_and_warns_once_per_run(script):
    at = AppTest.from_function(script)
    at.run()
    assert not at.exception
    assert at.markdown[0].value == "[None, None]"
    assert len(at.warning) == 1
    assert "Chrome" in at.warning[0].value
    # A rerun is a new run, so the warning is shown again (not suppressed by a stale flag).
    at.run()
    assert not at.exception
    assert len(at.warning) == 1


def test_real_kaleido_chrome_error_is_recognised():
    kaleido_errors = pytest.importorskip("kaleido.errors")
    assert _plot_export._is_chrome_missing(kaleido_errors.ChromeNotFoundError("x"))


def test_unrelated_exception_propagates(monkeypatch):
    fig = _small_fig()

    def broken_write_image(*args, **kwargs):
        raise ValueError("EPS export is not supported")

    monkeypatch.setattr(fig, "write_image", broken_write_image)
    with pytest.raises(ValueError):
        plotly_image_bytes(fig, "eps")
    monkeypatch.setattr(fig, "write_image", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("other")))
    with pytest.raises(RuntimeError):
        plotly_image_bytes(fig, "png")


def _form_script():
    import plotly.graph_objects as go
    import streamlit as st
    from _plot_export import deferred_download_button

    with st.form("f"):
        submitted = st.form_submit_button("Run")
    if submitted:
        fig = go.Figure(go.Bar(x=["a", "b"], y=[1, 2]))
        # Stub the renderer so the test does not need Chrome; the widget flow is what is under test.
        fig.write_image = lambda buf, format=None, **kw: buf.write(b"%PDF-stub")
        deferred_download_button(fig, "pdf", "x.pdf", "application/pdf", key_prefix="t")


def test_download_button_appears_in_submit_run():
    at = AppTest.from_function(_form_script)
    at.run()
    assert len(at.get("download_button")) == 0
    at.button[0].click().run()
    assert not at.exception
    assert len(at.button) == 1  # only the form submit button: no "Prepare" step
    buttons = at.get("download_button")
    assert len(buttons) == 1
    assert buttons[0].proto.label == "Download PDF"


def test_real_plotly_path_without_chrome_returns_none(monkeypatch):
    from choreographer.browsers.chromium import Chromium

    monkeypatch.setattr(Chromium, "find_browser", classmethod(lambda cls, **kw: None))
    assert plotly_image_bytes(_small_fig(), "png") is None
