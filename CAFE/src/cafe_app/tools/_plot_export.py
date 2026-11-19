"""Plotly image export helpers: write_image() needs Chrome (kaleido); without it they warn once per run and skip the export instead of raising."""
import io
import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

_WARNED_RUN_KEY = "_plot_export_chrome_warned_run"


def _is_chrome_missing(exc):
    # kaleido raises ChromeNotFoundError; plotly re-raises it (chained) as RuntimeError(PLOTLY_GET_CHROME_ERROR_MSG)
    for err in (exc, exc.__cause__, exc.__context__):
        if err is not None and any(cls.__name__ == "ChromeNotFoundError" for cls in type(err).__mro__):
            return True
    return isinstance(exc, RuntimeError) and "requires Google Chrome" in str(exc)


def _warn_chrome_missing_once():
    ctx = get_script_run_ctx()
    # ctx.cursors is a fresh dict on every (fragment) run, so holding it marks this run only
    run_marker = getattr(ctx, "cursors", None)
    if run_marker is not None and st.session_state.get(_WARNED_RUN_KEY) is run_marker:
        return
    if run_marker is not None:
        st.session_state[_WARNED_RUN_KEY] = run_marker
    st.warning(
        "Plot image export (PNG/PDF/SVG download) needs Google Chrome or Chromium, which was not found. "
        "Install Google Chrome normally, or run `plotly_get_chrome` in the CAFE environment, then run the plot again. "
        "The interactive plots and other downloads still work."
    )


def plotly_image_bytes(fig, file_format, **write_kwargs):
    """Return fig rendered as file_format bytes, or None (with one warning per run) if Chrome is missing."""
    buf = io.BytesIO()
    try:
        fig.write_image(buf, format=file_format, **write_kwargs)
    except Exception as exc:
        if not _is_chrome_missing(exc):
            raise
        _warn_chrome_missing_once()
        return None
    return buf.getvalue()


def deferred_download_button(fig, file_format, file_name, mime_type, key_prefix, *, label=None):
    """Render fig and show a download button that does not rerun the app; nothing if export is unavailable. key_prefix must be page-unique."""
    data = plotly_image_bytes(fig, file_format)
    if data is not None:
        st.download_button(
            label=label or f"Download {file_format.upper()}",
            data=data,
            file_name=file_name,
            mime=mime_type,
            key=f"_export_download_{key_prefix}",
            on_click="ignore",
        )
