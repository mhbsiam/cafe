"""Deferred-download helper for Plotly figure exports.

Streamlit's ``st.download_button`` needs its ``data=`` ready at render time, so
the expensive ``write_image()`` (kaleido) call normally runs eagerly on every
form submit.  This helper splits the interaction into two steps:

1. A "Prepare <FMT> download" button that, only when clicked, runs
   ``write_image()`` and stashes the bytes in ``st.session_state``.
2. A real ``st.download_button`` that appears once the bytes are prepared.

Because the bytes live in ``st.session_state``, they survive across reruns
(including fragment-scoped reruns), so the user only pays the kaleido cost
once, when they actually ask for the file.
"""
import io
import streamlit as st


def deferred_download_button(fig, file_format, file_name, mime_type, key_prefix,
                             *, write_fn=None, label=None):
    """Two-step 'Prepare download' -> 'Download' button.

    Defers the expensive ``write_image()`` (kaleido) call until the user
    actually asks for the file.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The figure to export.
    file_format : str
        e.g. ``"png"``, ``"pdf"``, ``"svg"``.
    file_name : str
        Download filename.
    mime_type : str
        MIME type for the download button.
    key_prefix : str
        Unique prefix for widget keys (must be unique across the page).
    write_fn : callable, optional
        Custom write function ``fn(fig, buffer, format)``.  Defaults to
        ``fig.write_image(buffer, format=format)``.
    label : str, optional
        Custom label for the prepare button.
    """
    if write_fn is None:
        write_fn = lambda f, buf, fmt: f.write_image(buf, format=fmt)
    buffer_key = f"_export_buffer_{key_prefix}"
    if st.button(label or f"Prepare {file_format.upper()} download",
                 key=f"_export_prepare_{key_prefix}"):
        buf = io.BytesIO()
        write_fn(fig, buf, file_format)
        buf.seek(0)
        st.session_state[buffer_key] = buf.getvalue()
    if buffer_key in st.session_state:
        st.download_button(
            label=f"Download {file_format.upper()}",
            data=st.session_state[buffer_key],
            file_name=file_name,
            mime=mime_type,
            key=f"_export_download_{key_prefix}",
        )
