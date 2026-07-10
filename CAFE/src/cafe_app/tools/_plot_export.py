"""Deferred-download helper for Plotly exports: a Prepare button runs write_image() once and stashes bytes, then a download button appears."""
import io
import streamlit as st


def deferred_download_button(fig, file_format, file_name, mime_type, key_prefix,
                             *, write_fn=None, label=None):
    """Two-step 'Prepare download' -> 'Download' button; defers write_image() until requested. key_prefix must be page-unique."""
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
