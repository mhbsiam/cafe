"""Per-session disk persistence for the loaded AnnData in local-only CAFE.

A hard browser refresh clears both ``st.session_state`` and the ``st.file_uploader`` widget, so the
in-memory ``adata`` is lost and the user must re-upload. This module keeps each session's AnnData on disk,
keyed by a token carried in the URL query param (``?s=<token>``), which is the one piece of state that survives a
refresh. On the next run we reload from disk, so data reappears without a re-upload.

Data is stored under ``~/.cafe/sessions`` and survives refresh AND app restart, overwritten on each new
load. The URL token is a local recovery key, not an authentication mechanism; this module assumes the
Streamlit server is bound to localhost.
"""
import json
import os
import re
import shutil
import tempfile
import time
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

import anndata as ad
import streamlit as st

from utils import validate_path

# uuid4().hex is 32 lowercase hex chars; reject anything else so a tampered ?s= can't escape the store base.
_TOKEN_RE = re.compile(r"^[0-9a-f]{32}$")
_ADATA_NAME = "adata.h5ad"
_META_NAME = "meta.json"
_DEFAULT_RETENTION_DAYS = 7
_PRUNE_INTERVAL_SECONDS = 60 * 60
_last_prune_at = None


def _owner_only(path, mode):
    """Best-effort POSIX permissions; harmless on platforms that do not support chmod."""
    try:
        path.chmod(mode)
    except OSError:
        pass


def store_base():
    """Root directory holding all per-session stores."""
    cafe_dir = Path.home() / ".cafe"
    cafe_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    _owner_only(cafe_dir, 0o700)
    base = cafe_dir / "sessions"
    base.mkdir(parents=True, exist_ok=True, mode=0o700)
    _owner_only(base, 0o700)
    return base


def token_dir(token):
    """Directory for ``token``, created lazily. Raises on a malformed/traversing token."""
    if not token or not _TOKEN_RE.match(token):
        raise ValueError("Invalid session token.")
    base = store_base()
    # Defense in depth on top of the regex: guarantee the resolved path stays inside the base.
    path = Path(validate_path(str(base / token), base=str(base)))
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    _owner_only(path, 0o700)
    return path


def _retention_days():
    try:
        return max(1, int(os.environ.get("CAFE_SESSION_RETENTION_DAYS", _DEFAULT_RETENTION_DAYS)))
    except ValueError:
        return _DEFAULT_RETENTION_DAYS


def _session_timestamp(session_dir):
    """Return the saved timestamp, falling back to the directory modification time."""
    meta_path = session_dir / _META_NAME
    try:
        saved_at = json.loads(meta_path.read_text()).get("saved_at")
        if saved_at:
            return datetime.fromisoformat(saved_at.replace("Z", "+00:00")).astimezone(timezone.utc)
    except (OSError, ValueError, TypeError):
        pass
    return datetime.fromtimestamp(session_dir.stat().st_mtime, tz=timezone.utc)


def prune_expired_sessions(*, max_age_days=None, now=None):
    """Delete abandoned token directories older than the configured retention period."""
    base = store_base()
    cutoff = (now or datetime.now(timezone.utc)) - timedelta(
        days=max_age_days if max_age_days is not None else _retention_days()
    )
    removed = 0
    for session_dir in base.iterdir():
        if not session_dir.is_dir() or not _TOKEN_RE.fullmatch(session_dir.name):
            continue
        try:
            if _session_timestamp(session_dir) < cutoff:
                shutil.rmtree(session_dir)
                removed += 1
        except OSError:
            continue
    return removed


def _maybe_prune_expired_sessions():
    """Limit cleanup scans to once per process per hour despite Streamlit reruns."""
    global _last_prune_at
    now = time.monotonic()
    if _last_prune_at is not None and now - _last_prune_at < _PRUNE_INTERVAL_SECONDS:
        return
    _last_prune_at = now
    prune_expired_sessions()


def ensure_token():
    """Ensure a stable session token exists in the URL and session_state; return it.

    The URL query param is what survives a refresh; session_state mirrors it for convenient access.
    """
    _maybe_prune_expired_sessions()
    token = st.query_params.get("s")
    if not token or not _TOKEN_RE.match(token):
        token = uuid.uuid4().hex
        st.query_params["s"] = token
    st.session_state["_cafe_token"] = token
    return token


def current_token():
    """The active token (mirrored into session_state by ensure_token)."""
    return st.session_state.get("_cafe_token")


def save_adata(token, adata, *, source, filename=None):
    """Atomically persist ``adata`` for ``token`` and record lightweight metadata."""
    d = token_dir(token)
    target = d / _ADATA_NAME
    # Write to a temp file in the same dir, then os.replace for an atomic swap (no half-written .h5ad).
    fd, tmp = tempfile.mkstemp(suffix=".h5ad", dir=str(d))
    os.close(fd)
    try:
        adata.write(tmp)
        os.replace(tmp, target)
        _owner_only(target, 0o600)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    meta = {
        "filename": filename,
        "source": source,
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
    }
    meta_path = d / _META_NAME
    meta_path.write_text(json.dumps(meta))
    _owner_only(meta_path, 0o600)


def load_adata(token):
    """Return the persisted AnnData for ``token``, or None if nothing is stored."""
    if not token or not _TOKEN_RE.match(token):
        return None
    target = store_base() / token / _ADATA_NAME
    if not target.exists():
        return None
    return ad.read_h5ad(target)


def has_adata(token):
    if not token or not _TOKEN_RE.match(token):
        return False
    return (store_base() / token / _ADATA_NAME).exists()


def load_meta(token):
    if not token or not _TOKEN_RE.match(token):
        return {}
    meta = store_base() / token / _META_NAME
    if not meta.exists():
        return {}
    try:
        return json.loads(meta.read_text())
    except (ValueError, OSError):
        return {}


def clear(token):
    """Remove the entire on-disk store for ``token`` (no-op if absent)."""
    if not token or not _TOKEN_RE.match(token):
        return
    shutil.rmtree(store_base() / token, ignore_errors=True)


def clear_current_data(*, extra_keys=()):
    """Delete the active disk store and reset dataset-related Streamlit state."""
    clear(current_token())
    for key in (*_CLEARABLE_KEYS, *extra_keys):
        st.session_state.pop(key, None)


def set_adata(adata, *, source, filename=None):
    """Set ``st.session_state.adata`` and persist it under the current session token.

    Use this in place of a bare ``st.session_state.adata = adata`` at every write site so the object
    survives a refresh. Persistence failures are surfaced as a warning but never block the in-memory update.
    """
    st.session_state["adata"] = adata
    token = current_token() or ensure_token()
    try:
        save_adata(token, adata, source=source, filename=filename)
    except Exception as exc:  # persistence is best-effort; the live session still has the object
        st.warning(f"Could not persist data for refresh recovery: {exc}")


# Session-state keys wiped alongside the on-disk store when the user clears data. Keeps a "Clear" reset from
# leaving stale pipeline flags that would desync the next upload.
_CLEARABLE_KEYS = (
    "adata", "_dp_started", "_viz_adata", "_viz_adata_key", "_viz_pending_upload",
    "qc_done", "pca_done", "pca_selected", "perform_pca",
    "batch_correction_done", "leiden_computed", "umap_computed", "uploaded_files",
)


def render_clear_control():
    """Sidebar button that clears the persisted store and resets in-memory data. Call after pg.run()."""
    token = current_token()
    if not token or not has_adata(token):
        return
    with st.sidebar:
        if st.button("🗑️ Clear loaded data", help="Remove the saved dataset from this session."):
            clear_current_data()
            st.rerun()


def restore_session_adata():
    """Restore ``adata`` from disk once per session, if the session has no live copy.

    Runs at most once per session (guarded by ``_cafe_restored``), so a deliberate in-app reset of
    ``adata`` to ``None`` is not resurrected on the next rerun. A hard refresh creates a fresh session
    (guard absent), which is exactly when we do want to reload from disk.
    """
    if st.session_state.get("_cafe_restored"):
        return
    st.session_state["_cafe_restored"] = True
    if st.session_state.get("adata") is not None:
        return
    token = current_token()
    if token and has_adata(token):
        try:
            st.session_state["adata"] = load_adata(token)
        except Exception as exc:
            st.warning(f"Could not restore saved data: {exc}")
