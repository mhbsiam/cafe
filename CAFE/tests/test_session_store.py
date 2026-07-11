"""Regression tests for local AnnData session persistence."""
import json
import stat
import uuid
from datetime import datetime, timedelta, timezone

import session_store


def _use_temp_store(monkeypatch, tmp_path):
    base = tmp_path / "sessions"
    base.mkdir(mode=0o700)
    monkeypatch.setattr(session_store, "store_base", lambda: base)
    return base


def test_clear_current_data_removes_disk_and_state(monkeypatch, tmp_path):
    base = _use_temp_store(monkeypatch, tmp_path)
    token = uuid.uuid4().hex
    session_dir = base / token
    session_dir.mkdir()
    (session_dir / "adata.h5ad").write_bytes(b"data")
    state = {"_cafe_token": token, "adata": object(), "qc_done": True, "custom": True}
    monkeypatch.setattr(session_store.st, "session_state", state)

    session_store.clear_current_data(extra_keys=("custom",))

    assert not session_dir.exists()
    assert "adata" not in state
    assert "qc_done" not in state
    assert "custom" not in state
    assert state["_cafe_token"] == token


def test_prune_expired_sessions_keeps_active_and_ignores_malformed(monkeypatch, tmp_path):
    base = _use_temp_store(monkeypatch, tmp_path)
    now = datetime.now(timezone.utc)
    expired = base / uuid.uuid4().hex
    active = base / uuid.uuid4().hex
    malformed = base / "not-a-token"
    for directory, saved_at in (
        (expired, now - timedelta(days=8)),
        (active, now - timedelta(days=1)),
        (malformed, now - timedelta(days=30)),
    ):
        directory.mkdir()
        (directory / "meta.json").write_text(json.dumps({"saved_at": saved_at.isoformat()}))

    assert session_store.prune_expired_sessions(max_age_days=7, now=now) == 1
    assert not expired.exists()
    assert active.exists()
    assert malformed.exists()


def test_prune_uses_directory_mtime_when_metadata_is_invalid(monkeypatch, tmp_path):
    base = _use_temp_store(monkeypatch, tmp_path)
    now = datetime.now(timezone.utc)
    expired = base / uuid.uuid4().hex
    expired.mkdir()
    (expired / "meta.json").write_text("not json")
    old = (now - timedelta(days=8)).timestamp()
    expired.touch()
    import os
    os.utime(expired, (old, old))

    assert session_store.prune_expired_sessions(max_age_days=7, now=now) == 1
    assert not expired.exists()


def test_cleanup_runs_immediately_then_is_throttled(monkeypatch):
    calls = []
    monkeypatch.setattr(session_store, "_last_prune_at", None)
    monkeypatch.setattr(session_store, "prune_expired_sessions", lambda: calls.append(True))
    monkeypatch.setattr(session_store.time, "monotonic", lambda: 10.0)

    session_store._maybe_prune_expired_sessions()
    session_store._maybe_prune_expired_sessions()

    assert calls == [True]


def test_token_directory_and_saved_files_are_owner_only(monkeypatch, tmp_path, small_dense_adata):
    base = _use_temp_store(monkeypatch, tmp_path)
    token = uuid.uuid4().hex

    session_store.save_adata(token, small_dense_adata, source="test")

    session_dir = base / token
    assert stat.S_IMODE(session_dir.stat().st_mode) == 0o700
    assert stat.S_IMODE((session_dir / "adata.h5ad").stat().st_mode) == 0o600
    assert stat.S_IMODE((session_dir / "meta.json").stat().st_mode) == 0o600
