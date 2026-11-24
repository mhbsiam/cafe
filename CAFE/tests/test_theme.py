"""Accessibility markup of the theme helpers (info cards and the wizard stepper)."""
from cafe_app import theme


def _capture(monkeypatch):
    out = []
    monkeypatch.setattr(theme.st, "markdown", lambda html, **_: out.append(html))
    return out


def test_info_card_kind_sets_style_role_and_label(monkeypatch):
    out = _capture(monkeypatch)
    theme.info_card("Duplicate samples", "body", kind="error")
    theme.info_card("Tip", "body", kind="info")
    theme.info_card("Odd", "body", kind="nonsense")
    error, info, fallback = out
    assert 'cafe-info-card--error" role="alert"' in error and ">Error</span>" in error
    assert 'cafe-info-card--info" role="status"' in info and "cafe-info-card__kind" not in info
    assert "cafe-info-card--info" in fallback


def test_stepper_marks_current_and_done_steps(monkeypatch):
    out = _capture(monkeypatch)
    theme.stepper(["Upload", "QC", "Cluster"], 1, {0})
    (html,) = out
    assert html.startswith('<ol class="cafe-stepper" aria-label="Progress">')
    assert "Upload<span class=\"cafe-sr-only\"> (completed)</span>" in html
    assert html.count('aria-current="step"') == 1
    assert 'aria-current="step"><span class="cafe-stepper__dot" aria-hidden="true">2</span>' in html
