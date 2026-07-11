"""Shared design tokens, global CSS, and UI helpers for CAFE (teal accent, Source Sans 3, reduced-motion aware)."""
import os
from dataclasses import dataclass

import streamlit as st

# Bundled image assets, resolved relative to this file (cwd-independent).
IMG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "img")


@dataclass(frozen=True)
class Tokens:
    """Single source of truth for CAFE design tokens."""

    # Color
    primary: str = "#0f5070"
    primary_hover: str = "#0a3d57"
    primary_light: str = "#e6f1f4"
    primary_lighter: str = "#f0f6f8"
    background: str = "#ffffff"
    surface: str = "#f8f9fa"
    surface_hover: str = "#e9ecef"
    border: str = "#d5e0e4"
    text: str = "#1a262a"
    text_muted: str = "#5c6e74"
    text_subtle: str = "#8fa3aa"
    success: str = "#2d6a4f"
    warning: str = "#b35900"
    error: str = "#c1121f"
    white: str = "#ffffff"

    # Typography
    font: str = '"Source Sans 3", "Source Sans Pro", sans-serif'

    # Spacing (base 4px grid)
    space_xs: str = "4px"
    space_sm: str = "8px"
    space_md: str = "16px"
    space_lg: str = "24px"
    space_xl: str = "32px"
    space_2xl: str = "48px"
    space_3xl: str = "64px"

    # Motion
    expo: str = "cubic-bezier(0.16, 1, 0.3, 1)"
    duration_micro: str = "0.18s"
    duration_fade: str = "0.35s"
    duration_tab: str = "0.28s"
    duration_hero: str = "0.6s"
    duration_progress: str = "0.45s"
    duration_sheen: str = "1.6s"


TOKENS = Tokens()

_EXPO = TOKENS.expo

# CSS variables + font import
_CSS_VARIABLES = f"""
<style>
@import url('https://fonts.bunny.net/css?family=source-sans-3:400,600,700&display=swap');

:root, [data-testid="stApp"] {{
  --cafe-primary: {TOKENS.primary};
  --cafe-primary-hover: {TOKENS.primary_hover};
  --cafe-primary-light: {TOKENS.primary_light};
  --cafe-primary-lighter: {TOKENS.primary_lighter};
  --cafe-background: {TOKENS.background};
  --cafe-surface: {TOKENS.surface};
  --cafe-surface-hover: {TOKENS.surface_hover};
  --cafe-border: {TOKENS.border};
  --cafe-text: {TOKENS.text};
  --cafe-text-muted: {TOKENS.text_muted};
  --cafe-text-subtle: {TOKENS.text_subtle};
  --cafe-success: {TOKENS.success};
  --cafe-warning: {TOKENS.warning};
  --cafe-error: {TOKENS.error};
  --cafe-white: {TOKENS.white};
  --cafe-font: {TOKENS.font};
  --cafe-space-xs: {TOKENS.space_xs};
  --cafe-space-sm: {TOKENS.space_sm};
  --cafe-space-md: {TOKENS.space_md};
  --cafe-space-lg: {TOKENS.space_lg};
  --cafe-space-xl: {TOKENS.space_xl};
  --cafe-space-2xl: {TOKENS.space_2xl};
  --cafe-space-3xl: {TOKENS.space_3xl};
  --cafe-expo: {_EXPO};
  --cafe-duration-micro: {TOKENS.duration_micro};
  --cafe-duration-fade: {TOKENS.duration_fade};
  --cafe-duration-tab: {TOKENS.duration_tab};
  --cafe-duration-hero: {TOKENS.duration_hero};
  --cafe-duration-progress: {TOKENS.duration_progress};
  --cafe-duration-sheen: {TOKENS.duration_sheen};
}}
</style>
"""

# Global typography override
_BASE_CSS = f"""
<style>
html, body, [data-testid="stApp"] {{
  font-family: var(--cafe-font) !important;
}}

[data-testid="stMarkdownContainer"] p,
[data-testid="stMarkdownContainer"] li,
[data-testid="stMarkdownContainer"] td,
[data-testid="stMarkdownContainer"] th,
[data-testid="stTextInput"] label,
[data-testid="stSelectbox"] label,
[data-testid="stSlider"] label,
[data-testid="stRadio"] label,
[data-testid="stCheckbox"] label,
[data-testid="stNumberInput"] label,
[data-testid="stDateInput"] label,
[data-testid="stTimeInput"] label,
[data-testid="stFileUploader"] label,
[data-testid="stTextArea"] label,
[data-testid="stColorPicker"] label,
[data-testid="stDataFrame"] td,
[data-testid="stDataFrame"] th,
[data-testid="stTable"] td,
[data-testid="stTable"] th,
[data-baseweb="select"] div,
[data-baseweb="radio"] span,
[data-baseweb="checkbox"] span,
[data-baseweb="tab"] span,
[data-baseweb="progress-bar"] div,
button[data-testid*="stBaseButton"] div {{
  font-family: var(--cafe-font) !important;
}}

[data-testid="stMain"] .block-container {{
  padding-top: var(--cafe-space-xl) !important;
  padding-bottom: var(--cafe-space-2xl) !important;
}}

[data-testid="stSidebar"] .block-container {{
  padding-top: var(--cafe-space-lg) !important;
}}
</style>
"""

# Premium micro-interactions
_PREMIUM_CSS = f"""
<style>
/* Lock the app to its configured light theme: every chart is a static
   matplotlib/Plotly image with a light background, so letting users switch
   to Streamlit's dark theme leaves charts as bright boxes on a dark page.
   Streamlit 1.4x exposes theme switching as a selectbox + "Edit active
   theme" button inside the Settings dialog; newer Streamlit versions use a
   dedicated stThemeSwitcher control instead. Hide both. */
[data-testid="stThemeSwitcher"] {{ display: none !important; }}
[role="dialog"] label:has(+ [data-testid="stSelectbox"]) {{ display: none !important; }}
[role="dialog"] [data-testid="stSelectbox"]:has(+ [data-testid="edit-theme"]) {{ display: none !important; }}
[data-testid="edit-theme"] {{ display: none !important; }}

/* All motion lives behind this query, so reduced-motion users get the app with
   no animation at all. This is the accessible default, so no separate override is needed. */
@media (prefers-reduced-motion: no-preference) {{

  /* --- Buttons: lift on hover, press down on click ------------------------ */
  [data-testid="stMain"] [data-testid*="stBaseButton-secondary"],
  [data-testid="stMain"] [data-testid*="stBaseButton-primary"] {{
    transition: transform var(--cafe-duration-micro) var(--cafe-expo),
                box-shadow var(--cafe-duration-micro) var(--cafe-expo),
                border-color var(--cafe-duration-micro) ease-out,
                background-color var(--cafe-duration-micro) ease-out;
  }}
  [data-testid="stMain"] [data-testid*="stBaseButton-secondary"]:hover,
  [data-testid="stMain"] [data-testid*="stBaseButton-primary"]:hover {{
    transform: translateY(-1px);
  }}
  [data-testid="stMain"] [data-testid*="stBaseButton-primary"]:hover {{
    /* Mirrors .streamlit/config.toml's primaryColor (#0f5070) */
    box-shadow: 0 4px 14px rgba(15, 80, 112, 0.35);
  }}
  [data-testid="stMain"] [data-testid*="stBaseButton-secondary"]:active,
  [data-testid="stMain"] [data-testid*="stBaseButton-primary"]:active {{
    transform: translateY(0) scale(0.98);
  }}

  /* --- Alerts / status: rise in instead of snapping into place ------------- */
  [data-testid="stMain"] [data-testid="stAlert"] {{
    animation: cafe-rise-in 0.32s var(--cafe-expo);
  }}

  /* --- Progress bar: ease the fill, and run a quiet sheen across it ------- */
  /* Heavy steps (neighbors / UMAP / Leiden) can sit at the same percentage
     for a while. A slow light sheen travelling over the fill signals "still
     working" without the jitter of a spinner. It provides reassurance, not decoration. */
  [data-testid="stMain"] [data-baseweb="progress-bar"] div div {{
    transition: width var(--cafe-duration-progress) var(--cafe-expo);
    background-image: linear-gradient(
      100deg,
      rgba(255, 255, 255, 0) 30%,
      rgba(255, 255, 255, 0.55) 50%,
      rgba(255, 255, 255, 0) 70%
    );
    background-size: 220% 100%;
    background-repeat: no-repeat;
    animation: cafe-sheen var(--cafe-duration-sheen) linear infinite;
  }}

  /* --- Results settle in: charts and tables fade like the plots already do  */
  /* One consistent "the result just arrived" cue across every output type,
     so a Plotly figure or a dataframe doesn't snap in while images glide. */
  [data-testid="stMain"] [data-testid="stPlotlyChart"],
  [data-testid="stMain"] [data-testid="stDataFrame"],
  [data-testid="stMain"] [data-testid="stTable"] {{
    animation: cafe-fade-in var(--cafe-duration-fade) ease-out;
  }}

  /* --- Sidebar navigation: slide + tint on hover ------------------------- */
  /* Makes the primary way of moving through the app feel responsive. */
  [data-testid="stSidebarNav"] a {{
    transition: transform var(--cafe-duration-micro) var(--cafe-expo),
                background-color var(--cafe-duration-micro) ease-out;
  }}
  [data-testid="stSidebarNav"] a:hover {{
    transform: translateX(3px);
  }}

  /* --- File uploader: the dropzone is a primary action, so make it react -- */
  [data-testid="stFileUploaderDropzone"] {{
    transition: transform var(--cafe-duration-micro) var(--cafe-expo),
                border-color var(--cafe-duration-micro) ease-out,
                box-shadow var(--cafe-duration-micro) ease-out,
                background-color var(--cafe-duration-micro) ease-out;
  }}
  [data-testid="stFileUploaderDropzone"]:hover {{
    transform: translateY(-1px);
    border-color: rgba(15, 80, 112, 0.55);
    box-shadow: 0 6px 18px rgba(15, 80, 112, 0.12);
  }}

  /* --- Tabs: glide the underline and fade the switched-in panel ----------- */
  [data-testid="stMain"] [data-baseweb="tab-highlight"] {{
    transition: all var(--cafe-duration-tab) var(--cafe-expo) !important;
  }}
  [data-testid="stMain"] [data-baseweb="tab"] {{
    transition: color var(--cafe-duration-micro) ease-out, background-color var(--cafe-duration-micro) ease-out;
  }}
  [data-testid="stMain"] [data-testid="stTabPanel"] {{
    animation: cafe-fade-in var(--cafe-duration-tab) ease-out;
  }}

  /* --- Expander: soften the header hover --------------------------------- */
  [data-testid="stMain"] [data-testid="stExpander"] summary {{
    transition: background-color var(--cafe-duration-micro) ease-out;
  }}
  [data-testid="stMain"] [data-testid="stExpander"] summary:hover {{
    background-color: var(--cafe-surface-hover);
  }}

  /* --- Images / plots: fade in so freshly rendered figures don't pop ------ */
  /* Charts are re-rendered whenever a parameter changes; fading the new image
     in is quiet confirmation that the figure updated. Opacity-only, so it's
     GPU-cheap even for large plots. */
  [data-testid="stMain"] [data-testid="stImage"] img {{
    animation: cafe-fade-in var(--cafe-duration-fade) ease-out;
  }}
}}

@keyframes cafe-rise-in {{
  from {{ opacity: 0; transform: translateY(-4px); }}
  to {{ opacity: 1; transform: translateY(0); }}
}}

@keyframes cafe-fade-in {{
  from {{ opacity: 0; }}
  to {{ opacity: 1; }}
}}

@keyframes cafe-sheen {{
  from {{ background-position: 220% 0; }}
  to {{ background-position: -120% 0; }}
}}

@keyframes cafe-hero-in {{
  from {{ opacity: 0; transform: translateY(14px); }}
  to {{ opacity: 1; transform: translateY(0); }}
}}
</style>
"""

# Landing-page entrance (one-shot per session)
_HERO_CSS = f"""
<style>
@media (prefers-reduced-motion: no-preference) {{
  [data-testid="stMain"] .block-container {{
    animation: cafe-hero-in var(--cafe-duration-hero) var(--cafe-expo) both;
  }}
}}
</style>
"""

# Landing-page ambience (dot field + scroll reveals)
_LANDING_CSS = f"""
<style>
[data-testid="stApp"]::before {{
  content: "";
  position: fixed;
  inset: 0;
  z-index: 0;
  pointer-events: none;
  background-image: radial-gradient(rgba(15, 80, 112, 0.10) 1px, transparent 1.4px);
  background-size: 26px 26px;
  opacity: 0.5;
  -webkit-mask-image: linear-gradient(180deg, #000 0%, transparent 55%);
          mask-image: linear-gradient(180deg, #000 0%, transparent 55%);
}}
@media (prefers-reduced-motion: no-preference) {{
  [data-testid="stApp"]::before {{
    animation: cafe-drift 26s linear infinite alternate;
  }}
}}
/* Keep real content above the texture. */
[data-testid="stMain"] .block-container {{ position: relative; z-index: 1; }}

@supports (animation-timeline: view()) {{
  @media (prefers-reduced-motion: no-preference) {{
    [data-testid="stMain"] .block-container .element-container {{
      animation: cafe-reveal linear both;
      animation-timeline: view();
      animation-range: entry 0% cover 22%;
    }}
  }}
}}

@keyframes cafe-drift {{
  from {{ background-position: 0 0; }}
  to {{ background-position: 90px 140px; }}
}}
@keyframes cafe-reveal {{
  from {{ opacity: 0; transform: translateY(20px); }}
  to {{ opacity: 1; transform: translateY(0); }}
}}
</style>
"""

# Stepper CSS (used by Data Processing wizard)
_STEPPER_CSS = """
<style>
.cafe-stepper {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: var(--cafe-space-sm);
  margin-bottom: var(--cafe-space-xl);
  padding: var(--cafe-space-md) 0;
  border-bottom: 1px solid var(--cafe-border);
  overflow-x: auto;
}
.cafe-stepper__step {
  display: flex;
  align-items: center;
  gap: var(--cafe-space-sm);
  white-space: nowrap;
  font-size: 0.875rem;
  color: var(--cafe-text-subtle);
  transition: color 0.18s ease-out;
}
.cafe-stepper__step--current {
  color: var(--cafe-text);
  font-weight: 600;
}
.cafe-stepper__step--done {
  color: var(--cafe-text-muted);
}
.cafe-stepper__dot {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 24px;
  height: 24px;
  border-radius: 50%;
  border: 2px solid var(--cafe-border);
  font-size: 0.75rem;
  font-weight: 700;
  color: var(--cafe-text-subtle);
  background: var(--cafe-white);
  transition: border-color 0.18s ease-out, background-color 0.18s ease-out, color 0.18s ease-out;
  flex-shrink: 0;
}
.cafe-stepper__step--current .cafe-stepper__dot {
  border-color: var(--cafe-primary);
  background: var(--cafe-primary);
  color: var(--cafe-white);
}
.cafe-stepper__step--done .cafe-stepper__dot {
  border-color: var(--cafe-primary);
  background: var(--cafe-primary);
  color: var(--cafe-white);
}
.cafe-stepper__label {
  line-height: 1.2;
}
</style>
"""

# Info-card / section helpers
_INFO_CARD_CSS = """
<style>
.cafe-info-card {
  padding: var(--cafe-space-md);
  border-radius: 8px;
  border: 1px solid var(--cafe-border);
  background: var(--cafe-surface);
  margin-bottom: var(--cafe-space-md);
}

.cafe-info-card__title {
  margin: 0 0 var(--cafe-space-xs) 0;
  font-size: 1rem;
  font-weight: 600;
  color: var(--cafe-text);
}
.cafe-info-card__title + .cafe-info-card__body {
  margin-top: var(--cafe-space-xs);
}
.cafe-info-card__body {
  margin: 0;
  font-size: 0.9375rem;
  line-height: 1.5;
  color: var(--cafe-text-muted);
}
.cafe-info-card__body p:last-child {
  margin-bottom: 0;
}
.cafe-section-header {
  margin-bottom: var(--cafe-space-md);
}
.cafe-section-header__title {
  margin: 0;
  font-size: 1.25rem;
  font-weight: 600;
  color: var(--cafe-text);
}
.cafe-section-header__subtitle {
  margin: var(--cafe-space-xs) 0 0 0;
  font-size: 0.9375rem;
  color: var(--cafe-text-muted);
}
</style>
"""


# Merge-builder CSS (used by Merge Clusters page)
_MERGE_BUILDER_CSS = """
<style>
/* ── Cluster chip ─────────────────────────────────────────────────────────── */
.cafe-chip {
  display: inline-flex;
  align-items: center;
  padding: 2px 10px;
  border-radius: 99px;
  font-size: 0.8125rem;
  font-weight: 600;
  line-height: 1.6;
  white-space: nowrap;
  font-family: var(--cafe-font);
}
.cafe-chip--primary {
  background: var(--cafe-primary-light);
  color: var(--cafe-primary);
}
.cafe-chip--muted {
  background: var(--cafe-surface);
  color: var(--cafe-text-muted);
  border: 1px solid var(--cafe-border);
}
.cafe-chip--assigned {
  background: var(--cafe-surface);
  color: var(--cafe-text-subtle);
  border: 1px solid var(--cafe-border);
  text-decoration: line-through;
  opacity: 0.65;
}

/* ── Cluster inventory bar ────────────────────────────────────────────────── */
.cafe-cluster-inventory {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  padding: var(--cafe-space-md);
  border: 1px solid var(--cafe-border);
  border-radius: 8px;
  background: var(--cafe-surface);
  margin-bottom: var(--cafe-space-md);
}
.cafe-cluster-inventory__label {
  width: 100%;
  font-size: 0.8125rem;
  font-weight: 600;
  color: var(--cafe-text-muted);
  text-transform: uppercase;
  letter-spacing: 0.06em;
  margin-bottom: 2px;
}

/* ── Merge group preview card ─────────────────────────────────────────────── */
.cafe-merge-group {
  display: flex;
  align-items: flex-start;
  gap: var(--cafe-space-sm);
  padding: var(--cafe-space-sm) var(--cafe-space-md);
  border: 1px solid var(--cafe-border);
  border-radius: 8px;
  background: var(--cafe-white);
  margin-bottom: var(--cafe-space-sm);
}
@media (prefers-reduced-motion: no-preference) {
  .cafe-merge-group {
    animation: cafe-rise-in 0.28s var(--cafe-expo) both;
  }
}
.cafe-merge-group__arrow {
  font-size: 1rem;
  color: var(--cafe-primary);
  padding-top: 2px;
  flex-shrink: 0;
  font-weight: 700;
}
.cafe-merge-group__target {
  font-size: 0.9375rem;
  font-weight: 700;
  color: var(--cafe-text);
  flex-shrink: 0;
  padding-top: 1px;
  min-width: 60px;
}
.cafe-merge-group__sources {
  display: flex;
  flex-wrap: wrap;
  gap: 5px;
}
.cafe-merge-group--empty {
  border-style: dashed;
  background: var(--cafe-surface);
}
.cafe-merge-group--empty .cafe-merge-group__target {
  color: var(--cafe-text-muted);
}

/* ── Format string preview ────────────────────────────────────────────────── */
.cafe-format-preview {
  padding: var(--cafe-space-sm) var(--cafe-space-md);
  border-radius: 8px;
  background: var(--cafe-surface);
  border: 1px solid var(--cafe-border);
  font-size: 0.875rem;
  font-family: "SFMono-Regular", "Consolas", "Liberation Mono", monospace;
  color: var(--cafe-text-muted);
  word-break: break-all;
  margin-top: var(--cafe-space-sm);
}
.cafe-format-preview__label {
  font-size: 0.75rem;
  font-weight: 600;
  color: var(--cafe-text-subtle);
  text-transform: uppercase;
  letter-spacing: 0.06em;
  margin-bottom: 4px;
  font-family: var(--cafe-font);
}
.cafe-format-preview__value {
  color: var(--cafe-primary);
  font-weight: 600;
}

/* ── Warning badge ────────────────────────────────────────────────────────── */
.cafe-warn-badge {
  display: inline-flex;
  align-items: center;
  gap: 5px;
  padding: 3px 10px;
  border-radius: 99px;
  font-size: 0.8125rem;
  font-weight: 600;
  background: #fff3e0;
  color: var(--cafe-warning);
  border: 1px solid #ffe0b2;
}
</style>
"""


# Public API: injection functions
def apply_theme():
    """Inject the global design system and micro-interaction layer."""
    st.markdown(_CSS_VARIABLES, unsafe_allow_html=True)
    st.markdown(_BASE_CSS, unsafe_allow_html=True)
    st.markdown(_PREMIUM_CSS, unsafe_allow_html=True)
    st.markdown(_INFO_CARD_CSS, unsafe_allow_html=True)
    st.markdown(_STEPPER_CSS, unsafe_allow_html=True)
    st.markdown(_MERGE_BUILDER_CSS, unsafe_allow_html=True)


def landing_ambient():
    """Inject the welcome-page backdrop + scroll reveals (safe to omit; page is then static)."""
    st.markdown(_LANDING_CSS, unsafe_allow_html=True)


def hero_entrance():
    """Play the landing-page entrance once per session (guarded so reruns don't replay it)."""
    if st.session_state.get("_cafe_hero_shown"):
        return
    st.session_state["_cafe_hero_shown"] = True
    st.markdown(_HERO_CSS, unsafe_allow_html=True)


# Public API: component helpers
def page_header(title: str, subtitle: str | None = None):
    """Render a consistent page header with optional subtitle."""
    st.title(title)
    if subtitle:
        st.caption(subtitle)


def section_header(title: str, subtitle: str | None = None, step: int | None = None):
    """Render a styled section header with optional step number."""
    prefix = f"{step}. " if step is not None else ""
    subtitle_html = (
        f'<p class="cafe-section-header__subtitle">{subtitle}</p>' if subtitle else ""
    )
    st.markdown(
        f'<div class="cafe-section-header">'
        f'<h3 class="cafe-section-header__title">{prefix}{title}</h3>'
        f"{subtitle_html}"
        f"</div>",
        unsafe_allow_html=True,
    )


def info_card(title: str, body: str, kind: str = "info"):
    """Render a styled info/warning/error/success card."""
    allowed = {"info", "warning", "error", "success"}
    kind = kind if kind in allowed else "info"
    st.markdown(
        f'<div class="cafe-info-card cafe-info-card--{kind}">'
        f'<h4 class="cafe-info-card__title">{title}</h4>'
        f'<div class="cafe-info-card__body">{body}</div>'
        f"</div>",
        unsafe_allow_html=True,
    )


def stepper(steps: list[str], current_index: int, done_indices: set[int] | None = None):
    """Render a horizontal wizard stepper (current_index active, done_indices ticked)."""
    done_indices = done_indices or set()
    items = []
    for i, label in enumerate(steps):
        if i in done_indices:
            state = "done"
            dot = "✓"
        elif i == current_index:
            state = "current"
            dot = str(i + 1)
        else:
            state = ""
            dot = str(i + 1)
        items.append(
            f'<div class="cafe-stepper__step cafe-stepper__step--{state}">'
            f'<span class="cafe-stepper__dot">{dot}</span>'
            f'<span class="cafe-stepper__label">{label}</span>'
            f"</div>"
        )
    st.markdown(
        f'<div class="cafe-stepper">{"".join(items)}</div>',
        unsafe_allow_html=True,
    )


def loading_status(title: str, steps: list[str], current_index: int | None = None):
    """Render a labeled progress block; current_index highlights the running step, None shows just the title."""
    if current_index is None:
        st.markdown(f"**{title}**")
        return
    progress = (current_index + 1) / len(steps) if steps else 0
    st.progress(progress, text=f"{title}: {steps[current_index]}")
