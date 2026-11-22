"""Shared design tokens, global CSS, and UI helpers for CAFE (mirrors the docs site: teal accent, IBM Plex Sans / Instrument Sans, reduced-motion aware)."""
import os
from dataclasses import dataclass

import streamlit as st

# Bundled image assets, resolved relative to this file (cwd-independent).
IMG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "img")


@dataclass(frozen=True)
class Tokens:
    """Single source of truth for CAFE design tokens."""

    # Color (light values of the docs site's tokens in src/styles/custom.css)
    primary: str = "#0c5a7a"
    primary_hover: str = "#083f57"
    primary_light: str = "#ddedf3"
    primary_lighter: str = "#eef6f9"
    background: str = "#ffffff"
    surface: str = "#f6f8f9"
    surface_hover: str = "#e2e9ec"
    border: str = "#d7e0e4"
    border_strong: str = "#c5d1d5"
    text: str = "#0b1a20"
    text_muted: str = "#53656c"
    text_subtle: str = "#627579"  # 4.5:1 on surface, AA for small text
    spectrum_1: str = "#0e9bbd"
    spectrum_2: str = "#3f6fdc"
    spectrum_3: str = "#7b4fd4"
    success: str = "#2d6a4f"
    warning: str = "#b35900"
    error: str = "#c1121f"
    white: str = "#ffffff"

    # Typography
    font: str = '"IBM Plex Sans", ui-sans-serif, system-ui, sans-serif'
    font_display: str = '"Instrument Sans", "IBM Plex Sans", ui-sans-serif, sans-serif'
    font_mono: str = '"JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, monospace'

    # Shape
    radius_sm: str = "8px"
    radius: str = "12px"
    radius_lg: str = "20px"
    shadow_sm: str = "0 1px 2px rgb(11 34 44 / 0.06), 0 1px 1px rgb(11 34 44 / 0.04)"
    shadow: str = "0 6px 16px -6px rgb(11 34 44 / 0.12), 0 2px 4px -2px rgb(11 34 44 / 0.06)"
    shadow_lg: str = "0 28px 56px -24px rgb(11 34 44 / 0.22), 0 8px 16px -8px rgb(11 34 44 / 0.08)"

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
@import url('https://fonts.bunny.net/css?family=ibm-plex-sans:400,400i,500,600|instrument-sans:500,600,700|jetbrains-mono:400,500&display=swap');

:root, [data-testid="stApp"] {{
  --cafe-primary: {TOKENS.primary};
  --cafe-primary-hover: {TOKENS.primary_hover};
  --cafe-primary-light: {TOKENS.primary_light};
  --cafe-primary-lighter: {TOKENS.primary_lighter};
  --cafe-background: {TOKENS.background};
  --cafe-surface: {TOKENS.surface};
  --cafe-surface-hover: {TOKENS.surface_hover};
  --cafe-border: {TOKENS.border};
  --cafe-border-strong: {TOKENS.border_strong};
  --cafe-spectrum-1: {TOKENS.spectrum_1};
  --cafe-spectrum-2: {TOKENS.spectrum_2};
  --cafe-spectrum-3: {TOKENS.spectrum_3};
  --cafe-font-display: {TOKENS.font_display};
  --cafe-font-mono: {TOKENS.font_mono};
  --cafe-radius-sm: {TOKENS.radius_sm};
  --cafe-radius: {TOKENS.radius};
  --cafe-radius-lg: {TOKENS.radius_lg};
  --cafe-shadow-sm: {TOKENS.shadow_sm};
  --cafe-shadow: {TOKENS.shadow};
  --cafe-shadow-lg: {TOKENS.shadow_lg};
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

/* --- Docs-site look: display headings, ink text, hairline surfaces ------- */
html, body {{ -webkit-font-smoothing: antialiased; }}
::selection {{ background: rgba(12, 90, 122, 0.22); }}

[data-testid="stMain"] h1,
[data-testid="stMain"] h2,
[data-testid="stMain"] h3,
[data-testid="stMain"] h4 {{
  font-family: var(--cafe-font-display) !important;
  color: var(--cafe-text);
  letter-spacing: -0.02em;
}}
[data-testid="stMain"] h1 {{ font-weight: 600 !important; }}

code, pre, kbd, [data-testid="stCode"] * {{ font-family: var(--cafe-font-mono) !important; }}
[data-testid="stMarkdownContainer"] code {{
  background: var(--cafe-surface);
  border: 1px solid var(--cafe-border);
  border-radius: 6px;
  padding: 0.1em 0.35em;
  color: var(--cafe-text);
  font-size: 0.85em;
}}

/* Sidebar: white rail with a hairline, quiet group labels, tinted active pill */
[data-testid="stSidebar"] {{
  background: var(--cafe-background) !important;
  border-right: 1px solid var(--cafe-border);
}}
[data-testid="stSidebarNav"] header,
[data-testid="stNavSectionHeader"] {{
  font-size: 0.6875rem !important;
  font-weight: 600 !important;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--cafe-text-subtle) !important;
}}
[data-testid="stSidebarNav"] a {{
  border-radius: var(--cafe-radius-sm);
}}
/* Navigation stays open: no collapse button. The expand control is kept so a
   sidebar auto-collapsed on a narrow window can still be reopened. */
[data-testid="stSidebarCollapseButton"] {{ display: none !important; }}
[data-testid="stSidebarNav"] a[aria-current="page"] {{
  background: var(--cafe-primary-light) !important;
}}
[data-testid="stSidebarNav"] a[aria-current="page"] span {{
  color: var(--cafe-primary) !important;
  font-weight: 600;
}}

/* Buttons: pills like the site's CTAs, visible offset focus ring */
button[data-testid*="stBaseButton"] {{ border-radius: 999px !important; }}
button[data-testid*="stBaseButton-primary"] {{ box-shadow: 0 6px 16px -8px rgba(12, 90, 122, 0.6); }}
button[data-testid*="stBaseButton"]:focus-visible {{
  outline: 2px solid var(--cafe-primary) !important;
  outline-offset: 2px;
}}

/* Number steppers: white pill field with round outlined − / + buttons (app-wide). */
[data-testid="stNumberInputContainer"] {{
  background: var(--cafe-background) !important;
  border: 1px solid var(--cafe-border-strong) !important;
  border-radius: 999px !important;
  padding-right: 4px;
}}
[data-testid="stNumberInputContainer"]:focus-within {{ border-color: var(--cafe-primary) !important; }}
[data-testid="stNumberInputField"] {{ background: transparent !important; font-variant-numeric: tabular-nums; }}
[data-testid="stNumberInputStepDown"],
[data-testid="stNumberInputStepUp"] {{
  width: 28px !important; height: 28px !important; margin: auto 2px;
  border: 1px solid var(--cafe-border) !important; border-radius: 999px !important;
  background: var(--cafe-surface) !important; color: var(--cafe-text) !important;
}}
[data-testid="stNumberInputStepDown"]:hover:not(:disabled),
[data-testid="stNumberInputStepUp"]:hover:not(:disabled) {{
  background: var(--cafe-primary-light) !important; border-color: var(--cafe-primary) !important;
  color: var(--cafe-primary) !important;
}}

/* Raised panels: alerts, expanders, uploader, tables share the site's card shape */
[data-testid="stMain"] [data-testid="stAlert"] > div,
[data-testid="stMain"] [data-testid="stExpander"] details,
[data-testid="stFileUploaderDropzone"],
[data-testid="stMain"] [data-testid="stDataFrame"],
[data-testid="stMain"] [data-testid="stTable"] {{
  border-radius: var(--cafe-radius) !important;
}}
[data-testid="stMain"] [data-testid="stExpander"] details {{
  border-color: var(--cafe-border) !important;
  box-shadow: var(--cafe-shadow-sm);
}}
[data-testid="stFileUploaderDropzone"] {{
  border: 1px dashed var(--cafe-border-strong) !important;
  background: var(--cafe-surface) !important;
}}
[data-testid="stMain"] hr {{ border-color: var(--cafe-border) !important; }}
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
    /* Mirrors .streamlit/config.toml's primaryColor (#0c5a7a) */
    box-shadow: 0 8px 20px -8px rgba(12, 90, 122, 0.55);
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
    border-color: rgba(12, 90, 122, 0.55);
    box-shadow: var(--cafe-shadow);
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

# Landing-page hero styles + scroll reveals (plain white page, no backdrop)
_LANDING_CSS = f"""
<style>
/* Welcome hero: same badge / display headline / pill links as the docs landing page. */
.cafe-hero {{ padding: var(--cafe-space-md) 0 var(--cafe-space-sm); }}
.cafe-hero__badge {{
  display: inline-flex; align-items: center; gap: 8px;
  margin: 0 0 var(--cafe-space-md) !important;
  padding: 4px 12px 4px 10px;
  border: 1px solid var(--cafe-border); border-radius: 999px;
  background: var(--cafe-background); box-shadow: var(--cafe-shadow-sm);
  font-size: 0.8125rem !important; font-weight: 500; color: var(--cafe-text);
}}
.cafe-hero__dot {{
  width: 8px; height: 8px; border-radius: 50%;
  background: var(--cafe-primary); box-shadow: 0 0 0 3px var(--cafe-primary-light);
}}
.cafe-hero__year {{
  padding-left: 8px; border-left: 1px solid var(--cafe-border);
  font-family: var(--cafe-font-mono); font-size: 0.75rem; color: var(--cafe-text-muted);
}}
[data-testid="stMain"] h1.cafe-hero__title {{
  margin: 0 0 var(--cafe-space-md); padding: 0;
  font-size: clamp(2.25rem, 5vw, 3.25rem) !important; line-height: 1.05;
  letter-spacing: -0.035em;
}}
.cafe-hero__lede {{
  max-width: 38rem; font-size: 1.0625rem !important; line-height: 1.6;
  color: var(--cafe-text-muted);
}}
.cafe-hero__actions {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: var(--cafe-space-lg) !important; }}
.cafe-hero__btn {{
  display: inline-flex; align-items: center; min-height: 40px; padding: 0 18px;
  border: 1px solid var(--cafe-border-strong); border-radius: 999px;
  background: var(--cafe-background); color: var(--cafe-text) !important;
  font-weight: 600; font-size: 0.9375rem; text-decoration: none !important;
  box-shadow: var(--cafe-shadow-sm);
}}
.cafe-hero__btn:hover {{ border-color: var(--cafe-primary); }}
.cafe-hero__btn:focus-visible {{ outline: 2px solid var(--cafe-primary); outline-offset: 2px; }}
.cafe-hero__btn--primary {{
  background: var(--cafe-primary); border-color: var(--cafe-primary); color: var(--cafe-white) !important;
  box-shadow: 0 8px 20px -10px rgba(12, 90, 122, 0.7);
}}
.cafe-hero__btn--primary:hover {{ background: var(--cafe-primary-hover); }}
.cafe-hero__btn--ghost {{ border-color: transparent; background: transparent; box-shadow: none; }}


@supports (animation-timeline: view()) {{
  @media (prefers-reduced-motion: no-preference) {{
    [data-testid="stMain"] .block-container .element-container {{
      animation: cafe-reveal linear both;
      animation-timeline: view();
      animation-range: entry 0% cover 22%;
    }}
  }}
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
  box-shadow: 0 0 0 4px var(--cafe-primary-light);
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
  padding: var(--cafe-space-md) var(--cafe-space-lg);
  border-radius: var(--cafe-radius);
  box-shadow: var(--cafe-shadow-sm);
  border: 1px solid var(--cafe-border);
  background: var(--cafe-surface);
  margin-bottom: var(--cafe-space-md);
}

.cafe-info-card__title {
  margin: 0 0 var(--cafe-space-xs) 0;
  font-family: var(--cafe-font-display);
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
  font-family: var(--cafe-font-display);
  letter-spacing: -0.01em;
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
  font-family: var(--cafe-font-mono);
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
    """Inject the welcome-page hero styles + scroll reveals."""
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
