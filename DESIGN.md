# CAFE — Design System

## Design strategy

**Restrained + one confident accent.** The interface is dominated by warm, tinted neutrals. The teal primary color (`#0f5070`) carries the action: buttons, active states, progress, and the stepper. Chroma is reduced at the extremes so nothing feels garish against the light background.

## Theme

**Light, calm, laboratory.** The physical scene is a researcher working at a laptop or monitor in a brightly lit lab or office, reviewing cell populations and preparing figures. A light theme keeps charts (which are always light-background static images) looking native to the page and avoids the eye strain of bright plots on dark backgrounds.

## Color palette

| Token | Hex | Usage |
|-------|-----|-------|
| `--cafe-primary` | `#0f5070` | Primary buttons, active states, stepper current step, progress bars, links |
| `--cafe-primary-hover` | `#0a3d57` | Button hover, pressed state |
| `--cafe-primary-light` | `#e6f1f4` | Tinted backgrounds for selected/highlighted items |
| `--cafe-primary-lighter` | `#f0f6f8` | Subtle surface tints |
| `--cafe-background` | `#ffffff` | Page background |
| `--cafe-surface` | `#f8f9fa` | Cards, expanders, form sections |
| `--cafe-surface-hover` | `#e9ecef` | Hover state on surface elements |
| `--cafe-border` | `#d5e0e4` | Dividers, borders, inactive stepper rings |
| `--cafe-text` | `#1a262a` | Primary text, headings |
| `--cafe-text-muted` | `#5c6e74` | Secondary text, captions, labels |
| `--cafe-text-subtle` | `#8fa3aa` | Disabled/future stepper text |
| `--cafe-success` | `#2d6a4f` | Success states, positive indicators |
| `--cafe-warning` | `#b35900` | Warnings, cautionary notes |
| `--cafe-error` | `#c1121f` | Errors, destructive actions |
| `--cafe-white` | `#ffffff` | Pure white for charts, dropdowns, contrast surfaces |

**Color strategy notes:**

- Neutrals are tinted toward the teal hue (roughly hue 200–210) so the page feels cohesive without being blue.
- Pure black and pure white are avoided in UI chrome; `#ffffff` is reserved for chart backgrounds and dropdowns where Streamlit expects it.
- Data visualizations may use non-brand palettes (e.g., `tab20c`, `RdBu_r`) when they encode categorical or continuous data. The brand applies to UI chrome, not to data encoding.

## Typography

**Font family:** Source Sans 3 (loaded via Bunny Fonts or Google Fonts).

- Clean, readable, and neutral enough for dense scientific interfaces.
- Not Inter, Roboto, Arial, or system defaults.
- Matplotlib and Plotly figures default to the same family where possible; otherwise they use Arial for backward compatibility with publication workflows.

**Type scale:**

| Level | Size | Line height | Weight | Usage |
|-------|------|-------------|--------|-------|
| H1 | `2rem` (32px) | 1.2 | 600 | Page titles |
| H2 | `1.5rem` (24px) | 1.3 | 600 | Section headers |
| H3 | `1.25rem` (20px) | 1.35 | 600 | Subsection headers |
| Body | `1rem` (16px) | 1.5 | 400 | Paragraphs, labels, tables |
| Caption | `0.8125rem` (13px) | 1.4 | 400 | Help text, captions, meta |
| Mono | `0.875rem` (14px) | 1.4 | 400 | Code, sample IDs, file names |

**Typography rules:**

- Body line length is capped at ~75ch in reading-heavy sections.
- Hierarchy is created through weight and size contrast, not all-caps or underline.
- Headings are left-aligned. Centered headings are reserved for the splash/landing page only.

## Spacing

Base unit: **4px**. Use a 4-point grid for all spacing.

| Token | Value | Usage |
|-------|-------|-------|
| `--space-xs` | 4px | Tight internal gaps, icon padding |
| `--space-sm` | 8px | Between related inline controls |
| `--space-md` | 16px | Default gap between form elements |
| `--space-lg` | 24px | Between sections within a card |
| `--space-xl` | 32px | Between major page sections |
| `--space-2xl` | 48px | Between page-level blocks |
| `--space-3xl` | 64px | Hero/landing section separation |

**Spacing rules:**

- Vary spacing for rhythm; do not use the same padding everywhere.
- Tighten related controls; separate unrelated sections generously.
- Avoid nested cards. Use surface tints and borders to group related content instead.

## Elevation

CAFE avoids heavy shadows. Depth is communicated through:

- **Surface tints** (`--cafe-surface`, `--cafe-surface-hover`) for grouping.
- **1px borders** (`--cafe-border`) for subtle separation.
- **One soft shadow** for primary button hover and focused cards: `0 4px 14px rgba(15, 80, 112, 0.12)`.

No glassmorphism, no multi-layer floating cards, no dark blurred overlays.

## Components

### Buttons

- **Primary:** teal background, white text, rounded `6px` corners, hover lift `translateY(-1px)` + teal shadow, active press `scale(0.98)`.
- **Secondary:** transparent background, teal border, teal text, hover surface tint.
- **Danger:** reserved for destructive actions; use the error token.

### Inputs / file uploader

- Use Streamlit defaults, but override the dropzone hover state with teal border and shadow.
- Keep labels concise and place help text below, not inside placeholders.

### Cards / sections

- A card is a surface-tinted container with a 1px border and `12px`–`16px` internal padding.
- Do not nest cards inside cards. Use vertical spacing and headings to separate content.

### Alerts / info cards

- **Info:** teal left border or teal icon on surface background.
- **Warning:** amber left border or icon.
- **Error:** red left border or icon.
- **Success:** green left border or icon.
- No side-stripe borders greater than 1px; use a full 1px border or a tinted background instead.

### Progress

- Teal fill with a slow white sheen animation for long-running steps.
- Step text updates to tell the user what is currently running.
- Fallback to a determinate progress bar when percentage is known.

### Stepper (Data Processing)

- Horizontal bar at the top of the page.
- Steps: Upload → QC → PCA → Batch → UMAP/Leiden → Export.
- Completed step: teal checkmark, muted label.
- Current step: filled teal circle, bold label.
- Future step: empty gray circle, subtle label.
- Responsive fallback: collapse to a vertical list on very narrow screens.

## Motion

**Easing:** `cubic-bezier(0.16, 1, 0.3, 1)` (ease-out-expo) for movement; plain `ease-out` for color/opacity.

**Durations:**

| Type | Duration |
|------|----------|
| Micro-interaction (hover, focus) | 0.18s |
| Fade-in (charts, tables) | 0.35s |
| Tab / panel switch | 0.28s |
| Hero entrance | 0.6s |
| Progress sheen | 1.6s linear loop |

**Rules:**

- Only animate `transform` and `opacity`. Never animate `width`, `height`, `margin`, `padding`, or `top`/`left`.
- No bounce, no elastic, no decorative parallax.
- All motion is gated behind `@media (prefers-reduced-motion: no-preference)`.

## Accessibility

- Text contrast ≥ 4.5:1 for normal text, ≥ 3:1 for large text.
- Focus indicators are visible on all custom interactive elements.
- Color is never the only indicator of state (e.g., stepper uses icons + labels).
- `prefers-reduced-motion` disables all animations.
- Form inputs have associated labels.

## Responsive approach

CAFE is a Streamlit app, so true mobile responsiveness is limited. The design uses a **best-effort responsive** strategy:

- **Desktop (≥1024px):** full sidebar, multi-column layouts where helpful, large charts.
- **Tablet (768px–1023px):** single-column form layouts, smaller charts, simplified stepper.
- **Small screens (<768px):** single column, stacked controls, minimum 44px touch targets, vertical stepper.
- Avoid 3-column control groups on narrow screens.
- Test manually at 375px, 768px, and 1440px viewports.

## File locations

- Design tokens and global CSS: `CAFE/bin/theme.py`
- App configuration: `CAFE/.streamlit/config.toml`
- Product context: `PRODUCT.md` (this repo)
- Design context: `DESIGN.md` (this repo)
- Persisted AI context: `.impeccable.md`
