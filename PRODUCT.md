# CAFE — Product Definition

## Product name

CAFE: Cell Analyzer for Flow Experiment

## Register

**Product UI.** CAFE is a no-code scientific tool. Design serves the analysis workflow; it is not a marketing/brand experience.

## Users

Immunology and flow-cytometry researchers, primarily bench scientists analyzing spectral flow cytometry (SFCM) data exported from FlowJo. They are not developers. They need publication-ready figures and reproducible results without writing code.

Typical user context:

- Works on a laptop or lab workstation, often with a mouse and a large screen.
- Uploads per-sample CSV files, runs dimensionality reduction and clustering, and exports static figures for papers.
- Values clarity, precision, and speed over visual spectacle.
- Needs confidence that the tool will not silently alter their data or figures.

## Use cases

1. **Upload and process CSVs.** Convert per-sample FlowJo exports into a single AnnData object.
2. **Quality control.** Remove doublets, debris, and unmixing outliers using MAD-based outlier detection.
3. **Batch correction.** Correct for technical batch effects while preserving biological groups.
4. **Dimensionality reduction and clustering.** Run PCA, UMAP, and Leiden clustering with adjustable parameters.
5. **Interactive exploration.** Pan, zoom, and select cells in an embedding to profile clusters and groups.
6. **Statistical comparison.** Compare cluster abundance and marker expression across experimental groups.
7. **Export publication figures.** Download 300 DPI PNG/PDF/SVG plots and a ZIP of all results.

## Product purpose

Make high-dimensional spectral flow cytometry analysis accessible to non-coders while preserving the rigor and figure quality expected by peer reviewers.

## Brand personality / tone

**Serious lab tool.** Precise, trustworthy, and restrained. The "wow" should come from capability — fast computation, clear results, reliable exports — not from visual flash.

- **Voice:** Direct, human, and sparing. Every word earns its place.
- **Mood:** Calm, focused, and professional. The interface should feel like a well-made instrument, not a dashboard template.
- **Motion:** Purposeful and minimal. Animations are used to signal state changes and progress, never to decorate.
- **Data integrity:** Interactive and animated elements must never replace or alter the downloadable static publication figures.

## Anti-references

Avoid these patterns:

- Generic SaaS dashboards with big-number hero metrics and gradient text.
- Dark-mode tools with glowing neon accents chosen "because science looks cool dark."
- Overly playful or friendly chatbot-style copy in a scientific tool.
- Glassmorphism, heavy drop shadows, or decorative blur.
- Identical icon + heading + text card grids.
- Modal-first interactions when inline or progressive disclosure would work.
- Category-reflex palettes (e.g., healthcare teal, finance navy, crypto neon).

## Strategic principles

1. **Science first.** The interface must make the data pipeline obvious, not hide it behind abstraction.
2. **Publication integrity.** All downloadable figures remain static matplotlib/Plotly exports at 300 DPI. Interactive views are additive companions, not replacements.
3. **Progressive disclosure.** Start simple; reveal advanced parameters only when the user asks for them.
4. **Responsiveness within constraints.** The app is built in Streamlit, so full mobile-native responsiveness is not achievable. The design should still be usable on tablets and small laptops by avoiding cramped multi-column layouts and tiny touch targets.
5. **Accessibility by default.** Respect `prefers-reduced-motion`, maintain strong contrast, and never rely on color alone to convey meaning.
6. **Consistency across surfaces.** Desktop app, hosted web app, and documentation site should share the same visual language.

## Constraints

- Streamlit app; locked light theme because charts are rendered with light backgrounds.
- Datasets typically <100,000 cells; Plotly ScatterGL is sufficient for interactivity.
- No automated test suite; validation is manual.
- All motion is gated behind `prefers-reduced-motion`.
- Python/Scanpy/AnnData backend must remain intact.
