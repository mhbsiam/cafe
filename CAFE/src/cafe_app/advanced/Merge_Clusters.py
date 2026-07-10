import io
import os
import random
import tempfile
import zipfile

import matplotlib.pyplot as plt
import scanpy as sc
import streamlit as st

from theme import IMG_DIR, apply_theme, page_header, section_header

st.set_page_config(layout="centered")
apply_theme()

st.logo(os.path.join(IMG_DIR, 's_logo.png'))

page_header(
    "Merge Clusters",
    subtitle="Combine Leiden subclusters into metaclusters and export the updated AnnData file.",
)

# ── Session-state bootstrap ───────────────────────────────────────────────────

if 'adata' not in st.session_state:
    st.session_state.adata = None

# merge_groups: list of {"target": str, "sources": list[str]}
if 'merge_groups' not in st.session_state:
    st.session_state.merge_groups = [{"target": "", "sources": []}]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _cluster_ids(adata) -> list[str]:
    """Return sorted unique cluster IDs (numeric sort where possible)."""
    raw = adata.obs['leiden'].astype(str).unique().tolist()
    return sorted(raw, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))


def _build_format_string(groups: list[dict]) -> str:
    """Convert the merge_groups list into the {new:old,old} format string."""
    parts = []
    for g in groups:
        t = g.get("target", "").strip()
        s = [x.strip() for x in g.get("sources", []) if x.strip()]
        if t and s:
            parts.append("{" + t + ":" + ",".join(s) + "}")
    return ", ".join(parts)


def _render_inventory_split(cluster_ids: list[str], claimed: set[str]) -> None:
    """Render a two-tier chip bar: remaining (teal) vs assigned (struck-through)."""
    remaining = [c for c in cluster_ids if c not in claimed]
    assigned  = [c for c in cluster_ids if c in claimed]

    n_rem = len(remaining)
    n_asgn = len(assigned)

    remaining_chips = "".join(
        f'<span class="cafe-chip cafe-chip--primary">{c}</span>'
        for c in remaining
    ) or '<span style="font-size:0.8125rem;color:var(--cafe-text-subtle)">none — all clusters assigned</span>'

    assigned_chips = "".join(
        f'<span class="cafe-chip cafe-chip--assigned">{c}</span>'
        for c in assigned
    ) if assigned else ""

    count_label = (
        f'<span class="cafe-cluster-inventory__label">'
        f'{n_rem} remaining'
        f'{f" &nbsp;·&nbsp; {n_asgn} assigned" if n_asgn else ""}'
        f'</span>'
    )

    st.markdown(
        f'<div class="cafe-cluster-inventory">'
        f'{count_label}'
        f'{remaining_chips}'
        f'{("&nbsp;" + assigned_chips) if assigned_chips else ""}'
        f'</div>',
        unsafe_allow_html=True,
    )

    if remaining:
        passthrough_list = ", ".join(remaining)
        st.caption(
            f"Clusters not assigned to any merge group keep their original Leiden number "
            f"unchanged in the output. "
            f"Currently unassigned: **{passthrough_list}**."
        )


def _render_merge_preview(groups: list[dict]) -> None:
    """Render the live merge-group summary cards and format string."""
    has_any = any(
        g.get("target", "").strip() or g.get("sources")
        for g in groups
    )
    if not has_any:
        return

    st.markdown(
        '<p style="font-size:0.8125rem;font-weight:600;color:var(--cafe-text-subtle);'
        'text-transform:uppercase;letter-spacing:0.06em;margin-bottom:8px;">'
        'Merge preview</p>',
        unsafe_allow_html=True,
    )

    for g in groups:
        t = g.get("target", "").strip()
        sources = [x.strip() for x in g.get("sources", []) if x.strip()]
        if not t and not sources:
            continue

        target_label = t if t else '<span style="color:var(--cafe-text-subtle)">unnamed</span>'
        source_chips = "".join(
            f'<span class="cafe-chip cafe-chip--muted">{s}</span>'
            for s in sources
        ) if sources else '<span style="font-size:0.8125rem;color:var(--cafe-text-subtle)">no sources yet</span>'

        empty_cls = "" if (t and sources) else " cafe-merge-group--empty"
        st.markdown(
            f'<div class="cafe-merge-group{empty_cls}">'
            f'<span class="cafe-merge-group__target">{target_label}</span>'
            f'<span class="cafe-merge-group__arrow">&#8592;</span>'
            f'<span class="cafe-merge-group__sources">{source_chips}</span>'
            f'</div>',
            unsafe_allow_html=True,
        )

    fmt = _build_format_string(groups)
    if fmt:
        st.markdown(
            f'<div class="cafe-format-preview">'
            f'<div class="cafe-format-preview__label">Generated format string</div>'
            f'<span class="cafe-format-preview__value">{fmt}</span>'
            f'</div>',
            unsafe_allow_html=True,
        )


def _render_umap(adata, title: str = "") -> None:
    """Render a clean UMAP matplotlib figure with cluster labels."""
    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(
        adata,
        color='leiden',
        palette="tab20c",
        ax=ax,
        show=False,
        legend_loc='on data',
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=11, color='#1a262a', pad=8) if title else ax.set_title('')
    st.pyplot(fig, width="stretch")
    plt.close(fig)


def _validate_merge_groups(groups: list[dict], cluster_ids: list[str]) -> list[str]:
    """Return a list of error strings; empty means valid."""
    errors = []
    seen_sources: dict[str, str] = {}
    for g in groups:
        t = g.get("target", "").strip()
        sources = [x.strip() for x in g.get("sources", []) if x.strip()]
        if not t and not sources:
            continue
        if not t:
            errors.append("One merge group has sources selected but no new cluster name.")
            continue
        if not sources:
            errors.append(f"Merge group \"{t}\" has a name but no source clusters selected.")
            continue
        for s in sources:
            if s not in cluster_ids:
                errors.append(f"Cluster \"{s}\" (in group \"{t}\") does not exist in the loaded data.")
            if s in seen_sources:
                errors.append(
                    f"Cluster \"{s}\" appears in both group \"{seen_sources[s]}\" and \"{t}\". "
                    f"Each source cluster can only be merged once."
                )
            seen_sources[s] = t
    return errors


def _reset_merge_builder() -> None:
    """Reset the builder to one empty row and drop persisted merge_target_*/merge_sources_* widget state."""
    st.session_state.merge_groups = [{"target": "", "sources": []}]
    stale = [
        k for k in list(st.session_state.keys())
        if k.startswith("merge_target_") or k.startswith("merge_sources_")
    ]
    for k in stale:
        del st.session_state[k]


# ── Phase 1: Upload ───────────────────────────────────────────────────────────

with st.form(key='upload_form'):
    uploaded_file = st.file_uploader("Upload AnnData (.h5ad) file", type="h5ad")
    submit_upload = st.form_submit_button("Load AnnData", type="primary")

if uploaded_file and submit_upload:
    with st.spinner("Loading AnnData…"):
        st.session_state.adata = sc.read_h5ad(uploaded_file)
        # Reset the builder (and any prior merge result) for the new file
        _reset_merge_builder()
        st.session_state.pop("merge_result", None)
    if 'X_umap' not in st.session_state.adata.obsm.keys():
        with st.spinner("Computing UMAP (not found in file)…"):
            if 'neighbors' not in st.session_state.adata.uns:
                sc.pp.neighbors(st.session_state.adata)
            sc.tl.umap(st.session_state.adata)
    n_cells, n_markers = st.session_state.adata.shape
    n_clusters = len(_cluster_ids(st.session_state.adata))
    st.success(
        f"Loaded {n_cells:,} cells \u00b7 {n_markers} markers \u00b7 {n_clusters} Leiden clusters"
    )

# ── Phase 2: Workspace ────────────────────────────────────────────────────────

if st.session_state.adata is not None:
    adata = st.session_state.adata
    all_clusters = _cluster_ids(adata)

    st.markdown("---")

    # Find claimed clusters
    all_claimed: set[str] = set()
    for g in st.session_state.merge_groups:
        all_claimed.update(s.strip() for s in g.get("sources", []) if s.strip())

    # Cluster inventory bar — split into remaining vs assigned
    _render_inventory_split(all_clusters, all_claimed)

    # Two-column workspace
    col_map, col_builder = st.columns([1, 1], gap="large")

    with col_map:
        section_header("Cluster map", "Use the cluster numbers on the map to build your merge groups.")
        _render_umap(adata)

    with col_builder:
        section_header("Merge groups", "Define which clusters to combine. Each group gets a new name.")

        # Options stay the full cluster list so keyed widgets never hold a dropped
        # value. Widget keys are the source of truth; merge_groups mirrors them.
        n_rows = len(st.session_state.merge_groups)
        updated_groups = []

        for i in range(n_rows):
            with st.container():
                c1, c2 = st.columns([1, 2])
                with c1:
                    target_val = st.text_input(
                        "New name",
                        key=f"merge_target_{i}",
                        placeholder=f"e.g. {i + 1}",
                        label_visibility="visible",
                    )
                with c2:
                    sources_val = st.multiselect(
                        "Source clusters",
                        options=all_clusters,
                        key=f"merge_sources_{i}",
                        label_visibility="visible",
                    )
            updated_groups.append({"target": target_val, "sources": sources_val})

        st.session_state.merge_groups = updated_groups

        btn_col1, btn_col2 = st.columns([1, 1])
        with btn_col1:
            if st.button("+ Add group", key="add_group"):
                new_idx = len(st.session_state.merge_groups)
                # Drop any stale widget state a prior removal may have left here
                st.session_state.pop(f"merge_target_{new_idx}", None)
                st.session_state.pop(f"merge_sources_{new_idx}", None)
                st.session_state.merge_groups.append({"target": "", "sources": []})
                st.rerun()
        with btn_col2:
            if len(st.session_state.merge_groups) > 1:
                if st.button("Remove last", key="remove_group"):
                    last_idx = len(st.session_state.merge_groups) - 1
                    st.session_state.pop(f"merge_target_{last_idx}", None)
                    st.session_state.pop(f"merge_sources_{last_idx}", None)
                    st.session_state.merge_groups.pop()
                    st.rerun()

        st.markdown("---")

        # Live preview panel
        _render_merge_preview(st.session_state.merge_groups)

    # ── Phase 3: Apply merge ──────────────────────────────────────────────────

    st.markdown("---")

    valid_groups = [
        g for g in st.session_state.merge_groups
        if g.get("target", "").strip() and g.get("sources")
    ]

    if not valid_groups:
        st.info("Define at least one merge group above, then click **Finalize Merging**.")
    else:
        errors = _validate_merge_groups(st.session_state.merge_groups, all_clusters)
        if errors:
            for e in errors:
                st.error(e)

        # Find unchanged clusters
        remapped_sources: set[str] = set()
        for g in valid_groups:
            remapped_sources.update(s.strip() for s in g["sources"])
        passthrough = [c for c in all_clusters if c not in remapped_sources]

        if passthrough:
            st.caption(
                f"**{len(passthrough)} cluster(s) will pass through unchanged** "
                f"(keeping their original Leiden number): "
                + ", ".join(f"**{c}**" for c in passthrough)
                + ". Only the source clusters listed in your merge groups will be remapped."
            )

        apply_disabled = bool(errors)

        if st.button("Finalize Merging", type="primary", disabled=apply_disabled):
            with st.spinner("Merging clusters and packaging output…"):
                cluster_mapping: dict[str, str] = {}
                for g in valid_groups:
                    t = g["target"].strip()
                    for s in g["sources"]:
                        cluster_mapping[s.strip()] = t

                leiden_copy = adata.obs['leiden'].astype(str).copy()
                leiden_copy.replace(cluster_mapping, inplace=True)
                # Keep leiden categorical so scanpy's palette path stays valid
                adata.obs['leiden'] = leiden_copy.astype('category')

                random_number = random.randint(1000, 9999)
                merge_log_lines = ["Cluster Merging Input:", _build_format_string(valid_groups), ""]
                for g in valid_groups:
                    t = g["target"].strip()
                    sources = ", ".join(s.strip() for s in g["sources"])
                    merge_log_lines.append(f"  {t} <- {sources}")
                merge_log = "\n".join(merge_log_lines)

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zf:
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
                        tmp_path = tmp.name
                    adata.write_h5ad(tmp_path)
                    zf.write(tmp_path, arcname=f"merged_adata_{random_number}.h5ad")
                    os.unlink(tmp_path)

                    with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode='w') as log_tmp:
                        log_tmp.write(merge_log)
                        log_tmp_path = log_tmp.name
                    zf.write(log_tmp_path, arcname=f"merge_log_{random_number}.txt")
                    os.unlink(log_tmp_path)

            # Persist the result so the download button survives later reruns, then reset the builder.
            st.session_state["merge_result"] = {
                "zip": zip_buffer.getvalue(),
                "filename": f"merged_data_{random_number}.zip",
                "n_groups": len(valid_groups),
                "n_sources": len(cluster_mapping),
            }
            _reset_merge_builder()
            st.rerun()

    # ── Persisted merge result (survives download-click reruns) ────────────────

    result = st.session_state.get("merge_result")
    if result:
        st.markdown("---")
        st.success(
            f"Merge applied — {result['n_groups']} group(s), "
            f"{result['n_sources']} source clusters remapped. "
            f"The cluster map above reflects the update."
        )
        st.download_button(
            label=f"Download merged data + log  ({result['filename']})",
            data=result["zip"],
            file_name=result["filename"],
            mime="application/zip",
        )
        if st.button("Start a new merge", key="clear_merge_result"):
            st.session_state.pop("merge_result", None)
            st.rerun()
