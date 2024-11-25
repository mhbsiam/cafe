#Streamlit
import streamlit as st

#System
import os
import time
import gc
import io
from io import BytesIO
import random
import zipfile

#Data handling
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pv
import anndata as ad
import scipy.sparse as sp
from scipy.sparse import issparse
import scanpy as sc

#Plotting
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

#Statistics
from scipy import stats
from scipy.stats import (
    ttest_ind,
    mannwhitneyu,
    brunnermunzel,
    shapiro,
    gaussian_kde,
    chi2_contingency,
    f_oneway,
    kruskal
)
import scipy.cluster.hierarchy as sch
from cliffs_delta import cliffs_delta

# Matplotlib
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.colors import to_rgba, to_hex
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter, FuncFormatter
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import cmasher as cmr

plt.rcParams['font.family'] = 'Arial'
pio.templates["custom"] = pio.templates["plotly"]
pio.templates["custom"].layout.font.family = "Arial"
pio.templates.default = "custom"

#result_dir = os.path.join(os.getcwd(), 'result')
#os.makedirs(result_dir, exist_ok=True)

image_path = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path)

# Display an image
image_path = os.path.join('bin', 'img', 'logo_v2.png')
st.image(image_path, caption='', use_column_width=True)


st.title("Visualization")
st.write('*The app includes visualization options such as UMAP plots, marker expression plots, heatmaps, dendogram, dotplots, and barplots with user-adjustable settings like plot dimensions, dot size, and colormap selection. It enables users to compare markers, clusters and saves the figures to the current directory/result folder.*')

option = st.radio(
    "Choose your option:",
    ["None", "Load AnnData"], help="Reload the page to start from the beginning"
)

if option == "Load AnnData":
    st.write("You have selected: Load AnnData")
    uploaded_anndata = st.file_uploader("Once you have the AnnData, upload here", type=["h5ad"])

    if uploaded_anndata is not None:
        adata = ad.read_h5ad(uploaded_anndata)
        st.write(f"Loaded AnnData file: {uploaded_anndata.name}")

        total_samples = len(adata.obs['SampleID'].unique())
        total_cells = adata.n_obs
        total_leiden_clusters = len(adata.obs['leiden'].unique())

        col1, col2, col3 = st.columns(3)
        col1.metric("Total Samples", total_samples)
        col2.metric("Total Cells", total_cells)
        col3.metric("Leiden Clusters", total_leiden_clusters)

        st.write(adata)
        datacount = adata.obs["SampleID"].value_counts()
        st.write("Cell counts per sample")
        st.table(datacount)

        group_column = 'Group'
        cluster_column = 'leiden'

        unique_groups = adata.obs[group_column].unique()
        unique_clusters = sorted(adata.obs[cluster_column].unique(), key=int)

        st.subheader("Summary Tables")

        with st.form("data_form"):

            expression_method = st.radio("Select Expression Calculation Method", ('Mean', 'Median'))

            submit_button = st.form_submit_button(label='Generate Tables')

            if submit_button:
                cluster_list = []
                num_cells_list = []
                group_count_dict = {group: [] for group in unique_groups}

                for cluster in unique_clusters:
                    cluster_data = adata[adata.obs[cluster_column] == cluster]
                    num_cells = cluster_data.shape[0]
                    cluster_list.append(cluster)
                    num_cells_list.append(num_cells)

                    for group in unique_groups:
                        group_cells = cluster_data.obs[cluster_data.obs[group_column] == group].shape[0]
                        group_count_dict[group].append(group_cells)

                cluster_table = pd.DataFrame({
                    'Cluster': cluster_list,
                    'Number of Cells': num_cells_list
                })

                for group in unique_groups:
                    cluster_table[f'Cells in {group}'] = group_count_dict[group]

                st.write("Cluster and Group Cell Counts")
                st.dataframe(cluster_table)

                def calculate_expression(data, method):
                    if method == 'Mean':
                        return data.X.mean()
                    elif method == 'Median':
                        return np.median(data.X)

                marker_names = adata.var_names
                expression_dict = {'Marker': marker_names}

                for group in unique_groups:
                    group_expressions = [
                        calculate_expression(adata[adata.obs[group_column] == group][:, marker], expression_method)
                        for marker in marker_names
                    ]
                    expression_dict[f'{expression_method} Expression in {group}'] = group_expressions

                expression_group_table = pd.DataFrame(expression_dict)
                st.write(f"Marker {expression_method} Expression by Group")
                st.dataframe(expression_group_table)

                expression_dict = {'Marker': marker_names}

                for cluster in unique_clusters:
                    cluster_expressions = [
                        calculate_expression(adata[adata.obs[cluster_column] == cluster][:, marker], expression_method)
                        for marker in marker_names
                    ]
                    expression_dict[f'Cluster {cluster}'] = cluster_expressions

                expression_cluster_table = pd.DataFrame(expression_dict)
                st.write(f"Marker {expression_method} Expression by Cluster")
                st.dataframe(expression_cluster_table)

        st.subheader("Overview")

        with st.form(key='All_Markers_Ranking'):
            expression_method = st.radio("Select Expression Calculation Method", ('Mean', 'Median'))
            test_selected = st.radio('Select a statistical test:', ('Parametric (T-test)', 'Non-parametric (Wilcoxon)'))
            submit_button = st.form_submit_button(label='Apply')

        if submit_button:
            method = 't-test' if test_selected == 'Parametric (T-test)' else 'wilcoxon'

            sc.tl.rank_genes_groups(adata, 'leiden', method=method)
            sc.tl.dendrogram(adata, groupby='leiden')

            fig, ax = plt.subplots(figsize=(15, 8))
            sc.pl.dendrogram(adata, groupby='leiden', ax=ax, show=False)
            st.pyplot(fig)

            svg_buffer = BytesIO()
            fig.savefig(svg_buffer, format="svg")
            svg_buffer.seek(0)
            st.download_button(label="Download Dendrogram", data=svg_buffer, file_name="dendrogram.svg", mime="image/svg+xml")
            plt.close(fig)

            fig1, ax1 = plt.subplots(figsize=(10, 7))
            sc.pl.umap(adata, color='leiden', legend_fontsize=8, size=5, show=False, palette='tab20c', ax=ax1, legend_loc='on data')
            ax1.axis('off')
            st.pyplot(fig1)
            plt.close(fig1)

            def calexpression(data, method):
                if method == 'Mean':
                    return np.mean(data)
                elif method == 'Median':
                    return np.median(data)

            cluster_marker_dict = {}

            unique_clusters_1 = adata.obs['leiden'].unique()

            for xyz in unique_clusters_1:
                cluster_mask = adata.obs['leiden'] == xyz
                cluster_data = adata[cluster_mask]

                marker_expressions = []
                markers = adata.var_names
                for marker in markers:
                    cluster_expression = cluster_data[:, marker].X.toarray().flatten()
                    marker_expression_value = calexpression(cluster_expression, expression_method)
                    marker_expressions.append(marker_expression_value)

                sorted_idx = np.argsort(marker_expressions)[::-1]
                sorted_markers = []
                for i in sorted_idx:
                    marker_name = markers[i]
                    if marker_expressions[i] > 1000:
                        marker_name += "+"
                    elif marker_expressions[i] < 0:
                        marker_name += "-"
                    sorted_markers.append(marker_name)

                cluster_marker_dict[xyz] = sorted_markers

            sorted_cluster_keys = sorted(cluster_marker_dict.keys(), key=int)
            max_markers = max(len(markers) for markers in cluster_marker_dict.values())

            df = pd.DataFrame({cluster: cluster_marker_dict[cluster] + [''] * (max_markers - len(cluster_marker_dict[cluster]))
                                for cluster in sorted_cluster_keys}).T

            df = df.reset_index()
            df.columns = ['Cluster'] + [f'Marker {i+1}' for i in range(max_markers)]
            df['Cluster'] = 'Cluster ' + df['Cluster'].astype(str)

            st.dataframe(df)
            st.warning(f"{expression_method} expression higher than 10^3 was coded as '+' while values <0 were coded as '-'. Interpret the result with caution.")

            def calculate_expression(matrix, method):
                if method == 'Mean':
                    return np.mean(matrix, axis=0)
                elif method == 'Median':
                    return np.median(matrix, axis=0)

            clusters = sorted(adata.obs['leiden'].unique(), key=int)
            marker_names = adata.uns['rank_genes_groups']['names'].dtype.names
            marker_list = list(dict.fromkeys([m for cluster in marker_names for m in adata.uns['rank_genes_groups']['names'][cluster]]))

            expression_data = pd.DataFrame(index=marker_list, columns=clusters)
            for cluster in clusters:
                cluster_data = adata[adata.obs['leiden'] == cluster]
                cluster_expr = cluster_data[:, marker_list].X.toarray() if issparse(cluster_data.X) else cluster_data[:, marker_list].X
                expression_data.loc[:, cluster] = calculate_expression(cluster_expr, expression_method)

            expression_data = expression_data.astype(float)

            scaled_expression = expression_data.apply(lambda x: (x - x.mean()) / x.std() if x.std() != 0 else x - x.mean(), axis=1).fillna(0)

            gene_linkage = sch.linkage(scaled_expression, method='average')
            ordered_genes = [scaled_expression.index[i] for i in sch.dendrogram(gene_linkage, no_plot=True)['leaves']]

            fig = px.imshow(
                scaled_expression.loc[ordered_genes, clusters],
                labels=dict(x="Cluster", y="Marker", color="Scaled Expression"),
                x=[f"Cluster {c}" for c in clusters], y=ordered_genes,
                color_continuous_scale='RdBu_r', aspect='auto'
            )
            fig.update_xaxes(side="bottom", tickangle=45)
            fig.update_layout(width=1200, height=800, margin=dict(l=100, r=20, t=50, b=100))

            st.plotly_chart(fig)
            st.warning(f"Matrix plot was rendered based on {expression_method} expression of markers. Expression values were scaled for visualization.")


            # Dot Plot
            fig2, ax2 = plt.subplots(figsize=(10, 8))
            sc.pl.dotplot(
                adata,
                var_names=ordered_genes,
                groupby='leiden',
                dendrogram=True,
                cmap='Reds',
                standard_scale='var',
                show=False,
                ax=ax2
            )
            plt.tight_layout()
            st.pyplot(fig2)
            st.warning(f"Dot plot was rendered based on Mean expression of markers.")

            # Download Button for Dot Plot
            svg_buffer = BytesIO()
            fig2.savefig(svg_buffer, format="svg")
            svg_buffer.seek(0)
            st.download_button(label="Download Dot Plot", data=svg_buffer, file_name="dotplot.svg", mime="image/svg+xml")
            plt.close(fig2)

        st.divider()

        colors_tab20b = list(plt.get_cmap('tab20b').colors)
        colors_tab20c = list(plt.get_cmap('tab20c').colors)
        colors_combined = colors_tab20b + colors_tab20c  # Total of 40 colors

        new_cmap = ListedColormap(colors_combined)

        st.subheader("UMAP Analysis:")

        tab1, tab2, tab3, tab4, tab5 = st.tabs(["Plot1", "Plot2", "Plot3", "Plot4", "Plot5"])

        with tab1:

            # Combined UMAP Plots
            with st.form(key='Plot1'):
                st.write("Combined UMAP Plots:")

                umap_width = st.slider("Select plot width", min_value=5, max_value=20, value=8, key="umap_width")
                umap_height = st.slider("Select plot height", min_value=5, max_value=20, value=6, key="umap_height")
                dot_size = st.slider("Select dot size", min_value=1, max_value=100, value=5, key="dot_size")

                cmap_option_20 = ['Set2', 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                                  'Set1', 'Set3', 'GnBu', 'tab10', 'tab20', 'tab20b', 'tab20c', 'tab40bc', 'RdBu_r']
                umap_cmap = st.selectbox("Select colormap", options=cmap_option_20,
                                         index=cmap_option_20.index('tab20c'), key="umap_cmap")

                umap_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG', 'JPEG'), key='UMAPfileformat')

                submit_umap = st.form_submit_button("Plot UMAPs")

            if submit_umap:
                def get_palette(cmap_input, num_categories):
                    cmap = plt.get_cmap(cmap_input)
                    if hasattr(cmap, 'colors') and len(cmap.colors) >= num_categories:
                        return cmap.colors[:num_categories]
                    else:
                        return [cmap(i / num_categories) for i in range(num_categories)]

                color_by = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

                random_number = random.randint(1000, 9999)
                zip_buffer = io.BytesIO()

                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    # UMAP by 'leiden'
                    fig1, ax1 = plt.subplots(figsize=(umap_width, umap_height))
                    num_leiden_clusters = adata.obs['leiden'].nunique()
                    palette_leiden = get_palette(umap_cmap, num_leiden_clusters)
                    sc.pl.umap(adata, color=color_by, legend_fontsize=8, size=dot_size,
                               show=False, title='', palette=palette_leiden, ax=ax1, legend_loc='on data')
                    ax1.axis('off')
                    st.pyplot(fig1)

                    fig_data = io.BytesIO()
                    fig1.savefig(fig_data, format=umap_file_format.lower(), dpi=300, bbox_inches='tight')
                    fig_data.seek(0)
                    st.download_button(
                        label=f"Download UMAP Leiden ({umap_file_format.upper()})",
                        data=fig_data,
                        file_name=f"UMAP_leiden_{random_number}.{umap_file_format.lower()}",
                        mime=f"image/{umap_file_format.lower()}",
                    )
                    zip_file.writestr(f"UMAP_leiden_{random_number}.{umap_file_format.lower()}", fig_data.getvalue())
                    plt.close(fig1)

                    # UMAP by 'Group'
                    fig2, ax2 = plt.subplots(figsize=(umap_width, umap_height))
                    num_groups = adata.obs['Group'].nunique()
                    palette_group = get_palette(umap_cmap, num_groups)
                    sc.pl.umap(adata, color='Group', legend_fontsize=8, size=dot_size,
                               show=False, title='', palette=palette_group, ax=ax2)
                    ax2.axis('off')
                    st.pyplot(fig2)

                    fig_data = io.BytesIO()
                    fig2.savefig(fig_data, format=umap_file_format.lower(), dpi=300, bbox_inches='tight')
                    fig_data.seek(0)
                    st.download_button(
                        label=f"Download UMAP Group ({umap_file_format.upper()})",
                        data=fig_data,
                        file_name=f"UMAP_group_{random_number}.{umap_file_format.lower()}",
                        mime=f"image/{umap_file_format.lower()}",
                    )
                    zip_file.writestr(f"UMAP_group_{random_number}.{umap_file_format.lower()}", fig_data.getvalue())
                    plt.close(fig2)

                    # Individual UMAPs for each group
                    unique_groups = adata.obs['Group'].unique()
                    all_leiden_clusters = adata.obs['leiden'].unique()
                    num_all_leiden_clusters = len(all_leiden_clusters)
                    consistent_palette = get_palette(umap_cmap, num_all_leiden_clusters)
                    leiden_color_map = dict(zip(all_leiden_clusters, consistent_palette))

                    for group_id in unique_groups:
                        adata_group = adata[adata.obs['Group'] == group_id]
                        group_palette = [leiden_color_map[leiden] for leiden in adata_group.obs['leiden'].cat.categories]

                        fig, ax = plt.subplots(figsize=(umap_width, umap_height))
                        sc.pl.umap(
                            adata_group,
                            color='leiden',
                            palette=group_palette,
                            legend_fontsize=12,
                            size=dot_size,
                            show=False,
                            ax=ax,
                            title=f"{group_id}",
                            legend_loc='none'
                        )
                        ax.set_xlabel('')
                        ax.set_ylabel('')
                        st.pyplot(fig)

                        fig_data = io.BytesIO()
                        fig.savefig(fig_data, format=umap_file_format.lower(), dpi=300, bbox_inches='tight')
                        fig_data.seek(0)
                        st.download_button(
                            label=f"Download UMAP {group_id} ({umap_file_format.upper()})",
                            data=fig_data,
                            file_name=f"UMAP_{group_id}_{random_number}.{umap_file_format.lower()}",
                            mime=f"image/{umap_file_format.lower()}",
                        )
                        zip_file.writestr(f"UMAP_{group_id}_{random_number}.{umap_file_format.lower()}", fig_data.getvalue())
                        plt.close(fig)
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All UMAP Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"UMAP_plots_{random_number}.zip",
                    mime="application/zip",
                )


        def set_colorbar_format(colorbar, format_type):
            formatter = ScalarFormatter(useMathText=True)
            if format_type == 'Exact values':
                formatter.set_scientific(False)
                formatter.set_useOffset(False)
            else:
                formatter.set_scientific(True)
            colorbar.formatter = formatter
            colorbar.update_ticks()

        with tab2:

            with st.form(key='Plot2'):
                st.write("Marker Expression UMAPs:")

                split_by_group = st.radio("Split by group?", ('No', 'Yes'), index=0, key='split_by_group')

                key_suffix = '_split' if split_by_group == 'Yes' else '_no_split'

                marker_expr_width = st.slider("Select plot width", min_value=5, max_value=20, value=18, key=f"marker_expr_width{key_suffix}")
                marker_expr_height = st.slider("Select plot height", min_value=5, max_value=20, value=14, key=f"marker_expr_height{key_suffix}")

                marker_dot_size = st.slider("Select dot size", min_value=0.01, max_value=20.00, value=2.00, key=f"marker_dot_size{key_suffix}")

                vmin_slider = st.slider("Select vmin percentile", min_value=0, max_value=100, value=0, key=f"vmin_slider{key_suffix}", help="Controls the lower percentile of the color scale. A lower vmin focuses the color map on lower values.")
                vmax_slider = st.slider("Select vmax percentile", min_value=0, max_value=100, value=100, key=f"vmax_slider{key_suffix}", help="Controls the upper percentile of the color scale. A higher vmax expands the color map towards higher values.")

                show_colorbar = st.radio("Show colorbars for each plot?", ('Yes', 'No'), index=0, key=f'show_colorbar{key_suffix}')

                num_of_col = st.slider("Select number of columns in plot", min_value=1, max_value=7, value=4, key=f'num_of_col{key_suffix}')

                cmap_option = st.selectbox(
                    "Select a color map",
                    ['Gray to Red', 'Gray to light green', 'Gray to teal', 'Gray to light blue', 'Gray to cyan',
                    'viridis_r', 'plasma_r', 'inferno_r', 'magma_r', 'cividis_r',
                    'Rainforest_r', 'Ocean_r', 'Swamp_r', 'Torch_r', 'Arctic_r', 'Jungle_r', 'Amber_r', 'Eclipse_r', 'Sepia_r',
                    'Savanna_r', 'Neutral_r', 'Sapphire_r', 'Watermelon', 'Waterlily', 'Holly', 'Seasons', 'Fusion', 'Viola'],
                    key=f"cmap_option{key_suffix}"
                )

                available_markers = adata.var_names.tolist()  # Extract marker names from AnnData object
                selected_markers = st.multiselect("Select markers to plot", available_markers, default=available_markers[:5], key=f'selected_markers{key_suffix}')  # Default to first 5 markers

                number_format = 'Exact values'
                marker_expr_file_format = st.radio("Select file format to save the Marker Expression figure", ('PNG', 'PDF', 'SVG', 'JPEG'), key=f'marker_expr_file_format{key_suffix}')
                submit_marker_expression = st.form_submit_button("Plot Marker Expression UMAPs")

            if submit_marker_expression:

                random_number = random.randint(1000, 9999)
                markers_to_plot = [marker for marker in selected_markers if marker in available_markers]

                if not markers_to_plot:
                    st.error("None of the selected markers are found in the AnnData object.")
                else:
                    if cmap_option == 'Gray to Red':
                        colors = [(0.8, 0.8, 0.8), (1, 0, 0)]  # Gray to red
                        cm = LinearSegmentedColormap.from_list('Gray to Red', colors, N=100)
                    elif cmap_option == 'Gray to teal':
                        colors = [(0.8, 0.8, 0.8), (0.18, 0.55, 0.34)]  # Gray to teal
                        cm = LinearSegmentedColormap.from_list('Gray to teal', colors, N=100)
                    elif cmap_option == 'Gray to light green':
                        colors = [(0.8, 0.8, 0.8), (0.282, 0.561, 0.192)]  # Gray to light green
                        cm = LinearSegmentedColormap.from_list('Gray to light green', colors, N=100)
                    elif cmap_option == 'Gray to light blue':
                        colors = [(0.8, 0.8, 0.8), (0, 0.298, 0.427)]  # Gray to light blue
                        cm = LinearSegmentedColormap.from_list('Gray to light blue', colors, N=100)
                    elif cmap_option == 'Gray to cyan':
                        colors = [(0.8, 0.8, 0.8), (0, 0.227, 0.427)]  # Gray to cyan
                        cm = LinearSegmentedColormap.from_list('Gray to cyan', colors, N=100)
                    elif cmap_option in ['Rainforest_r', 'Ocean_r', 'Swamp_r', 'Torch_r', 'Arctic_r',
                                         'Jungle_r', 'Amber_r', 'Eclipse_r', 'Sepia_r', 'Savanna_r',
                                         'Neutral_r', 'Sapphire_r', 'Watermelon', 'Waterlily', 'Holly', 'Seasons', 'Fusion', 'Viola']:
                        cm = getattr(cmr, cmap_option.lower())
                    else:
                        cm = plt.get_cmap(cmap_option)

                    if split_by_group == 'No':
                        n_cols = num_of_col
                        n_rows = int(np.ceil(len(markers_to_plot) / n_cols))

                        fig, axes = plt.subplots(n_rows, n_cols, figsize=(marker_expr_width, marker_expr_height))
                        axes = axes.flatten()

                        for i, marker in enumerate(markers_to_plot):
                            ax = axes[i]
                            marker_data = adata[:, marker].X.flatten()
                            vmin = np.percentile(marker_data, vmin_slider)
                            vmax = np.percentile(marker_data, vmax_slider)

                            if show_colorbar == 'Yes':
                                sc.pl.umap(adata, color=marker, cmap=cm, vmin=vmin, vmax=vmax, size=marker_dot_size,
                                           show=False, title=marker, ax=ax, colorbar_loc='right')
                            else:
                                sc.pl.umap(adata, color=marker, cmap=cm, vmin=vmin, vmax=vmax, size=marker_dot_size,
                                           show=False, title=marker, ax=ax, colorbar_loc=None)

                            ax.set_xticks([])
                            ax.set_yticks([])
                            ax.set_xlabel('')
                            ax.set_ylabel('')
                            ax.set_title(marker, fontsize=10)

                        for j in range(len(markers_to_plot), n_rows * n_cols):
                            fig.delaxes(axes[j])

                        plt.tight_layout()
                        st.pyplot(fig)

                        fig_data = io.BytesIO()
                        fig.savefig(fig_data, format=marker_expr_file_format.lower(), dpi=300, bbox_inches='tight')
                        fig_data.seek(0)
                        plt.close(fig)

                        st.download_button(
                            label=f"Download Marker Expression UMAP ({marker_expr_file_format.upper()})",
                            data=fig_data,
                            file_name=f"UMAP_MarkerExpression_{random_number}.{marker_expr_file_format.lower()}",
                            mime=f"image/{marker_expr_file_format.lower()}",
                        )

                    else:
                        group_column = 'Group'
                        if group_column not in adata.obs.columns:
                            st.error(f"The group column '{group_column}' is not found in the AnnData object.")
                        else:
                            groups = adata.obs[group_column].unique()

                            progress = st.progress(0)
                            for idx, group in enumerate(groups):
                                st.write(f"Plotting UMAPs for group: {group}")

                                adata_group = adata[adata.obs[group_column] == group]
                                n_cols = num_of_col
                                n_rows = int(np.ceil(len(markers_to_plot) / n_cols))

                                fig, axes = plt.subplots(n_rows, n_cols, figsize=(marker_expr_width, marker_expr_height))
                                axes = axes.flatten()

                                for i, marker in enumerate(markers_to_plot):
                                    ax = axes[i]
                                    marker_data = adata[:, marker].X.flatten()
                                    vmin = np.percentile(marker_data, vmin_slider)
                                    vmax = np.percentile(marker_data, vmax_slider)

                                    if show_colorbar == 'Yes':
                                        sc.pl.umap(adata_group, color=marker, cmap=cm, vmin=vmin, vmax=vmax, size=marker_dot_size,
                                                   show=False, title=f"{marker} ({group})", ax=ax, colorbar_loc='right')
                                    else:
                                        sc.pl.umap(adata_group, color=marker, cmap=cm, vmin=vmin, vmax=vmax, size=marker_dot_size,
                                                   show=False, title=f"{marker} ({group})", ax=ax, colorbar_loc=None)

                                    ax.set_xticks([])
                                    ax.set_yticks([])
                                    ax.set_xlabel('')
                                    ax.set_ylabel('')
                                    ax.set_title(f"{marker} ({group})", fontsize=10)

                                for j in range(len(markers_to_plot), n_rows * n_cols):
                                    fig.delaxes(axes[j])

                                plt.tight_layout()
                                st.pyplot(fig)

                                fig_data = io.BytesIO()
                                fig.savefig(fig_data, format=marker_expr_file_format.lower(), dpi=300, bbox_inches='tight')
                                fig_data.seek(0)
                                plt.close(fig)

                                st.download_button(
                                    label=f"Download Marker Expression UMAP for {group} ({marker_expr_file_format.upper()})",
                                    data=fig_data,
                                    file_name=f"UMAP_MarkerExpression_{group}_{random_number}.{marker_expr_file_format.lower()}",
                                    mime=f"image/{marker_expr_file_format.lower()}",
                                )

                                progress.progress((idx + 1) / len(groups))


        def get_colormap(cmap_option):
            custom_cmaps = {
                'Gray to Red': [(0.8, 0.8, 0.8), (1, 0, 0)],
                'Gray to teal': [(0.8, 0.8, 0.8), (0.18, 0.55, 0.34)],
                'Gray to light green': [(0.8, 0.8, 0.8), (0.282, 0.561, 0.192)],
                'Gray to light blue': [(0.8, 0.8, 0.8), (0, 0.298, 0.427)],
                'Gray to cyan': [(0.8, 0.8, 0.8), (0, 0.227, 0.427)]
            }
            if cmap_option in custom_cmaps:
                return LinearSegmentedColormap.from_list(cmap_option, custom_cmaps[cmap_option], N=100)
            try:
                return getattr(cmr, cmap_option.lower(), plt.get_cmap(cmap_option))
            except (AttributeError, ValueError):
                st.error(f"Invalid colormap `{cmap_option}`.")
                st.stop()

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab3:

            with st.form(key='Plot3'):
                unique_clusters = sorted(adata.obs[cluster_col].unique())
                selected_clusters = st.multiselect(
                    f"Select one or more {'cell types' if cluster_col == 'cell_type' else 'clusters'}",
                    options=unique_clusters,
                    default=unique_clusters,
                    key="selected_clusters_unique"
                )

                unique_markers = adata.var_names
                selected_marker = st.selectbox(
                    "Select a marker",
                    options=unique_markers,
                    key="selected_marker_unique"
                )

                umap_file_format = st.radio(
                    "Select file format to save the figures",
                    ('PNG', 'PDF', 'SVG', 'JPEG'),
                    key='umap_file_format_unique'
                )

                cmap_option = st.selectbox(
                    "Select a color map",
                    [
                        'Gray to Red', 'Gray to light green', 'Gray to teal', 'Gray to light blue', 'Gray to cyan',
                        'viridis_r', 'plasma_r', 'inferno_r', 'magma_r', 'cividis_r', 'RdBu_r', 'Set2', 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                                            'Set1', 'Set3', 'GnBu', 'tab10', 'tab20', 'tab20b', 'tab20c'
                    ],
                    key="cmap_option_unique"
                )

                umap_width = st.slider("Select plot width (in inches)", 5, 20, 8, key="umap_width_unique")
                umap_height = st.slider("Select plot height (in inches)", 5, 20, 6, key="umap_height_unique")
                dot_size = st.slider("Select dot size", 1, 100, 5, key="dot_size_unique")

                vmin_slider = st.slider("Select lower percentile (vmin)", 0, 100, 1, key="vmin_slider_unique")
                vmax_slider = st.slider("Select upper percentile (vmax)", 0, 100, 99, key="vmax_slider_unique")

                submit_umap = st.form_submit_button("Plot UMAP")

            if submit_umap:
                cluster_mask = adata.obs[cluster_col].isin(selected_clusters)
                adata_subset = adata[cluster_mask, :]  # Subsetting the data to include only the selected clusters

                if selected_marker not in adata_subset.var_names:
                    st.error(f"The selected marker `{selected_marker}` is not present in the dataset.")
                else:
                    # Generate a random number for unique file naming
                    random_number = random.randint(1000, 9999)

                    colormap = get_colormap(cmap_option)
                    marker_data = adata_subset[:, selected_marker].X.flatten()  # Flattening sparse matrices if necessary

                    vmin = np.percentile(marker_data, vmin_slider)
                    vmax = np.percentile(marker_data, vmax_slider)

                    fig1, ax1 = plt.subplots(figsize=(umap_width, umap_height))

                    sc.pl.umap(
                        adata_subset,
                        color=selected_marker,
                        cmap=colormap,
                        vmin=vmin,
                        vmax=vmax,
                        legend_fontsize=8,
                        size=dot_size,
                        show=False,
                        title=f'{selected_marker}',
                        ax=ax1,
                        legend_loc='none'
                    )

                    for spine in ax1.spines.values():
                        spine.set_visible(False)

                    ax1.set_xlabel('')
                    ax1.set_ylabel('')
                    ax1.axis('off')

                    st.pyplot(fig1)

                    fig_data = io.BytesIO()
                    fig1.savefig(fig_data, format=umap_file_format.lower(), dpi=300, bbox_inches='tight')
                    fig_data.seek(0)
                    plt.close(fig1)

                    st.download_button(
                        label=f"Download UMAP ({umap_file_format.upper()})",
                        data=fig_data,
                        file_name=f"UMAP_{cluster_col}__marker_{selected_marker}_{random_number}.{umap_file_format.lower()}",
                        mime=f"image/{umap_file_format.lower()}",
                    )


        color_profiles = {
            'Gray to Red': [(0.8, 0.8, 0.8), (1, 0, 0)],
            'Gray to teal': [(0.8, 0.8, 0.8), (0.18, 0.55, 0.34)],
            'Gray to light green': [(0.8, 0.8, 0.8), (0.282, 0.561, 0.192)],
            'Gray to light blue': [(0.8, 0.8, 0.8), (0, 0.298, 0.427)],
            'Gray to cyan': [(0.8, 0.8, 0.8), (0, 0.227, 0.427)]
        }

        def create_colormap(cmap_name):
            return LinearSegmentedColormap.from_list('custom_cmap', color_profiles[cmap_name], N=100)

        def set_colorbar_formatting(colorbar):
            colorbar.formatter = ticker.ScalarFormatter(useOffset=False, useMathText=False)
            colorbar.formatter.set_scientific(False)  # Disable scientific notation
            colorbar.update_ticks()

        with tab4:

            with st.form(key='Plot4'):
                unique_markers = adata.var_names
                selected_marker = st.selectbox(
                    "Select a marker to plot",
                    options=unique_markers,
                    key="selected_marker_unique_24"
                )

                unique_groups = adata.obs['Group'].unique()
                selected_group1 = st.selectbox("Select first group", options=unique_groups, key="group1_unique_24")
                selected_group2 = st.selectbox("Select second group", options=unique_groups, key="group2_unique_24")

                cmap_option = st.selectbox(
                    "Select a color map",
                    ['Gray to Red', 'Gray to teal', 'Gray to light green', 'Gray to light blue', 'Gray to cyan'],
                    key="cmap_option_unique_24"
                )

                umap_width = st.slider("Select plot width (in inches)", 5, 20, 8, key="umap_width_unique_24")
                umap_height = st.slider("Select plot height (in inches)", 5, 20, 6, key="umap_height_unique_24")
                dot_size = st.slider("Select dot size", 1, 100, 5, key="dot_size_unique_24")

                submit_umap = st.form_submit_button("Plot UMAP")

            if submit_umap:
                group1_data = adata[adata.obs['Group'] == selected_group1, :]
                group2_data = adata[adata.obs['Group'] == selected_group2, :]

                combined_data = adata[(adata.obs['Group'] == selected_group1) | (adata.obs['Group'] == selected_group2), :]
                vmin = combined_data[:, selected_marker].X.min()
                vmax = combined_data[:, selected_marker].X.max()

                colormap = create_colormap(cmap_option)

                fig, axes = plt.subplots(1, 2, figsize=(umap_width * 2, umap_height))

                sc.pl.umap(
                    group1_data,
                    color=selected_marker,
                    cmap=colormap,
                    size=dot_size,
                    show=False,
                    ax=axes[0],
                    title=f'{selected_group1}: Marker {selected_marker}',
                    vmin=vmin,
                    vmax=vmax
                )
                axes[0].axis('off')

                sc.pl.umap(
                    group2_data,
                    color=selected_marker,
                    cmap=colormap,
                    size=dot_size,
                    show=False,
                    ax=axes[1],
                    title=f'{selected_group2}: Marker {selected_marker}',
                    vmin=vmin,
                    vmax=vmax
                )
                axes[1].axis('off')

                for ax in axes:
                    colorbar = ax.collections[-1].colorbar
                    set_colorbar_formatting(colorbar)

                st.pyplot(fig)

                fig_data = io.BytesIO()
                random_number = random.randint(1000, 9999)
                fig.savefig(fig_data, format="png", dpi=300, bbox_inches="tight")
                fig_data.seek(0)
                plt.close(fig)

                st.download_button(
                    label="Download Group Comparison UMAP (PNG)",
                    data=fig_data,
                    file_name=f"Group_Comparison_UMAP_{selected_group1}_vs_{selected_group2}_{selected_marker}_{random_number}.png",
                    mime="image/png",
                )


        with tab5:

            colors_tab20b = list(plt.get_cmap('tab20b').colors)
            colors_tab20c = list(plt.get_cmap('tab20c').colors)
            colors_combined = colors_tab20b + colors_tab20c

            new_cmap = ListedColormap(colors_combined)

            with st.form(key='Plot5'):
                st.write("UMAP of each sample:")

                dot_size = st.slider("Select dot size", min_value=1, max_value=100, value=15, key="sample_umap_dot_size_2")
                plot_width = st.slider("Select plot width", min_value=1, max_value=20, value=5, key="sample_umap_plot_width_2")
                plot_height = st.slider("Select plot height", min_value=1, max_value=20, value=4, key="sample_umap_plot_height_2")

                group_table_file_format = st.radio("Select file format to save the Frequency Table figure", ('PNG', 'PDF', 'SVG'), key='UMAP_by_sample_2')

                submit_frequency_table = st.form_submit_button("Plot UMAP of each sample")

            if submit_frequency_table:
                frequency_table = adata.obs[['SampleID', 'Group']].value_counts().reset_index(name='Frequency')
                updated_frequency_table_dropped = frequency_table.drop(columns=['Frequency'])
                dictionary = dict(zip(updated_frequency_table_dropped['SampleID'], updated_frequency_table_dropped['Group']))
                ordered_samples = list(dictionary.keys())
                random_number = random.randint(1000, 9999)

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for sample_id in ordered_samples:
                        adata_group = adata[adata.obs['SampleID'] == sample_id]

                        fig, ax = plt.subplots(figsize=(plot_width, plot_height))
                        sc.pl.umap(
                            adata_group,
                            color='leiden',
                            legend_fontsize=6,
                            size=dot_size,
                            show=False,
                            ax=ax,
                            title=f"{sample_id}",
                            legend_loc='none'
                        )

                        for spine in ax.spines.values():
                            spine.set_visible(False)

                        ax.set_xlabel('')
                        ax.set_ylabel('')
                        st.pyplot(fig)

                        fig_data = io.BytesIO()
                        fig.savefig(fig_data, format=group_table_file_format.lower(), dpi=300, bbox_inches='tight')
                        fig_data.seek(0)
                        plt.close(fig)

                        st.download_button(
                            label=f"Download UMAP for {sample_id} ({group_table_file_format.upper()})",
                            data=fig_data,
                            file_name=f"UMAP_{sample_id}_{random_number}.{group_table_file_format.lower()}",
                            mime=f"image/{group_table_file_format.lower()}",
                        )

                        zip_file.writestr(
                            f"UMAP_{sample_id}_{random_number}.{group_table_file_format.lower()}",
                            fig_data.getvalue()
                        )

                zip_buffer.seek(0)
                st.download_button(
                    label="Download All UMAPs as ZIP",
                    data=zip_buffer,
                    file_name=f"All_UMAPs_{random_number}.zip",
                    mime="application/zip",
                )

        st.divider()

        st.subheader("Cells and Clusters Analysis:")

        tab6, tab7, tab8, tab9 = st.tabs(["Plot6", "Plot7", "Plot8", "Plot9"])

        heatmap_cmap_options = [
            'Blues', 'plasma', 'inferno', 'magma', 'cividis', 'viridis_r',
            'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
            'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
            'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'coolwarm', 'Spectral_r'
        ]

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab6:

            with st.form(key='Plot6'):
                st.write("Heatmap Plot:")

                heatmap_width = st.slider("Select plot width", min_value=1, max_value=20, value=8, key="heatmap_width")
                heatmap_height = st.slider("Select plot height", min_value=0.1, max_value=5.0, value=3.0, key="heatmap_height")
                heatmap_cmap = st.selectbox("Select the colormap", heatmap_cmap_options, index=0, key="heatmap_cmap")

                heatmap_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG', 'JPEG'), key='heatmap_file_format')

                submit_combined_heatmap = st.form_submit_button("Plot Combined Heatmaps")

            if submit_combined_heatmap:
                random_number = random.randint(1000, 9999)

                frequency_table = adata.obs[['SampleID', 'Group']].value_counts().reset_index(name='Frequency')
                updated_frequency_table_dropped = frequency_table.drop(columns=['Frequency'])
                dictionary = dict(zip(updated_frequency_table_dropped['SampleID'], updated_frequency_table_dropped['Group']))

                combined_df = pd.DataFrame(adata.obs)
                unique_groups = list(set(dictionary.values()))
                combined_df['Group'] = pd.Categorical(combined_df['Group'], categories=unique_groups, ordered=True)

                heatmap_data = combined_df.groupby(['Group', cluster_col], observed=True).size().unstack(fill_value=0)

                plt.figure(figsize=(heatmap_width, heatmap_height))
                sns.heatmap(heatmap_data, cmap=heatmap_cmap, annot=False, fmt='d', cbar=True, linewidths=0.5)
                plt.title(f'Number of Cells in Each {cluster_col.capitalize()} by Group')
                plt.xlabel('')
                plt.ylabel('')
                plt.xticks(rotation=90)
                plt.yticks(rotation=0)
                plt.tight_layout()
                st.pyplot(plt.gcf())

                heatmap_data_buffer = io.BytesIO()
                plt.savefig(heatmap_data_buffer, format=heatmap_file_format.lower(), dpi=300, bbox_inches='tight')
                heatmap_data_buffer.seek(0)
                plt.close()

                st.download_button(
                    label=f"Download Number of Cells Heatmap ({heatmap_file_format.upper()})",
                    data=heatmap_data_buffer,
                    file_name=f"Numberofcells_{cluster_col.capitalize()}_byGroup_{random_number}.{heatmap_file_format.lower()}",
                    mime=f"image/{heatmap_file_format.lower()}",
                )

                relative_freq_data = heatmap_data.div(heatmap_data.sum(axis=1), axis=0)

                plt.figure(figsize=(heatmap_width, heatmap_height))
                sns.heatmap(relative_freq_data, cmap=heatmap_cmap, annot=False, fmt='.2f', cbar=True, linewidths=0.5)
                plt.title('Relative Frequency of Cells in Each Leiden Cluster by Group')
                plt.xlabel('')
                plt.ylabel('')
                plt.xticks(rotation=90)
                plt.yticks(rotation=0)
                plt.tight_layout()
                st.pyplot(plt.gcf())

                relative_freq_buffer = io.BytesIO()
                plt.savefig(relative_freq_buffer, format=heatmap_file_format.lower(), dpi=300, bbox_inches='tight')
                relative_freq_buffer.seek(0)
                plt.close()

                st.download_button(
                    label=f"Download Relative Frequency Heatmap ({heatmap_file_format.upper()})",
                    data=relative_freq_buffer,
                    file_name=f"RelativeFrequency_{cluster_col.capitalize()}_byGroup_{random_number}.{heatmap_file_format.lower()}",
                    mime=f"image/{heatmap_file_format.lower()}",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab7:

            with st.form(key='Plot7'):
                st.write("Combined Bar Chart Plot:")

                chart_width = st.slider("Select chart width", min_value=400, max_value=1500, value=400, key="chart_width")
                chart_height = st.slider("Select chart height", min_value=400, max_value=1000, value=600, key="chart_height")
                bar_width = st.slider("Select bar width", min_value=0.1, max_value=1.0, value=0.4, key="bar_width")
                spacing = st.slider("Select spacing", min_value=0.1, max_value=1.0, value=0.1, key="spacing")
                font_size = st.slider("Select the font size for x-y axis", min_value=10, max_value=20, value=14, key="fontsize_slider")

                color_palettes = ['Alphabet', 'Plotly', 'D3', 'G10', 'T10', 'Dark24', 'Light24', 'Magenta', 'Pastel', 'Bold', 'Safe', 'Dark2', 'Set2']
                selected_palette = st.selectbox("Select a categorical color palette for cell types:", options=color_palettes)

                bar_chart_file_format = st.radio("Select file format to save the figures", ('SVG', 'PDF'), key='bar_chart_file_format')

                submit_combined_bar_chart = st.form_submit_button("Plot Combined Bar Charts")

            if submit_combined_bar_chart:
                random_number = random.randint(1000, 9999)

                frequency_table = adata.obs[['SampleID', 'Group']].value_counts().reset_index(name='Frequency')
                updated_frequency_table_dropped = frequency_table.drop(columns=['Frequency'])
                dictionary = dict(zip(updated_frequency_table_dropped['SampleID'], updated_frequency_table_dropped['Group']))

                combined_df = pd.DataFrame(adata.obs)
                group_order = list(set(dictionary.values()))
                combined_df['Group'] = pd.Categorical(combined_df['Group'], categories=group_order, ordered=True)
                cluster_group_counts = combined_df.groupby([cluster_col, 'Group']).size().unstack(fill_value=0)
                relative_cell_percentages = cluster_group_counts.div(cluster_group_counts.sum(axis=0), axis=1) * 100
                num_clusters = len(combined_df[cluster_col].unique())
                colors = px.colors.qualitative.__getattribute__(selected_palette)

                fig_relative = go.Figure()
                for i, cluster in enumerate(relative_cell_percentages.index):
                    cluster_name = cluster.replace('cell_type ', '')
                    fig_relative.add_trace(
                        go.Bar(
                            x=group_order,
                            y=relative_cell_percentages.loc[cluster],
                            name=f'{cluster_name}',
                            marker_color=colors[i % len(colors)],
                            width=bar_width
                        )
                    )

                fig_relative.update_layout(
                    barmode='stack',
                    title=f"Relative Cell Percentages per Group ({cluster_col.capitalize()})",
                    xaxis_title="Group",
                    yaxis_title="Percentage of Cells",
                    xaxis={'categoryorder': 'total descending'},
                    font=dict(size=font_size),
                    width=chart_width,
                    height=chart_height,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.plotly_chart(fig_relative)

                fig_relative_buffer = io.BytesIO()
                fig_relative.write_image(fig_relative_buffer, format=bar_chart_file_format.lower())
                fig_relative_buffer.seek(0)

                st.download_button(
                    label=f"Download Relative Cell Percentages Bar Chart ({bar_chart_file_format.upper()})",
                    data=fig_relative_buffer,
                    file_name=f"Relative_Cell_Percentages_{cluster_col}_{random_number}.{bar_chart_file_format.lower()}",
                    mime=f"image/{bar_chart_file_format.lower()}",
                )

                fig_count = go.Figure()
                for i, cluster in enumerate(cluster_group_counts.index):
                    cluster_name = cluster.replace('cell_type ', '')
                    fig_count.add_trace(
                        go.Bar(
                            x=group_order,
                            y=cluster_group_counts.loc[cluster],
                            name=f'{cluster_name}',
                            marker_color=colors[i % len(colors)],
                            width=bar_width
                        )
                    )

                fig_count.update_layout(
                    barmode='stack',
                    title=f"Cell Counts per Group ({cluster_col.capitalize()})",
                    xaxis_title="Group",
                    yaxis_title="Number of Cells",
                    xaxis={'categoryorder': 'total descending'},
                    font=dict(size=font_size),
                    width=chart_width,
                    height=chart_height,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.plotly_chart(fig_count)

                fig_count_buffer = io.BytesIO()
                fig_count.write_image(fig_count_buffer, format=bar_chart_file_format.lower())
                fig_count_buffer.seek(0)

                st.download_button(
                    label=f"Download Cell Counts Bar Chart ({bar_chart_file_format.upper()})",
                    data=fig_count_buffer,
                    file_name=f"Cell_Counts_{cluster_col}_{random_number}.{bar_chart_file_format.lower()}",
                    mime=f"image/{bar_chart_file_format.lower()}",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    zip_file.writestr(
                        f"Relative_Cell_Percentages_{cluster_col}_{random_number}.{bar_chart_file_format.lower()}",
                        fig_relative_buffer.getvalue()
                    )
                    zip_file.writestr(
                        f"Cell_Counts_{cluster_col}_{random_number}.{bar_chart_file_format.lower()}",
                        fig_count_buffer.getvalue()
                    )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Bar Charts as ZIP",
                    data=zip_buffer,
                    file_name=f"Bar_Charts_{cluster_col}_{random_number}.zip",
                    mime="application/zip",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab8:

            with st.form(key='Plot8'):
                st.write("Total Cell Percentage and Count Plots by Group:")

                plot_width = st.slider("Select plot width", min_value=5, max_value=20, value=12, key="plot_width")
                plot_height = st.slider("Select plot height", min_value=1, max_value=15, value=3, key="plot_height")

                file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG', 'JPEG'), key='CellPlotsByGroup')

                submit_cell_percentage_count = st.form_submit_button("Plot Cell Percentage and Count")

            if submit_cell_percentage_count:
                random_number = random.randint(1000, 9999)
                frequency_table = adata.obs[['SampleID', 'Group']].value_counts().reset_index(name='Frequency')
                updated_frequency_table_dropped = frequency_table.drop(columns=['Frequency'])
                dictionary = dict(zip(updated_frequency_table_dropped['SampleID'], updated_frequency_table_dropped['Group']))
                group_order = list(set(dictionary.values()))
                combined_df = pd.DataFrame(adata.obs)

                combined_df['Group'] = pd.Categorical(combined_df['Group'], categories=group_order, ordered=True)
                cluster_group_counts = combined_df.groupby([cluster_col, 'Group']).size().unstack(fill_value=0)
                cluster_sums = cluster_group_counts.sum(axis=1)
                relative_cluster_percentages = cluster_group_counts.div(cluster_sums, axis=0) * 100
                color_indices = [0, 2, 4, 6]
                colors = plt.cm.tab20(np.linspace(0, 1, 20))[color_indices]

                fig1, ax1 = plt.subplots(figsize=(plot_width, plot_height))
                bar_width = 0.5

                for cluster_index, cluster in enumerate(relative_cluster_percentages.index):
                    bottom = 0
                    for group_index, group in enumerate(group_order):
                        height = relative_cluster_percentages.loc[cluster, group] if group in relative_cluster_percentages.columns else 0
                        ax1.bar(cluster_index, height, bottom=bottom, width=bar_width, edgecolor='white', color=colors[group_index], label=group if bottom == 0 else "")
                        bottom += height

                ax1.set_xticks(range(len(relative_cluster_percentages.index)))
                ax1.set_xticklabels(relative_cluster_percentages.index, rotation=45, ha='right')
                ax1.set_ylabel('Percentage')
                ax1.set_xlabel(cluster_col.capitalize())
                ax1.spines['top'].set_visible(False)
                ax1.spines['right'].set_visible(False)
                legend_elements = [Patch(facecolor=colors[i], label=group_order[i]) for i in range(len(group_order))]
                ax1.legend(handles=legend_elements, title='Group', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
                plt.tight_layout()
                st.pyplot(fig1)

                percentage_buffer = io.BytesIO()
                fig1.savefig(percentage_buffer, format=file_format.lower(), dpi=300, bbox_inches='tight')
                percentage_buffer.seek(0)
                plt.close(fig1)

                st.download_button(
                    label=f"Download Percentage Plot ({file_format.upper()})",
                    data=percentage_buffer,
                    file_name=f"CellPercentageByGroup_{cluster_col.capitalize()}_{random_number}.{file_format.lower()}",
                    mime=f"image/{file_format.lower()}",
                )

                fig2, ax2 = plt.subplots(figsize=(plot_width, plot_height))
                bar_width = 0.5

                for cluster_index, cluster in enumerate(cluster_group_counts.index):
                    bottom = 0
                    for group_index, group in enumerate(group_order):
                        count = cluster_group_counts.loc[cluster, group] if group in cluster_group_counts.columns else 0
                        ax2.bar(cluster_index, count, bottom=bottom, width=bar_width, edgecolor='white', color=colors[group_index], label=group if bottom == 0 else "")
                        bottom += count

                ax2.set_xticks(range(len(cluster_group_counts.index)))
                ax2.set_xticklabels(cluster_group_counts.index, rotation=45, ha='right')
                ax2.set_ylabel('Count')
                ax2.set_xlabel(cluster_col.capitalize())
                ax2.spines['top'].set_visible(False)
                ax2.spines['right'].set_visible(False)
                legend_elements = [Patch(facecolor=colors[i], label=group_order[i]) for i in range(len(group_order))]
                ax2.legend(handles=legend_elements, title='Group', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
                plt.tight_layout()
                st.pyplot(fig2)

                count_buffer = io.BytesIO()
                fig2.savefig(count_buffer, format=file_format.lower(), dpi=300, bbox_inches='tight')
                count_buffer.seek(0)
                plt.close(fig2)

                st.download_button(
                    label=f"Download Count Plot ({file_format.upper()})",
                    data=count_buffer,
                    file_name=f"CellCountByGroup_{cluster_col.capitalize()}_{random_number}.{file_format.lower()}",
                    mime=f"image/{file_format.lower()}",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    zip_file.writestr(
                        f"CellPercentageByGroup_{cluster_col.capitalize()}_{random_number}.{file_format.lower()}",
                        percentage_buffer.getvalue()
                    )
                    zip_file.writestr(
                        f"CellCountByGroup_{cluster_col.capitalize()}_{random_number}.{file_format.lower()}",
                        count_buffer.getvalue()
                    )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download Both Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"CellPlotsByGroup_{cluster_col.capitalize()}_{random_number}.zip",
                    mime="application/zip",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'
        with tab9:

            with st.form(key='Plot9'):
                st.write(f"Compare Relative Abundance of Cells in Selected {cluster_col.capitalize()}")

                selected_cluster_1 = st.selectbox(f"Select {cluster_col.capitalize()} 1", options=adata.obs[cluster_col].unique())
                selected_cluster_2 = st.selectbox(f"Select {cluster_col.capitalize()} 2", options=adata.obs[cluster_col].unique())

                submit_button = st.form_submit_button(label='Submit')

            if submit_button:
                random_number = random.randint(1000, 9999)

                overall_means = adata.X.mean(axis=0)

                cluster_1_data = adata[adata.obs[cluster_col] == selected_cluster_1, :].X
                cluster_2_data = adata[adata.obs[cluster_col] == selected_cluster_2, :].X

                total_cells_cluster_1 = cluster_1_data.shape[0]
                total_cells_cluster_2 = cluster_2_data.shape[0]

                num_cells_above_mean_1 = (cluster_1_data > overall_means).sum(axis=0)
                num_cells_above_mean_2 = (cluster_2_data > overall_means).sum(axis=0)

                rel_abundance_1 = num_cells_above_mean_1 / total_cells_cluster_1
                rel_abundance_2 = num_cells_above_mean_2 / total_cells_cluster_2

                fig_comparison = go.Figure()

                fig_comparison.add_trace(go.Bar(
                    x=adata.var_names,
                    y=rel_abundance_1,
                    name=f'{selected_cluster_1}',
                    marker_color='slategray'
                ))

                fig_comparison.add_trace(go.Bar(
                    x=adata.var_names,
                    y=rel_abundance_2,
                    name=f'{selected_cluster_2}',
                    marker_color='tomato'
                ))

                fig_comparison.update_layout(
                    title=f'Relative Abundance for All Markers Between {cluster_col.capitalize()} {selected_cluster_1} and {selected_cluster_2}',
                    xaxis_title='Marker',
                    yaxis_title='Relative Abundance',
                    barmode='group',
                    xaxis=dict(tickangle=-90),
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.plotly_chart(fig_comparison)

                fig_buffer = io.BytesIO()
                fig_comparison.write_image(fig_buffer, format="svg")
                fig_buffer.seek(0)

                st.download_button(
                    label="Download Comparison Plot (SVG)",
                    data=fig_buffer,
                    file_name=f"Relative_Cells_{cluster_col.capitalize()}_{selected_cluster_1}_and_{selected_cluster_2}_{random_number}.svg",
                    mime="image/svg+xml",
                )

                st.write(f"Total number of cells in {cluster_col.capitalize()} {selected_cluster_1}: {total_cells_cluster_1}")
                st.write(f"Total number of cells in {cluster_col.capitalize()} {selected_cluster_2}: {total_cells_cluster_2}")

        st.divider()

        st.subheader("Marker Expression Analysis:")

        tab10, tab11, tab12, tab13, tab14, tab15, tab16, tab17, tab18 = st.tabs(["Plot10", "Plot11", "Plot12", "Plot13", "Plot14", "Plot15", "Plot16", "Plot17", "Plot18"])

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab10:

            with st.form(key='Plot10'):
                st.write(f"All Markers vs All {cluster_col.capitalize()} Dot Plot:")

                dotplot_width = st.slider("Select plot width for Dot Plot", min_value=1, max_value=20, value=10, key="dotplot_width")
                dotplot_height = st.slider("Select plot height for Dot Plot", min_value=1, max_value=20, value=8, key="dotplot_height")
                show_dendrogram = st.radio("Show dendrogram?", ('Yes', 'No'), index=0, key="dendrogram_toggle")
                test_selected_1 = st.radio('Select a statistical test:', ('Parametric (T-test)', 'Non-parametric (Wilcoxon)'))
                scale_option = st.radio("Scale each marker independently?", ('Yes', 'No'), index=0, key="scale_toggle")
                font_size = st.slider("Select font size for labels", min_value=5, max_value=40, value=10, key="font_size")

                cmap_option_list = ['Reds', 'RdBu_r', 'viridis_r', 'plasma', 'inferno', 'magma', 'cividis', 'coolwarm', 'Spectral', 'YlGnBu', 'Blues']
                cmap_selected = st.selectbox("Select a colormap for the plot", options=cmap_option_list, index=cmap_option_list.index('Reds'), key="cmap_choice")

                dotplot_file_format = st.radio("Select file format to save the Custom Dot Plot figure", ('PNG', 'PDF', 'SVG', 'JPEG'), key='Dendogram_cluster_v2')

                submit_dotplot_dendrogram = st.form_submit_button("Plot Dot Plot and Dendrogram")

            if submit_dotplot_dendrogram:
                plt.rcParams.update({'font.size': font_size})

                method = 't-test' if test_selected_1 == 'Parametric (T-test)' else 'wilcoxon'
                sc.tl.rank_genes_groups(adata, cluster_col, method=method)

                if show_dendrogram == 'Yes':
                    sc.tl.dendrogram(adata, groupby=cluster_col)

                marker_list = [adata.uns['rank_genes_groups']['names'][i][0] for i in range(len(adata.uns['rank_genes_groups']['names']))]

                fig, ax = plt.subplots(figsize=(dotplot_width, dotplot_height))
                standard_scale = 'var' if scale_option == 'Yes' else None

                sc.pl.dotplot(
                    adata,
                    var_names=marker_list,
                    groupby=cluster_col,
                    dendrogram=(show_dendrogram == 'Yes'),
                    standard_scale=standard_scale,
                    cmap=cmap_selected,
                    show=False,
                    ax=ax
                )

                ax.tick_params(axis='x', labelsize=font_size)
                ax.tick_params(axis='y', labelsize=font_size)

                plt.tight_layout()
                st.pyplot(fig)

                random_number = random.randint(1000, 9999)
                dotplot_buffer = io.BytesIO()
                plt.savefig(dotplot_buffer, format=dotplot_file_format.lower(), dpi=300, bbox_inches='tight')
                dotplot_buffer.seek(0)
                plt.close(fig)

                st.download_button(
                    label=f"Download Dot Plot Dendrogram ({dotplot_file_format.upper()})",
                    data=dotplot_buffer,
                    file_name=f"Dendogram_{cluster_col.capitalize()}_{random_number}.{dotplot_file_format.lower()}",
                    mime=f"image/{dotplot_file_format.lower()}",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'
        marker_list = adata.var_names
        adata.obs[cluster_col] = adata.obs[cluster_col].astype(str)
        clusters = sorted(adata.obs[cluster_col].unique())
        expression_data = pd.DataFrame(index=marker_list, columns=clusters)

        def calculate_expression(matrix, method):
            if method == 'Mean':
                return np.mean(matrix, axis=0)
            elif method == 'Median':
                return np.median(matrix, axis=0)

        with tab11:

            with st.form(key="Plot11"):
                st.write(f"All Markers vs All {cluster_col.capitalize()} Matrixplot")

                expression_method = st.radio("Select Expression Calculation Method", options=["Mean", "Median"])

                width_per_cluster = st.slider("Width per Cluster", min_value=10, max_value=100, value=50, step=5)
                height_per_marker = st.slider("Height per Marker", min_value=10, max_value=100, value=40, step=5)
                font_size = st.slider("Font Size", min_value=6, max_value=20, value=12, step=1)
                color_options = st.selectbox("Select Color Scale", options=["Reds", "viridis_r", "RdBu_r", "cividis", "inferno", "plasma", "Blues", "Greens"])

                submit_button = st.form_submit_button(label="Generate Heatmap")

            if submit_button:
                for cluster in clusters:
                    cluster_data = adata[adata.obs[cluster_col] == cluster]
                    if cluster_data.shape[0] > 0:
                        if issparse(cluster_data.X):
                            cluster_expr = cluster_data[:, marker_list].X.toarray()
                        else:
                            cluster_expr = cluster_data[:, marker_list].X
                        cluster_expressions = calculate_expression(cluster_expr, expression_method)
                        expression_data[cluster] = cluster_expressions
                    else:
                        st.warning(f"{cluster_col.capitalize()} '{cluster}' has no data, skipping it.")

                expression_data.dropna(how='all', axis=1, inplace=True)

                if expression_data.isnull().all(axis=1).any():
                    st.warning("Some markers have no valid data across clusters.")

                def safe_standardize(x):
                    std = x.std()
                    return x - x.mean() if std == 0 else (x - x.mean()) / std

                scaled_expression = expression_data.apply(safe_standardize, axis=1).fillna(0)

                num_clusters = len(scaled_expression.columns)
                num_markers = len(scaled_expression.index)

                dynamic_width = max(400, min(num_clusters * width_per_cluster, 1600))
                dynamic_height = max(300, min(num_markers * height_per_marker, 2000))

                fig = px.imshow(
                    scaled_expression,
                    labels=dict(x=f"{cluster_col.capitalize()}", y="Marker", color="Scaled Expression"),
                    x=[f"{cluster}" for cluster in scaled_expression.columns],
                    y=scaled_expression.index,
                    color_continuous_scale=color_options,
                    aspect='auto',
                )

                fig.update_xaxes(side="bottom", tickangle=-45)
                fig.update_layout(
                    width=dynamic_width,
                    height=dynamic_height,
                    font=dict(size=font_size),
                    margin=dict(l=100, r=20, t=50, b=100),
                )

                st.plotly_chart(fig)

                random_number = random.randint(1000, 9999)
                matrixplot_buffer = io.BytesIO()
                fig.write_image(matrixplot_buffer, format="svg")
                matrixplot_buffer.seek(0)

                st.download_button(
                    label=f"Download Matrix Plot ({expression_method} Expression, SVG)",
                    data=matrixplot_buffer,
                    file_name=f"Matrixplot_{expression_method}_scaled_{cluster_col}_{random_number}.svg",
                    mime="image/svg+xml",
                )


        from scipy.stats import gaussian_kde

        color_options = [
            'teal', 'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black', 'blanchedalmond', 'blue',
            'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'cornsilk',
            'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta',
            'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray',
            'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick',
            'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow',
            'grey', 'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen',
            'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightgrey',
            'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue',
            'lightyellow', 'lime', 'limegreen', 'linen', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid',
            'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred',
            'midnightblue', 'mintcream', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange',
            'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru',
            'pink', 'plum', 'powderblue', 'purple', 'rebeccapurple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown',
            'seagreen', 'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue',
            'tan', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen'
        ]

        with tab12:

            with st.form(key='Plot12'):
                st.write(f"One Marker vs All Sample/Group/{cluster_col.capitalize()} Histogram Plot:")

                chart_width = st.slider("Select total plot width", min_value=5.0, max_value=20.0, value=10.0, step=0.1, key="kde_chart_width")
                chart_height = st.slider("Select total plot height", min_value=3.0, max_value=20.0, value=10.0, step=0.1, key="kde_chart_height")
                font_size = st.slider("Select the font size for labels", min_value=6, max_value=20, value=10, step=1, key="kde_font_size")
                plot_color = st.selectbox("Select color for KDE plots", options=color_options, key="kde_color_selector")
                num_cols = st.slider("Select number of columns for subplots", min_value=1, max_value=5, value=4, step=1, key="kde_num_cols")

                marker = st.selectbox("Select a marker to plot", options=adata.var_names, key="kde_marker_selector")
                group_by = st.selectbox("Group by", options=adata.obs.columns, key="kde_groupby_selector")
                x_axis_format = st.radio("Select x-axis value format", ('Exact Values', 'Scientific Notation'), key="x_axis_format")
                kde_plot_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG', 'JPEG'), key='kde_plot_file_format')

                submit_kde_plot = st.form_submit_button("Plot Histogram")

            if submit_kde_plot:
                unique_groups = adata.obs[group_by].unique()
                num_groups = len(unique_groups)
                num_rows = (num_groups // num_cols) + (num_groups % num_cols > 0)
                fig, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(chart_width, chart_height))
                axs = axs.flatten()

                for i, group in enumerate(unique_groups):
                    group_data = adata[adata.obs[group_by] == group, marker].X.toarray().flatten()

                    kde = gaussian_kde(group_data)
                    x_vals = np.linspace(group_data.min(), group_data.max(), 1000)
                    kde_vals = kde(x_vals) * len(group_data)

                    axs[i].plot(x_vals, kde_vals, color=plot_color, lw=2, label=f'{group}', alpha=0.7)
                    axs[i].fill_between(x_vals, kde_vals, color=plot_color, alpha=0.3)

                    axs[i].set_title(f'{group}', fontsize=font_size)
                    axs[i].set_xlabel('Fluorescence Intensity', fontsize=font_size)
                    axs[i].set_ylabel('Number of events', fontsize=font_size)
                    axs[i].tick_params(labelsize=font_size)

                    if x_axis_format == 'Exact Values':
                        axs[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
                    elif x_axis_format == 'Scientific Notation':
                        axs[i].xaxis.set_major_formatter(ticker.ScalarFormatter())
                        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(6, 6))

                    axs[i].legend(fontsize=font_size)

                for j in range(i + 1, len(axs)):
                    fig.delaxes(axs[j])

                plt.tight_layout()
                st.pyplot(fig)

                random_number = random.randint(1000, 9999)
                kde_plot_buffer = io.BytesIO()
                plt.savefig(kde_plot_buffer, format=kde_plot_file_format.lower(), dpi=300, bbox_inches='tight')
                kde_plot_buffer.seek(0)
                plt.close(fig)

                st.download_button(
                    label=f"Download KDE Plot ({kde_plot_file_format.upper()})",
                    data=kde_plot_buffer,
                    file_name=f"KDE_Plot_{marker}_{group_by}_{random_number}.{kde_plot_file_format.lower()}",
                    mime=f"image/{kde_plot_file_format.lower()}",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab13:

            with st.form(key='Plot_13'):
                st.write(f"One Marker vs Group Based Histogram Plot ({cluster_col.capitalize()})")

                chart_width = st.slider("Select total plot width", min_value=5.0, max_value=40.0, value=20.0, step=0.1, key="kde_chart_width_grp")
                chart_height = st.slider("Select total plot height", min_value=3.0, max_value=40.0, value=20.0, step=0.1, key="kde_chart_height__grp")
                font_size = st.slider("Select the font size for labels", min_value=6, max_value=20, value=10, step=1, key="kde_font_size_grp")
                plot_color = st.selectbox("Select color for KDE plots", options=['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black'], key="kde_color_selector_grp")
                num_cols = st.slider("Select number of columns for subplots", min_value=1, max_value=10, value=5, step=1, key="kde_num_cols_grp")
                marker = st.selectbox("Select a marker to plot", options=adata.var_names, key="kde_marker_selector_grp")
                group_by2 = 'Group'
                x_axis_format = st.radio("Select x-axis value format", ('Exact Values', 'Scientific Notation'), key="x_axis_format_grp")
                kde_plot_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG', 'JPEG'), key='kde_plot_file_format_grp')
                submit_kde_plot = st.form_submit_button("Plot Histogram")

            if submit_kde_plot:
                unique_groups = adata.obs[group_by2].unique()
                unique_clusters = adata.obs[cluster_col].unique()

                num_groups = len(unique_groups)
                num_clusters = len(unique_clusters)
                total_subplots = num_groups * num_clusters

                num_rows = (total_subplots // num_cols) + (total_subplots % num_cols > 0)

                fig, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(chart_width, chart_height))
                axs = axs.flatten()

                idx = 0
                plot_buffers = []

                for group in unique_groups:
                    for cluster in unique_clusters:
                        subset = adata[(adata.obs[group_by2] == group) & (adata.obs[cluster_col] == cluster), marker]
                        if subset.n_obs == 0:
                            continue

                        group_cluster_data = subset.X.toarray().flatten()

                        if len(group_cluster_data) > 1:
                            kde = gaussian_kde(group_cluster_data)
                            x_vals = np.linspace(group_cluster_data.min(), group_cluster_data.max(), 1000)
                            kde_vals = kde(x_vals) * len(group_cluster_data)
                        else:
                            x_vals = group_cluster_data
                            kde_vals = np.ones_like(group_cluster_data)

                        axs[idx].plot(x_vals, kde_vals, color=plot_color, lw=2, label=f'{group} - {cluster_col.capitalize()} {cluster}', alpha=0.7)
                        axs[idx].fill_between(x_vals, kde_vals, color=plot_color, alpha=0.3)

                        axs[idx].set_title(f'{group} - {cluster_col.capitalize()} {cluster}', fontsize=font_size)
                        axs[idx].set_xlabel('Fluorescence Intensity', fontsize=font_size)
                        axs[idx].set_ylabel('Number of events', fontsize=font_size)
                        axs[idx].tick_params(labelsize=font_size)

                        if x_axis_format == 'Exact Values':
                            axs[idx].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
                        elif x_axis_format == 'Scientific Notation':
                            axs[idx].xaxis.set_major_formatter(ticker.ScalarFormatter())
                            axs[idx].ticklabel_format(style='sci', axis='x', scilimits=(6, 6))

                        axs[idx].legend(fontsize=font_size)
                        idx += 1

                for j in range(idx, len(axs)):
                    fig.delaxes(axs[j])

                plt.tight_layout()
                st.pyplot(fig)

                random_number = random.randint(1000, 9999)
                kde_plot_buffer = io.BytesIO()
                plt.savefig(kde_plot_buffer, format=kde_plot_file_format.lower(), dpi=300, bbox_inches='tight')
                kde_plot_buffer.seek(0)
                plot_buffers.append(("KDE_Histogram", kde_plot_buffer))
                plt.close(fig)

                st.download_button(
                    label=f"Download KDE Histogram ({kde_plot_file_format.upper()})",
                    data=kde_plot_buffer,
                    file_name=f"KDE_Histogram_{marker}_{cluster_col}_{group_by2}_{random_number}.{kde_plot_file_format.lower()}",
                    mime=f"image/{kde_plot_file_format.lower()}",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{marker}_{cluster_col}_{group_by2}_{random_number}.{kde_plot_file_format.lower()}",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All KDE Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"All_KDE_Plots_{marker}_{cluster_col}_{group_by2}_{random_number}.zip",
                    mime="application/zip",
                )


        with tab14:

            with st.form(key='Plot_14'):
                st.write("One Marker vs All Sample/Group/Leiden Cluster Bar Plot:")

                chart_width = st.slider("Select total plot width", min_value=2.0, max_value=20.0, value=10.0, step=0.1, key="bar_chart_width_50")
                chart_height = st.slider("Select total plot height", min_value=3.0, max_value=20.0, value=6.0, step=0.1, key="bar_chart_height_50")
                font_size = st.slider("Select the font size for labels", min_value=6, max_value=20, value=10, step=1, key="bar_font_size_50")
                marker = st.selectbox("Select a marker to plot", options=adata.var_names, key="bar_marker_selector_50")
                group_by = st.selectbox("Group by", options=adata.obs.columns, key="bar_groupby_selector_50")
                agg_function = st.radio("Select aggregation function", options=['Mean', 'Median'], key="agg_function_selector_50")
                bar_color = st.color_picker("Pick a color for the bar plot", '#005677', key='bar_color_50')
                bar_plot_file_format = st.radio("Select file format to save the figures", ('SVG', 'PDF'), key='bar_plot_file_format_50')

                submit_bar_plot = st.form_submit_button("Plot Bar Plot")

            if submit_bar_plot:
                random_number = random.randint(1000, 9999)

                # Bar Plot
                if agg_function == 'Mean':
                    agg_expression = adata.to_df()[marker].groupby(adata.obs[group_by]).mean()
                else:
                    agg_expression = adata.to_df()[marker].groupby(adata.obs[group_by]).median()

                fig = go.Figure()

                fig.add_trace(go.Bar(
                    x=agg_expression.index,
                    y=agg_expression.values,
                    name=f'{agg_function} Expression',
                    marker_color=bar_color
                ))

                fig.update_layout(
                    title=f'{agg_function} {marker} Expression across {group_by}',
                    xaxis_title=group_by,
                    yaxis_title=f'{agg_function} Expression',
                    xaxis_tickangle=-90,
                    height=chart_height * 100,
                    width=chart_width * 100,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                )

                st.subheader(f'Marker profile based on {agg_function} expression')
                st.plotly_chart(fig)

                bar_plot_buffer = io.BytesIO()
                fig.write_image(bar_plot_buffer, format=bar_plot_file_format.lower())
                bar_plot_buffer.seek(0)

                st.download_button(
                    label=f"Download Bar Plot ({bar_plot_file_format.upper()})",
                    data=bar_plot_buffer,
                    file_name=f"BarPlot_{marker}_{group_by}_{agg_function}_{random_number}.{bar_plot_file_format.lower()}",
                    mime=f"image/{bar_plot_file_format.lower()}",
                )

                # Violin Plot
                set2_colors = px.colors.qualitative.Set2
                fig_violin = go.Figure()

                for i, group in enumerate(adata.obs[group_by].unique()):
                    fig_violin.add_trace(go.Violin(
                        y=adata.to_df()[marker][adata.obs[group_by] == group],
                        x=np.repeat(group, len(adata.to_df()[marker][adata.obs[group_by] == group])),
                        name=str(group),
                        box_visible=True,
                        meanline_visible=True,
                        line_color=set2_colors[i % len(set2_colors)],
                        fillcolor=set2_colors[i % len(set2_colors)]
                    ))

                fig_violin.update_layout(
                    title=f'Violin Plot of {marker} Expression across {group_by}',
                    xaxis_title=group_by,
                    yaxis_title=f'{marker} Expression',
                    height=chart_height * 100,
                    width=chart_width * 100,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    showlegend=False,
                    xaxis_tickangle=-90,
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.plotly_chart(fig_violin)

                violin_plot_buffer = io.BytesIO()
                fig_violin.write_image(violin_plot_buffer, format=bar_plot_file_format.lower())
                violin_plot_buffer.seek(0)

                st.download_button(
                    label=f"Download Violin Plot ({bar_plot_file_format.upper()})",
                    data=violin_plot_buffer,
                    file_name=f"ViolinPlot_{marker}_{group_by}_{random_number}.{bar_plot_file_format.lower()}",
                    mime=f"image/{bar_plot_file_format.lower()}",
                )

                # Box Plot
                fig_box = go.Figure()

                for i, group in enumerate(adata.obs[group_by].unique()):
                    fig_box.add_trace(go.Box(
                        y=adata.to_df()[marker][adata.obs[group_by] == group],
                        name=str(group),
                        marker_color=set2_colors[i % len(set2_colors)]
                    ))

                fig_box.update_layout(
                    title=f'Box (Whisker) Plot of {marker} Expression across {group_by}',
                    xaxis_title=group_by,
                    yaxis_title=f'{marker} Expression',
                    height=chart_height * 100,
                    width=chart_width * 100,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    showlegend=False,
                    xaxis_tickangle=-90,
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.plotly_chart(fig_box)

                box_plot_buffer = io.BytesIO()
                fig_box.write_image(box_plot_buffer, format=bar_plot_file_format.lower())
                box_plot_buffer.seek(0)

                st.download_button(
                    label=f"Download Box Plot ({bar_plot_file_format.upper()})",
                    data=box_plot_buffer,
                    file_name=f"BoxPlot_{marker}_{group_by}_{random_number}.{bar_plot_file_format.lower()}",
                    mime=f"image/{bar_plot_file_format.lower()}",
                )

                # ZIP file containing all plots
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    zip_file.writestr(
                        f"BarPlot_{marker}_{group_by}_{agg_function}_{random_number}.{bar_plot_file_format.lower()}",
                        bar_plot_buffer.getvalue()
                    )
                    zip_file.writestr(
                        f"ViolinPlot_{marker}_{group_by}_{random_number}.{bar_plot_file_format.lower()}",
                        violin_plot_buffer.getvalue()
                    )
                    zip_file.writestr(
                        f"BoxPlot_{marker}_{group_by}_{random_number}.{bar_plot_file_format.lower()}",
                        box_plot_buffer.getvalue()
                    )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"AllPlots_{marker}_{group_by}_{random_number}.zip",
                    mime="application/zip",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab15:

            with st.form(key='Plot_15'):
                st.write(f"One {cluster_col.capitalize()} vs All Markers Expression:")

                cluster_selected = st.selectbox(f'Select a {cluster_col}:', adata.obs[cluster_col].unique())
                test_selected = st.radio('Select a statistical test:', ('Parametric (T-test)', 'Non-parametric (Mann-Whitney U test/Wilcoxon)'))
                expression_method = st.radio('Select expression calculation method:', ('Mean', 'Median'))
                submit_button = st.form_submit_button(label='Run Analysis')

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                # Data Calculation
                selected_cluster_mask = adata.obs[cluster_col] == cluster_selected
                selected_cluster_data = adata[selected_cluster_mask]
                other_clusters_data = adata[~selected_cluster_mask]

                markers = adata.var_names
                fold_changes = []
                p_values = []
                test_statistics = []
                expr_selected = []
                expr_others = []

                def calculate_expression(data, method):
                    if method == 'Mean':
                        return np.mean(data)
                    elif method == 'Median':
                        return np.median(data)

                for marker in markers:
                    cluster_expression = selected_cluster_data[:, marker].X.toarray().flatten()
                    other_expression = other_clusters_data[:, marker].X.toarray().flatten()

                    fold_change = (
                        np.log2(calculate_expression(cluster_expression, expression_method) + 1) -
                        np.log2(calculate_expression(other_expression, expression_method) + 1)
                    )
                    fold_changes.append(fold_change)

                    if test_selected == 'Parametric (T-test)':
                        test_stat, p_value = ttest_ind(cluster_expression, other_expression, equal_var=False)
                    else:  # Mann-Whitney U test
                        test_stat, p_value = mannwhitneyu(cluster_expression, other_expression, alternative='two-sided')

                    p_values.append(p_value)
                    test_statistics.append(test_stat)
                    expr_selected.append(calculate_expression(cluster_expression, expression_method))
                    expr_others.append(calculate_expression(other_expression, expression_method))

                results_df = pd.DataFrame({
                    'Marker': markers,
                    'Fold Change': fold_changes,
                    f'{expression_method} Selected': expr_selected,
                    f'{expression_method} Others': expr_others,
                    'Test Statistic': test_statistics,
                    'P-value': p_values
                })

                # Bar Plot: Marker Expression
                fig = go.Figure()

                fig.add_trace(go.Bar(
                    x=results_df['Marker'],
                    y=results_df[f'{expression_method} Selected'],
                    name=f'Selected {cluster_col.capitalize()}',
                    marker_color='royalblue'
                ))

                fig.add_trace(go.Bar(
                    x=results_df['Marker'],
                    y=results_df[f'{expression_method} Others'],
                    name=f'Other {cluster_col.capitalize()}s',
                    marker_color='gray'
                ))

                fig.update_layout(
                    title=f'{expression_method} Marker Expression: {cluster_col.capitalize()} {cluster_selected} vs Other {cluster_col.capitalize()}s',
                    xaxis_title='Markers',
                    yaxis_title=f'{expression_method} Expression',
                    barmode='group',
                    xaxis_tickangle=-45,
                    height=600,
                    width=1000,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.subheader(f'Marker profile based on {expression_method} expression')
                st.plotly_chart(fig)

                bar_plot_buffer = io.BytesIO()
                fig.write_image(bar_plot_buffer, format="png")
                bar_plot_buffer.seek(0)
                plot_buffers.append(("Marker_Expression_Bar", bar_plot_buffer))

                st.download_button(
                    label="Download Marker Expression Bar Plot (PNG)",
                    data=bar_plot_buffer,
                    file_name=f"Marker_Expression_Bar_{cluster_selected}_{random_number}.png",
                    mime="image/png",
                )

                # Bar Plot: Fold Change
                results_df = results_df.dropna(subset=['Fold Change']).sort_values(by='Fold Change', ascending=False)
                fig2 = go.Figure()

                fig2.add_trace(go.Bar(
                    x=results_df['Marker'],
                    y=results_df['Fold Change'],
                    name='Fold Change',
                    marker_color='cadetblue'
                ))

                fig2.add_shape(
                    type='line',
                    x0=0, x1=1, y0=1, y1=1,
                    line=dict(
                        color='red',
                        width=2,
                        dash='dash',
                    ),
                    xref='paper',
                    yref='y',
                )

                fig2.update_layout(
                    title=f'Fold Change of Markers: {cluster_col.capitalize()} {cluster_selected} vs Other {cluster_col.capitalize()}s',
                    xaxis_title='Markers',
                    yaxis_title='Log2 Fold Change',
                    xaxis_tickangle=-45,
                    height=600,
                    width=1000,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.subheader('Marker profile based on fold change')
                st.plotly_chart(fig2)

                fold_change_buffer = io.BytesIO()
                fig2.write_image(fold_change_buffer, format="png")
                fold_change_buffer.seek(0)
                plot_buffers.append(("Fold_Change_Bar", fold_change_buffer))

                st.download_button(
                    label="Download Fold Change Bar Plot (PNG)",
                    data=fold_change_buffer,
                    file_name=f"Fold_Change_Bar_{cluster_selected}_{random_number}.png",
                    mime="image/png",
                )

                # Ranking Scores
                method = 't-test' if test_selected == 'Parametric (T-test)' else 'wilcoxon'
                sc.tl.rank_genes_groups(adata, cluster_col, method=method)

                markers = adata.uns['rank_genes_groups']['names'][cluster_selected]
                scores = adata.uns['rank_genes_groups']['scores'][cluster_selected]

                marker_df = pd.DataFrame({
                    'Marker': markers,
                    'scores': scores
                }).sort_values(by='scores', ascending=False)

                fig3 = go.Figure()

                fig3.add_trace(go.Bar(
                    x=marker_df['Marker'],
                    y=marker_df['scores'],
                    name=f'scores for {cluster_selected}',
                    marker_color='lightcoral'
                ))

                fig3.update_layout(
                    title=f'Markers vs scores for {cluster_col.capitalize()} {cluster_selected}',
                    xaxis_title='Markers',
                    yaxis_title='scores',
                    xaxis_tickangle=-45,
                    height=600,
                    width=1000,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.plotly_chart(fig3)

                ranking_buffer = io.BytesIO()
                fig3.write_image(ranking_buffer, format="png")
                ranking_buffer.seek(0)
                plot_buffers.append(("Ranking_Bar", ranking_buffer))

                st.download_button(
                    label="Download Ranking Scores Bar Plot (PNG)",
                    data=ranking_buffer,
                    file_name=f"Ranking_Bar_{cluster_selected}_{random_number}.png",
                    mime="image/png",
                )

                # Create ZIP File
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{cluster_selected}_{random_number}.png",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"All_Plots_{cluster_selected}_{random_number}.zip",
                    mime="application/zip",
                )


        def get_expression_by_group_and_cluster(adata, selected_marker):
            clusters = sorted(adata.obs['leiden'].unique().astype(int))
            groups = adata.obs['Group'].unique()

            expression_data = []

            for group in groups:
                for cluster in clusters:
                    cluster_group_data = adata[(adata.obs['leiden'] == str(cluster)) & (adata.obs['Group'] == group)]  # Ensure cluster is a string
                    marker_values = cluster_group_data[:, selected_marker].X.toarray().flatten()
                    expression_data.append({
                        'Group': group,
                        'Cluster': cluster,
                        'Marker Expression': marker_values
                    })

            return expression_data, groups

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab16:

            with st.form(key='Plot_16'):
                st.write(f"Grouped Violin Plot for Expression of a Selected Marker by Group and {cluster_col.capitalize()}")

                selected_marker = st.selectbox("Select a marker to plot", options=adata.var_names)
                file_format = st.selectbox("Select file format to save the plot", options=['SVG', 'PDF'])

                submit_group_violin_plot = st.form_submit_button("Generate Group Violin Plot")

            if submit_group_violin_plot:
                random_number = random.randint(1000, 9999)

                def get_expression_by_group_and_cluster(adata, marker):
                    expression_data = []
                    groups = adata.obs['Group'].unique()
                    clusters = adata.obs[cluster_col].unique()

                    for group in groups:
                        for cluster in clusters:
                            group_cluster_data = adata[(adata.obs['Group'] == group) & (adata.obs[cluster_col] == cluster), marker].X
                            if issparse(group_cluster_data):
                                group_cluster_data = group_cluster_data.toarray()
                            expression_data.append({
                                'Group': group,
                                'Cluster': cluster,
                                'Marker Expression': group_cluster_data.flatten()
                            })
                    return expression_data, groups

                expression_data, groups = get_expression_by_group_and_cluster(adata, selected_marker)

                fig = go.Figure()

                color_palette = ['blue', 'orange', 'green', 'purple', 'red']
                group_colors = {group: color_palette[i % len(color_palette)] for i, group in enumerate(groups)}

                for data in expression_data:
                    fig.add_trace(go.Violin(
                        x=[f"{data['Cluster']}"] * len(data['Marker Expression']),
                        y=data['Marker Expression'],
                        legendgroup=data['Group'],
                        scalegroup=data['Group'],
                        name=data['Group'],
                        line_color=group_colors[data['Group']],
                        box_visible=True,
                        meanline_visible=True,
                        showlegend=True if data['Cluster'] == clusters[0] else False
                    ))

                fig.update_layout(
                    title=f"Grouped Violin Plot of {selected_marker} Expression by Group and {cluster_col.capitalize()}",
                    xaxis_title=f"{cluster_col.capitalize()}",
                    yaxis_title="Marker Expression",
                    violinmode='group',
                    yaxis_zeroline=False
                )

                st.plotly_chart(fig)

                # Save the figure to a buffer
                violin_plot_buffer = io.BytesIO()
                fig.write_image(violin_plot_buffer, format=file_format.lower())
                violin_plot_buffer.seek(0)

                # Add a download button for the plot
                st.download_button(
                    label=f"Download Violin Plot ({file_format.upper()})",
                    data=violin_plot_buffer,
                    file_name=f"ViolinPlot_{selected_marker}_{cluster_col}_{random_number}.{file_format.lower()}",
                    mime=f"image/{file_format.lower()}",
                )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab17:

            with st.form(key='Plot_17'):
                st.write(f"One {cluster_col.capitalize()} vs Combination of Markers")

                font_size = st.slider("Select the font size for labels", min_value=6, max_value=20, value=10, step=1, key="bar_font_size_")
                selected_cluster = st.selectbox(f"Select a {cluster_col.capitalize()} cluster", options=adata.obs[cluster_col].unique(), key="cluster_selector_")
                selected_markers = st.multiselect("Select markers to plot", options=adata.var_names, key="multi_marker_selector_")
                expression_method = st.radio("Select expression calculation method", ('Mean', 'Median'), key="expression_method_bar_")
                bar_plot_file_format = st.radio("Select file format to save the figure", ('SVG', 'PDF'), key='bar_plot_file_format_')
                selected_cluster_color = st.color_picker(f"Pick a color for the selected {cluster_col.capitalize()} cluster", '#005677', key='selected_cluster_color_')
                others_color = st.color_picker("Pick a color for the others group", '#91e0ff', key='others_color_')

                submit_bar_plot = st.form_submit_button("Plot Bar Plot")

            if submit_bar_plot:
                random_number = random.randint(1000, 9999)

                marker_expression_df = pd.DataFrame(adata[:, selected_markers].X.toarray(), columns=selected_markers)
                marker_expression_df[cluster_col] = adata.obs[cluster_col].values

                def calculate_expression(data, method):
                    if method == 'Mean':
                        return data.mean(numeric_only=True)
                    elif method == 'Median':
                        return data.median(numeric_only=True)

                expression_selected_cluster = calculate_expression(
                    marker_expression_df[marker_expression_df[cluster_col] == selected_cluster],
                    expression_method
                )
                expression_others = calculate_expression(
                    marker_expression_df[marker_expression_df[cluster_col] != selected_cluster],
                    expression_method
                )

                bar_data = pd.DataFrame({
                    'Markers': selected_markers * 2,
                    f'{expression_method} Expression': pd.concat([expression_selected_cluster, expression_others]),
                    'Group': [f'Selected {cluster_col.capitalize()}'] * len(selected_markers) + ['Others'] * len(selected_markers)
                })

                bar_fig = px.bar(
                    bar_data, x='Markers', y=f'{expression_method} Expression', color='Group',
                    barmode='group', title=f'{expression_method} Expression of Markers',
                    labels={'Markers': 'Markers', f'{expression_method} Expression': f'{expression_method} Expression'},
                    color_discrete_map={
                        f'Selected {cluster_col.capitalize()}': selected_cluster_color,
                        'Others': others_color
                    }
                )

                bar_fig.update_layout(
                    autosize=True,
                    font=dict(size=font_size),
                    margin=dict(l=40, r=40, t=40, b=40),
                    title_x=0.5
                )

                st.plotly_chart(bar_fig)

                bar_plot_buffer = io.BytesIO()
                bar_fig.write_image(bar_plot_buffer, format=bar_plot_file_format.lower())
                bar_plot_buffer.seek(0)

                st.download_button(
                    label=f"Download Bar Plot ({bar_plot_file_format.upper()})",
                    data=bar_plot_buffer,
                    file_name=f"BarPlot_{selected_cluster}_{cluster_col}_{random_number}.{bar_plot_file_format.lower()}",
                    mime=f"image/{bar_plot_file_format.lower()}",
                )


        with tab18:

            with st.form(key='Plot_18'):
                st.write("Group vs Combination of Markers")

                font_size = st.slider("Select the font size for labels", min_value=6, max_value=20, value=10, step=1, key="bar_font_size_group_13")

                selected_markers = st.multiselect("Select markers to plot", options=adata.var_names, key="multi_marker_selector_group_13")

                expression_method = st.radio("Select expression calculation method", ('Mean', 'Median'), key="expression_method_bar_group_13")

                bar_plot_file_format = st.radio("Select file format to save the figure", ('SVG', 'PDF'), key='bar_plot_file_format_group_13')

                if 'Group' in adata.obs:
                    group_options = adata.obs['Group'].unique()
                    if len(group_options) >= 2:
                        group_1_color = st.color_picker(f"Pick color for {group_options[0]}", "#1f77b4", key='group_1_color')
                        group_2_color = st.color_picker(f"Pick color for {group_options[1]}", "#ff7f0e", key='group_2_color')
                    else:
                        st.error("You need at least two groups for this color selection.")

                submit_bar_plot_group = st.form_submit_button("Plot Bar Plot")

            if submit_bar_plot_group:
                if 'Group' in adata.obs:
                    random_number = random.randint(1000, 9999)

                    # Prepare data for the bar plot
                    marker_expression_df = pd.DataFrame(adata[:, selected_markers].X.toarray(), columns=selected_markers)
                    marker_expression_df['Group'] = adata.obs['Group'].values

                    def calculate_expression(data, method):
                        if method == 'Mean':
                            return data.groupby('Group').mean(numeric_only=True)
                        elif method == 'Median':
                            return data.groupby('Group').median(numeric_only=True)

                    expression_by_group = calculate_expression(marker_expression_df, expression_method)
                    bar_data = expression_by_group.reset_index().melt(id_vars='Group', var_name='Markers', value_name=f'{expression_method} Expression')

                    group_colors = {group_options[0]: group_1_color, group_options[1]: group_2_color}

                    # Create the bar plot
                    bar_fig = px.bar(
                        bar_data, x='Markers', y=f'{expression_method} Expression', color='Group',
                        color_discrete_map=group_colors,
                        barmode='group', title=f'{expression_method} Expression of Markers Among Groups',
                        labels={'Markers': 'Markers', f'{expression_method} Expression': f'{expression_method} Expression'}
                    )

                    bar_fig.update_layout(
                        autosize=True,
                        font=dict(size=font_size),
                        margin=dict(l=40, r=40, t=40, b=40),
                        title_x=0.5,
                        plot_bgcolor='white',
                        paper_bgcolor='white',
                        yaxis=dict(
                            gridcolor='lightgray',
                        )
                    )

                    st.plotly_chart(bar_fig)

                    # Save the plot to a buffer
                    bar_plot_buffer = io.BytesIO()
                    bar_fig.write_image(bar_plot_buffer, format="pdf")
                    bar_plot_buffer.seek(0)

                    # Download button for the bar plot
                    st.download_button(
                        label="Download Bar Plot (PDF)",
                        data=bar_plot_buffer,
                        file_name=f"BarPlot_GroupComparison_{random_number}.pdf",
                        mime="application/pdf",
                    )

        st.divider()

        st.subheader("Additional Plots")
        tab19, tab20, tab21, tab22, tab23, tab24, tab25, tab26, tab27, tab28, tab29, tab30 = st.tabs(["Plot19", "Plot20", "Plot21", "Plot22", "Plot23", "Plot24", "Plot25", "Plot26", "Plot27", "Plot28", "Plot29", "Plot30"])

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab19:

            with st.form(key='Plot_19'):
                st.write(f"Compare Markers Across Three {cluster_col.capitalize()} Clusters")

                selected_cluster_1 = st.selectbox(f"Select {cluster_col.capitalize()} Cluster 1", options=adata.obs[cluster_col].unique(), key='cluster_1')
                selected_cluster_2 = st.selectbox(f"Select {cluster_col.capitalize()} Cluster 2", options=adata.obs[cluster_col].unique(), key='cluster_2')
                selected_cluster_3 = st.selectbox(f"Select {cluster_col.capitalize()} Cluster 3", options=adata.obs[cluster_col].unique(), key='cluster_3')

                central_tendency = st.radio("Select Central Tendency Calculation Method", ["Mean", "Median"], key='central_tendency')

                submit_button = st.form_submit_button(label='Submit')

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []
                data_buffers = []

                # Central Tendency Bar Plot
                marker_names = adata.var_names
                expression_cluster_1 = []
                expression_cluster_2 = []
                expression_cluster_3 = []

                for marker in marker_names:
                    if central_tendency == "Median":
                        expression_cluster_1.append(np.median(adata[adata.obs[cluster_col] == selected_cluster_1, marker].X))
                        expression_cluster_2.append(np.median(adata[adata.obs[cluster_col] == selected_cluster_2, marker].X))
                        expression_cluster_3.append(np.median(adata[adata.obs[cluster_col] == selected_cluster_3, marker].X))
                    else:
                        expression_cluster_1.append(np.mean(adata[adata.obs[cluster_col] == selected_cluster_1, marker].X))
                        expression_cluster_2.append(np.mean(adata[adata.obs[cluster_col] == selected_cluster_2, marker].X))
                        expression_cluster_3.append(np.mean(adata[adata.obs[cluster_col] == selected_cluster_3, marker].X))

                data = pd.DataFrame({
                    'Marker': marker_names,
                    f'{selected_cluster_1}': expression_cluster_1,
                    f'{selected_cluster_2}': expression_cluster_2,
                    f'{selected_cluster_3}': expression_cluster_3
                })

                fig = go.Figure()

                fig.add_trace(go.Bar(
                    x=data['Marker'],
                    y=data[f'{selected_cluster_1}'],
                    name=f'{selected_cluster_1}',
                    marker_color='slategray'
                ))

                fig.add_trace(go.Bar(
                    x=data['Marker'],
                    y=data[f'{selected_cluster_2}'],
                    name=f'{selected_cluster_2}',
                    marker_color='lightsalmon'
                ))

                fig.add_trace(go.Bar(
                    x=data['Marker'],
                    y=data[f'{selected_cluster_3}'],
                    name=f'{selected_cluster_3}',
                    marker_color='lightseagreen'
                ))

                fig.update_layout(
                    title=f'{central_tendency} Expression of All Markers Across Selected {cluster_col.capitalize()} Clusters',
                    xaxis_title='Markers',
                    yaxis_title=f'{central_tendency} Expression',
                    barmode='group',
                    xaxis_tickangle=-45,
                    width=1200,
                    height=600,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.subheader(f'Marker Profile Based on Central Tendency')
                st.plotly_chart(fig)

                central_tendency_buffer = io.BytesIO()
                fig.write_image(central_tendency_buffer, format="png")
                central_tendency_buffer.seek(0)
                plot_buffers.append(("Central_Tendency_Bar_Plot", central_tendency_buffer))

                st.download_button(
                    label=f"Download Central Tendency Bar Plot (PNG)",
                    data=central_tendency_buffer,
                    file_name=f"CentralTendency_BarPlot_{random_number}.png",
                    mime="image/png",
                )

                method = 't-test' if test_selected == 'Parametric (T-test)' else 'wilcoxon'
                sc.tl.rank_genes_groups(adata, cluster_col, method=method)

                ranked_markers_1 = adata.uns['rank_genes_groups']['names'][selected_cluster_1]
                ranked_scores_1 = adata.uns['rank_genes_groups']['scores'][selected_cluster_1]

                ranked_markers_2 = adata.uns['rank_genes_groups']['names'][selected_cluster_2]
                ranked_scores_2 = adata.uns['rank_genes_groups']['scores'][selected_cluster_2]

                ranked_markers_3 = adata.uns['rank_genes_groups']['names'][selected_cluster_3]
                ranked_scores_3 = adata.uns['rank_genes_groups']['scores'][selected_cluster_3]

                rank_df = pd.DataFrame({
                    f'Ranked Markers {selected_cluster_1}': ranked_markers_1,
                    f'Ranked Scores {selected_cluster_1}': ranked_scores_1,
                    f'Ranked Markers {selected_cluster_2}': ranked_markers_2,
                    f'Ranked Scores {selected_cluster_2}': ranked_scores_2,
                    f'Ranked Markers {selected_cluster_3}': ranked_markers_3,
                    f'Ranked Scores {selected_cluster_3}': ranked_scores_3
                })

                fig3 = go.Figure()

                fig3.add_trace(go.Bar(
                    x=rank_df[f'Ranked Markers {selected_cluster_1}'],
                    y=rank_df[f'Ranked Scores {selected_cluster_1}'],
                    name=f'{selected_cluster_1}',
                    marker_color='slategray'
                ))

                fig3.add_trace(go.Bar(
                    x=rank_df[f'Ranked Markers {selected_cluster_2}'],
                    y=rank_df[f'Ranked Scores {selected_cluster_2}'],
                    name=f'{selected_cluster_2}',
                    marker_color='lightsalmon'
                ))

                fig3.add_trace(go.Bar(
                    x=rank_df[f'Ranked Markers {selected_cluster_3}'],
                    y=rank_df[f'Ranked Scores {selected_cluster_3}'],
                    name=f'{selected_cluster_3}',
                    marker_color='lightseagreen'
                ))

                fig3.update_layout(
                    title=f'Ranked Markers Comparison Across {cluster_col.capitalize()} Clusters',
                    xaxis_title='Markers',
                    yaxis_title='Scores',
                    barmode='group',
                    xaxis_tickangle=-45,
                    width=1200,
                    height=600,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray',
                    )
                )

                st.subheader(f'Ranked Markers Comparison')
                st.plotly_chart(fig3)

                ranked_markers_buffer = io.BytesIO()
                fig3.write_image(ranked_markers_buffer, format="png")
                ranked_markers_buffer.seek(0)
                plot_buffers.append(("Ranked_Markers_Bar_Plot", ranked_markers_buffer))

                st.download_button(
                    label=f"Download Ranked Markers Bar Plot (PNG)",
                    data=ranked_markers_buffer,
                    file_name=f"RankedMarkers_BarPlot_{random_number}.png",
                    mime="image/png",
                )

                csv_buffer = io.BytesIO()
                rank_df.to_csv(csv_buffer, index=False)
                csv_buffer.seek(0)
                data_buffers.append(("Ranked_Markers_Data", csv_buffer))

                st.download_button(
                    label="Download Ranked Markers Data (CSV)",
                    data=csv_buffer,
                    file_name=f"RankedMarkers_Data_{random_number}.csv",
                    mime="text/csv",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.png",
                            plot_buffer.getvalue(),
                        )
                    for data_name, data_buffer in data_buffers:
                        zip_file.writestr(
                            f"{data_name}_{random_number}.csv",
                            data_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Results as ZIP",
                    data=zip_buffer,
                    file_name=f"All_Results_{random_number}.zip",
                    mime="application/zip",
                )


        with tab20:

            with st.form(key='Plot_20'):
                st.write("Identify Clusters with Marker Combination")

                def calexpression(data, method):
                    if method == 'Mean':
                        return np.mean(data)
                    elif method == 'Median':
                        return np.median(data)

                cluster_marker_dict = {}

                unique_clusters_1 = adata.obs['leiden'].unique()

                for xyz in unique_clusters_1:
                    cluster_mask = adata.obs['leiden'] == xyz
                    cluster_data = adata[cluster_mask]

                    marker_expressions = []
                    markers = adata.var_names
                    for marker in markers:
                        cluster_expression = cluster_data[:, marker].X.toarray().flatten()
                        marker_expression_value = calexpression(cluster_expression, expression_method)
                        marker_expressions.append(marker_expression_value)

                    sorted_idx = np.argsort(marker_expressions)[::-1]
                    sorted_markers = []
                    for i in sorted_idx:
                        marker_name = markers[i]
                        if marker_expressions[i] > 1000:
                            marker_name += "+"
                        elif marker_expressions[i] < 0:
                            marker_name += "-"
                        sorted_markers.append(marker_name)

                    cluster_marker_dict[xyz] = sorted_markers

                sorted_cluster_keys = sorted(cluster_marker_dict.keys(), key=int)

                max_markers = max(len(markers) for markers in cluster_marker_dict.values())

                formatted_df = pd.DataFrame({cluster: cluster_marker_dict[cluster] + [''] * (max_markers - len(cluster_marker_dict[cluster]))
                                            for cluster in sorted_cluster_keys}).T

                formatted_df = formatted_df.reset_index()
                formatted_df.columns = ['Cluster'] + [f'Marker {i+1}' for i in range(max_markers)]
                formatted_df['Cluster'] = 'Cluster ' + formatted_df['Cluster'].astype(str)

                st.dataframe(formatted_df)

                unique_markers = pd.unique(formatted_df.iloc[:, 1:].values.ravel())
                unique_markers = [marker for marker in unique_markers if marker != '']
                selected_markers = st.multiselect("Select markers to filter clusters:", unique_markers)

                submit_button = st.form_submit_button(label='Apply')

            if submit_button:
                if selected_markers:
                    selected_clusters = formatted_df.apply(
                        lambda row: all(marker in row.values for marker in selected_markers), axis=1
                    )

                    clusters_with_markers = formatted_df['Cluster'][selected_clusters].tolist()

                    st.write(f"Clusters where selected markers {selected_markers} are present:")
                    st.write(clusters_with_markers)


        with tab21:

            with st.form(key='Plot_21'):
                st.write("One Cluster vs All Markers Expression:")

                expression_method = st.radio('Select expression calculation method:', ('Mean', 'Median'))

                num_cols = st.slider("Select number of columns for subplots", min_value=1, max_value=10, value=4)

                plot_width = st.slider("Select plot width (in pixels)", min_value=600, max_value=1500, value=1000)
                plot_height = st.slider("Select plot height (in pixels)", min_value=400, max_value=1200, value=600)

                bar_color = st.color_picker("Pick a color for the bars", "#1f77b4")

                submit_button = st.form_submit_button(label='Run Analysis')

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                selected_cluster_mask = adata.obs['leiden'] == cluster_selected
                selected_cluster_data = adata[selected_cluster_mask]
                other_clusters_data = adata[~selected_cluster_mask]

                markers = adata.var_names
                fold_changes = []
                p_values = []
                test_statistics = []
                expr_selected = []
                expr_others = []

                def calculate_expression(data, method):
                    if method == 'Mean':
                        return np.mean(data)
                    elif method == 'Median':
                        return np.median(data)

                for marker in markers:
                    cluster_expression = selected_cluster_data[:, marker].X.toarray().flatten()
                    other_expression = other_clusters_data[:, marker].X.toarray().flatten()

                    fold_change = (
                        np.log2(calculate_expression(cluster_expression, expression_method) + 1) -
                        np.log2(calculate_expression(other_expression, expression_method) + 1)
                    )
                    fold_changes.append(fold_change)

                    if expression_method == 'Parametric (T-test)':
                        test_stat, p_value = ttest_ind(cluster_expression, other_expression, equal_var=False)
                    else:  # Mann-Whitney U test
                        test_stat, p_value = mannwhitneyu(cluster_expression, other_expression, alternative='two-sided')

                    p_values.append(p_value)
                    test_statistics.append(test_stat)
                    expr_selected.append(calculate_expression(cluster_expression, expression_method))
                    expr_others.append(calculate_expression(other_expression, expression_method))

                sorted_idx = np.argsort(expr_selected)[::-1]
                markers_sorted = [markers[i] for i in sorted_idx]
                expr_selected_sorted = [expr_selected[i] for i in sorted_idx]
                expr_others_sorted = [expr_others[i] for i in sorted_idx]

                unique_clusters = adata.obs['leiden'].unique()
                num_clusters = len(unique_clusters)
                num_rows = (num_clusters // num_cols) + (num_clusters % num_cols > 0)

                fig = make_subplots(rows=num_rows, cols=num_cols, subplot_titles=[f'Cluster {c}' for c in unique_clusters])

                for idx, cluster in enumerate(unique_clusters):
                    row = (idx // num_cols) + 1
                    col = (idx % num_cols) + 1

                    cluster_mask = adata.obs['leiden'] == cluster
                    cluster_data = adata[cluster_mask]

                    cluster_expr = [
                        calculate_expression(cluster_data[:, marker].X.toarray().flatten(), expression_method)
                        for marker in markers_sorted
                    ]

                    fig.add_trace(go.Bar(
                        x=markers_sorted,
                        y=cluster_expr,
                        name=f'Cluster {cluster}',
                        marker_color=bar_color
                    ), row=row, col=col)

                fig.update_layout(
                    height=plot_height,
                    width=plot_width,
                    title_text=f'Marker Expression Across Clusters',
                    showlegend=False,
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )

                st.subheader(f'Marker profile based on {expression_method} expression across all clusters')

                st.plotly_chart(fig)

                cluster_expression_buffer = io.BytesIO()
                fig.write_image(cluster_expression_buffer, format="png")
                cluster_expression_buffer.seek(0)
                plot_buffers.append(("Cluster_Expression_Plot", cluster_expression_buffer))

                st.download_button(
                    label=f"Download Marker Expression Plot (PNG)",
                    data=cluster_expression_buffer,
                    file_name=f"MarkerExpression_Plot_{random_number}.png",
                    mime="image/png",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.png",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"All_MarkerExpression_Plots_{random_number}.zip",
                    mime="application/zip",
                )

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'
        with tab22:

            with st.form(key='Plot_22'):
                group1 = st.selectbox('Select Group 1:', adata.obs['Group'].unique())
                group2 = st.selectbox('Select Group 2:', adata.obs['Group'].unique())

                cmap_option_30 = ['Set2', 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                                    'Set1', 'Set3', 'GnBu', 'tab10', 'tab20', 'tab20b', 'tab20c', 'RdBu_r']
                plot_cmap = st.selectbox("Select colormap", options=cmap_option_30,
                                            index=cmap_option_30.index('Set2'), key="plot_cmap")
                plot_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG', 'JPEG'), key='plot_file_format')

                submit_button = st.form_submit_button(label='Submit')

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                filtered_data = adata.obs[adata.obs['Group'].isin([group1, group2])]

                cell_counts = filtered_data.groupby(['SampleID', 'Group', cluster_col]).size().reset_index(name='cell_count')
                hue_order = [group1, group2]

                fig, ax = plt.subplots(figsize=(10, 6))
                sns.boxplot(
                    x=cluster_col, y='cell_count', hue='Group', data=cell_counts, palette=plot_cmap,
                    hue_order=hue_order, ax=ax, dodge=True
                )
                sns.stripplot(
                    x=cluster_col, y='cell_count', hue='Group', data=cell_counts, palette=plot_cmap,
                    hue_order=hue_order, ax=ax, dodge=True, jitter=True, color='black', marker='o', linewidth=0.5, legend=False
                )

                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()

                st.pyplot(fig)

                plot_buffer = io.BytesIO()
                plt.savefig(plot_buffer, format=plot_file_format.lower(), dpi=300, bbox_inches='tight')
                plot_buffer.seek(0)
                plt.close(fig)

                plot_buffers.append(("CellCount_By_Group", plot_buffer))

                st.download_button(
                    label=f"Download Cell Count Plot ({plot_file_format.upper()})",
                    data=plot_buffer,
                    file_name=f"CellCount_celltype_byGroup_{random_number}.{plot_file_format.lower()}",
                    mime=f"image/{plot_file_format.lower()}",
                )


        with tab23:

            with st.form(key='Plot_23'):
                st.subheader("Sankey Diagram for Leiden Clusters and Groups")

                def hex_to_rgba(hex_color, alpha=0.6):
                    rgba = to_rgba(hex_color)
                    return f'rgba({int(rgba[0]*255)}, {int(rgba[1]*255)}, {int(rgba[2]*255)}, {alpha})'

                def create_colored_alluvial_diagram_with_highlight(dataframe, leiden_col, group_col, control_group, cmap, highlight_color, control_group_color):

                    data_counts = dataframe.groupby([leiden_col, group_col]).size().reset_index(name='Count')

                    all_labels = list(data_counts[leiden_col].unique()) + list(data_counts[group_col].unique())
                    leiden_labels = list(data_counts[leiden_col].unique())
                    group_labels = list(data_counts[group_col].unique())

                    source = [all_labels.index(l) for l in data_counts[leiden_col]]
                    target = [all_labels.index(g) for g in data_counts[group_col]]
                    value = data_counts['Count']

                    treatment_group = [g for g in group_labels if g != control_group][0]

                    link_colors = [control_group_color if all_labels[t] == control_group else highlight_color for t in target]

                    num_leiden_clusters = len(leiden_labels)
                    num_groups = len(group_labels)

                    palette_leiden = [hex_to_rgba(to_hex(c)) for c in cmap[:num_leiden_clusters]]
                    palette_groups = [hex_to_rgba(to_hex(c)) for c in cmap[num_leiden_clusters:num_leiden_clusters + num_groups]]

                    node_colors = palette_leiden + palette_groups

                    fig = go.Figure(go.Sankey(
                        node=dict(
                            pad=15,
                            thickness=20,
                            line=dict(color="black", width=1.5),
                            label=all_labels,
                            color=node_colors
                        ),
                        link=dict(
                            source=source,
                            target=target,
                            value=value,
                            color=link_colors
                        )
                    ))

                    fig.update_layout(
                        title_text=f"Alluvial Diagram: Highlighting flows to {treatment_group}",
                        font=dict(size=16),
                        margin=dict(l=50, r=50, t=50, b=50),
                        hoverlabel=dict(bgcolor='white', font=dict(size=14, color='black'))
                    )

                    return fig

                colors_tab20b = list(plt.get_cmap('tab20b').colors)
                colors_tab20c = list(plt.get_cmap('tab20c').colors)
                colors_combined = colors_tab20b + colors_tab20c
                new_cmap = ListedColormap(colors_combined)

                if 'leiden' in adata.obs.columns and 'Group' in adata.obs.columns:

                    available_groups = adata.obs['Group'].unique()
                    control_group = st.selectbox("Select Control Group", options=available_groups)

                    color_options = {
                        "Green": 'rgba(0, 128, 0, 0.5)',
                        "Blue": 'rgba(0, 0, 255, 0.5)',
                        "Red": 'rgba(255, 0, 0, 0.5)',
                        "Purple": 'rgba(128, 0, 128, 0.5)'
                    }
                    selected_color = st.selectbox("Select Highlight Color for Treatment Group", options=list(color_options.keys()), index=0)
                    highlight_color = color_options[selected_color]

                    control_group_color = st.color_picker("Pick a color for the Control Group", '#808080')
                    control_group_color_rgba = hex_to_rgba(control_group_color, alpha=0.2)

                    submit_sankey = st.form_submit_button(label='Generate Sankey Diagram')

            if submit_sankey:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                alluvial_fig = create_colored_alluvial_diagram_with_highlight(
                    adata.obs, 'leiden', 'Group', control_group, colors_combined, highlight_color, control_group_color_rgba
                )

                st.plotly_chart(alluvial_fig)

                sankey_buffer = io.BytesIO()
                alluvial_fig.write_image(sankey_buffer, format="png")
                sankey_buffer.seek(0)
                plot_buffers.append(("Sankey_Plot", sankey_buffer))

                st.download_button(
                    label=f"Download Sankey Plot (PNG)",
                    data=sankey_buffer,
                    file_name=f"SankeyPlot_{random_number}.png",
                    mime="image/png",
                )


        with tab24:

            with st.form("Plot_24"):
                selected_markers = st.multiselect("Select up to 10 markers", adata.var_names, max_selections=10)
                threshold = st.selectbox("Select expression threshold", [1, 100, 1000, 10000])
                cell_type_column = "cell_type" if "cell_type" in adata.obs.columns else "leiden"

                cell_types = adata.obs[cell_type_column].unique().tolist()
                selected_cell_types = st.multiselect("Select cell types/clusters to include in the diagram", cell_types)

                plot_width = st.slider("Select plot width", min_value=500, max_value=1500, value=1000)
                plot_height = st.slider("Select plot height", min_value=300, max_value=1200, value=600)

                submitted = st.form_submit_button("Submit")

            if submitted:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                filtered_data = {}
                for marker in selected_markers:
                    marker_data = adata[:, marker].X > threshold
                    adata_filtered = adata[marker_data, :]
                    filtered_data[marker] = adata_filtered.obs[cell_type_column].value_counts()

                labels = list(selected_markers) + selected_cell_types
                sources = []
                targets = []
                values = []

                cmap = plt.get_cmap("Set2")
                cell_type_colors = {cell_type: cmap(i / len(selected_cell_types)) for i, cell_type in enumerate(selected_cell_types)}
                link_colors = []

                for i, marker in enumerate(selected_markers):
                    for j, cell_type in enumerate(selected_cell_types):
                        value = filtered_data[marker].get(cell_type, 0)
                        if value > 0:
                            sources.append(i)
                            targets.append(len(selected_markers) + j)
                            values.append(value)

                            rgba_color = cell_type_colors[cell_type]
                            hex_color = f"rgba({int(rgba_color[0]*255)}, {int(rgba_color[1]*255)}, {int(rgba_color[2]*255)}, 0.8)"
                            link_colors.append(hex_color)

                fig = go.Figure(go.Sankey(
                    node=dict(
                        pad=15,
                        thickness=20,
                        line=dict(color="black", width=0.5),
                        label=labels,
                    ),
                    link=dict(
                        source=sources,
                        target=targets,
                        value=values,
                        color=link_colors
                    )
                ))

                fig.update_layout(
                    width=plot_width,
                    height=plot_height,
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )

                st.plotly_chart(fig)

                sankey_buffer = io.BytesIO()
                fig.write_image(sankey_buffer, format="pdf")
                sankey_buffer.seek(0)
                plot_buffers.append(("Sankey_Comparison", sankey_buffer))

                st.download_button(
                    label=f"Download Sankey Plot (PDF)",
                    data=sankey_buffer,
                    file_name=f"SankeyComparison_{random_number}.pdf",
                    mime="application/pdf",
                )


        with tab25:

            marker_options = list(adata.var_names)

            with st.form(key="Plot_25"):
                st.write("Select two markers for KDE plot:")

                marker1 = st.selectbox("Marker 1", options=marker_options, key="marker1")
                marker2 = st.selectbox("Marker 2", options=marker_options, key="marker2")

                hue_option = st.selectbox("Hue by", options=["None", "leiden", "cell_type", "Group", "SampleID"], index=0)

                x_min = float(adata[:, marker1].X.min())
                y_min = float(adata[:, marker2].X.min())
                x_max = st.slider("X-axis maximum", min_value=x_min, max_value=float(adata[:, marker1].X.max()), value=float(adata[:, marker1].X.max()))
                y_max = st.slider("Y-axis maximum", min_value=y_min, max_value=float(adata[:, marker2].X.max()), value=float(adata[:, marker2].X.max()))

                submit_button = st.form_submit_button(label="Generate KDE Plot")

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                try:
                    x = adata[:, marker1].X.flatten()
                    y = adata[:, marker2].X.flatten()
                except KeyError:
                    st.error(f"One or both of the markers '{marker1}' or '{marker2}' are not present in the dataset.")
                    st.stop()

                plot_data = pd.DataFrame({marker1: x, marker2: y})
                if hue_option != "None" and hue_option in adata.obs.columns:
                    plot_data[hue_option] = adata.obs[hue_option].values
                else:
                    hue_option = None  # Handle invalid or missing hue_option

                fig1, ax1 = plt.subplots(figsize=(8, 6))
                sns.kdeplot(x=plot_data[marker1], y=plot_data[marker2], ax=ax1)
                ax1.set_xlabel(marker1)
                ax1.set_ylabel(marker2)
                ax1.set_xlim(x_min, x_max)
                ax1.set_ylim(y_min, y_max)
                ax1.set_title(f"2D KDE plot of {marker1} vs {marker2}")
                st.pyplot(fig1)

                kde_plot1_buffer = io.BytesIO()
                fig1.savefig(kde_plot1_buffer, format="png", dpi=300, bbox_inches="tight")
                kde_plot1_buffer.seek(0)
                plt.close(fig1)
                plot_buffers.append(("2D_KDE_Plot", kde_plot1_buffer))

                st.download_button(
                    label=f"Download 2D KDE Plot ({marker1} vs {marker2})",
                    data=kde_plot1_buffer,
                    file_name=f"2D_KDE_Plot_{marker1}_vs_{marker2}_{random_number}.png",
                    mime="image/png",
                )

                if hue_option:
                    fig2, ax2 = plt.subplots(figsize=(8, 6))
                    sns.kdeplot(
                        data=plot_data,
                        x=marker1,
                        y=marker2,
                        hue=hue_option,
                        fill=True,
                        common_norm=False,
                        palette="tab20c",
                        alpha=0.5,
                        linewidth=0,
                        ax=ax2
                    )
                    ax2.set_xlabel(marker1)
                    ax2.set_ylabel(marker2)
                    ax2.set_xlim(x_min, x_max)
                    ax2.set_ylim(y_min, y_max)
                    ax2.set_title(f"Enhanced 2D KDE plot of {marker1} vs {marker2} by {hue_option}")
                    st.pyplot(fig2)

                    kde_plot2_buffer = io.BytesIO()
                    fig2.savefig(kde_plot2_buffer, format="png", dpi=300, bbox_inches="tight")
                    kde_plot2_buffer.seek(0)
                    plt.close(fig2)
                    plot_buffers.append(("Enhanced_2D_KDE_Plot", kde_plot2_buffer))

                    st.download_button(
                        label=f"Download Enhanced KDE Plot ({marker1} vs {marker2} by {hue_option})",
                        data=kde_plot2_buffer,
                        file_name=f"Enhanced_KDE_Plot_{marker1}_vs_{marker2}_by_{hue_option}_{random_number}.png",
                        mime="image/png",
                    )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.png",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All KDE Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"KDE_Plots_{random_number}.zip",
                    mime="application/zip",
                )


        with tab26:

            heatmap_cmap_options = [
                'Blues', 'plasma', 'inferno', 'magma', 'cividis', 'viridis',
                'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
                'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
                'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'coolwarm', 'Spectral'
            ]

            with st.form(key='Plot_26'):
                st.write("Plot Total Number of Cells vs Clusters:")

                chart_width = st.slider("Select chart width", min_value=1.0, max_value=10.0, value=4.5, step=0.1, key="cell_vs_cluster_chart_width")
                chart_height = st.slider("Select chart height", min_value=1.0, max_value=10.0, value=3.0, step=0.1, key="cell_vs_cluster_chart_height")
                font_size = st.slider("Select the font size for labels", min_value=1, max_value=20, value=5, step=1, key="cell_vs_cluster_font_size")

                plot_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG'), key='cell_vs_cluster_plot_file_format')

                submit_plot = st.form_submit_button("Plot Cells vs Clusters")

            if submit_plot:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                cell_counts_per_cluster = adata.obs['leiden'].value_counts().reset_index()
                cell_counts_per_cluster.columns = ['Cluster', 'Cell Count']

                fig = px.bar(
                    cell_counts_per_cluster, x='Cluster', y='Cell Count',
                    title='Total Number of Cells in Each Cluster',
                    labels={'Cluster': 'Cluster', 'Cell Count': 'Number of Cells'},
                    width=int(chart_width * 100), height=int(chart_height * 100)
                )

                fig.update_layout(
                    xaxis_title_font=dict(size=font_size),
                    yaxis_title_font=dict(size=font_size),
                    font=dict(size=font_size)
                )

                st.plotly_chart(fig)

                plot_buffer = io.BytesIO()
                fig.write_image(plot_buffer, format=plot_file_format.lower(), engine='kaleido')
                plot_buffer.seek(0)
                plot_buffers.append(("ClusterVsCells_BarPlot", plot_buffer))

                st.download_button(
                    label=f"Download Cluster vs Cells Plot ({plot_file_format.upper()})",
                    data=plot_buffer,
                    file_name=f"ClusterVsCellsPlot_{random_number}.{plot_file_format.lower()}",
                    mime=f"image/{plot_file_format.lower()}",
                )


        with tab27:

            with st.form(key='Plot_27'):
                st.write("Compare Relative Abundance of Cells in Selected Clusters")

                marker_x = st.selectbox("Select Marker 1", options=adata.var_names)
                marker_y = st.selectbox("Select Marker 2", options=adata.var_names)

                selected_cluster_1 = st.selectbox("Select Leiden Cluster 1", options=adata.obs['leiden'].unique())
                selected_cluster_2 = st.selectbox("Select Leiden Cluster 2", options=adata.obs['leiden'].unique())

                submit_button = st.form_submit_button(label='Submit')

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                # Calculate overall means
                overall_mean_x = adata[:, marker_x].X.mean()
                overall_mean_y = adata[:, marker_y].X.mean()

                # Get data for the selected clusters
                cluster_data_1 = adata[adata.obs['leiden'] == selected_cluster_1, :]
                cluster_data_2 = adata[adata.obs['leiden'] == selected_cluster_2, :]

                total_cells_cluster_1 = cluster_data_1.n_obs
                total_cells_cluster_2 = cluster_data_2.n_obs

                marker_x_expression_1 = cluster_data_1[:, marker_x].X.flatten()
                marker_y_expression_1 = cluster_data_1[:, marker_y].X.flatten()

                marker_x_expression_2 = cluster_data_2[:, marker_x].X.flatten()
                marker_y_expression_2 = cluster_data_2[:, marker_y].X.flatten()

                # Calculate relative abundances
                num_cells_above_mean_x_1 = (marker_x_expression_1 > overall_mean_x).sum()
                num_cells_above_mean_y_1 = (marker_y_expression_1 > overall_mean_y).sum()

                num_cells_above_mean_x_2 = (marker_x_expression_2 > overall_mean_x).sum()
                num_cells_above_mean_y_2 = (marker_y_expression_2 > overall_mean_y).sum()

                rel_abundance_x_1 = num_cells_above_mean_x_1 / total_cells_cluster_1
                rel_abundance_y_1 = num_cells_above_mean_y_1 / total_cells_cluster_1

                rel_abundance_x_2 = num_cells_above_mean_x_2 / total_cells_cluster_2
                rel_abundance_y_2 = num_cells_above_mean_y_2 / total_cells_cluster_2

                bar_data = pd.DataFrame({
                    'Marker': [marker_x, marker_y] * 2,
                    'Cluster': [selected_cluster_1] * 2 + [selected_cluster_2] * 2,
                    'RelativeAbundance': [rel_abundance_x_1, rel_abundance_y_1, rel_abundance_x_2, rel_abundance_y_2]
                })

                fig = go.Figure()

                fig.add_trace(go.Bar(
                    x=bar_data[bar_data['Cluster'] == selected_cluster_1]['Marker'],
                    y=bar_data[bar_data['Cluster'] == selected_cluster_1]['RelativeAbundance'],
                    name=f'Cluster {selected_cluster_1}',
                    marker_color='slategray'
                ))

                fig.add_trace(go.Bar(
                    x=bar_data[bar_data['Cluster'] == selected_cluster_2]['Marker'],
                    y=bar_data[bar_data['Cluster'] == selected_cluster_2]['RelativeAbundance'],
                    name=f'Cluster {selected_cluster_2}',
                    marker_color='tomato'
                ))

                fig.update_layout(
                    title=f'Relative Abundance of Cells with Expression Above Overall Mean in Clusters {selected_cluster_1} and {selected_cluster_2}',
                    xaxis_title='Marker',
                    yaxis_title='Relative Abundance',
                    barmode='group',
                    legend_title_text='Leiden Cluster',
                    plot_bgcolor='rgba(0,0,0,0)'
                )

                st.plotly_chart(fig)

                bar_plot_buffer = io.BytesIO()
                fig.write_image(bar_plot_buffer, format="png", engine="kaleido")
                bar_plot_buffer.seek(0)
                plot_buffers.append(("RelativeAbundance_Comparison", bar_plot_buffer))

                # Second bar plot
                overall_means = adata.X.mean(axis=0)

                cluster_1_data = adata[adata.obs['leiden'] == selected_cluster_1, :].X
                cluster_2_data = adata[adata.obs['leiden'] == selected_cluster_2, :].X

                num_cells_above_mean_1 = (cluster_1_data > overall_means).sum(axis=0)
                num_cells_above_mean_2 = (cluster_2_data > overall_means).sum(axis=0)

                rel_abundance_1 = num_cells_above_mean_1 / total_cells_cluster_1
                rel_abundance_2 = num_cells_above_mean_2 / total_cells_cluster_2

                fig_comparison = go.Figure()

                fig_comparison.add_trace(go.Bar(
                    x=adata.var_names,
                    y=rel_abundance_1,
                    name=f'Cluster {selected_cluster_1}',
                    marker_color='slategray'
                ))

                fig_comparison.add_trace(go.Bar(
                    x=adata.var_names,
                    y=rel_abundance_2,
                    name=f'Cluster {selected_cluster_2}',
                    marker_color='tomato'
                ))

                fig_comparison.update_layout(
                    title=f'Relative Abundance for All Markers Between Cluster {selected_cluster_1} and Cluster {selected_cluster_2}',
                    xaxis_title='Marker',
                    yaxis_title='Relative Abundance',
                    barmode='group',
                    plot_bgcolor='rgba(0,0,0,0)'
                )

                st.plotly_chart(fig_comparison)

                comparison_plot_buffer = io.BytesIO()
                fig_comparison.write_image(comparison_plot_buffer, format="png", engine="kaleido")
                comparison_plot_buffer.seek(0)
                plot_buffers.append(("AllMarkers_RelativeAbundance", comparison_plot_buffer))

                st.download_button(
                    label=f"Download Relative Abundance Plot (Cluster Comparison)",
                    data=bar_plot_buffer,
                    file_name=f"RelativeAbundance_Comparison_{random_number}.png",
                    mime="image/png",
                )

                st.download_button(
                    label=f"Download Relative Abundance Plot (All Markers)",
                    data=comparison_plot_buffer,
                    file_name=f"AllMarkers_RelativeAbundance_{random_number}.png",
                    mime="image/png",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.png",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"RelativeAbundance_Plots_{random_number}.zip",
                    mime="application/zip",
                )

                st.write(f"Total number of cells in Cluster {selected_cluster_1}: {total_cells_cluster_1}")
                st.write(f"Total number of cells in Cluster {selected_cluster_2}: {total_cells_cluster_2}")


        with tab28:

            with st.form(key='Plot_28'):
                st.write("Correlation Matrix Plot:")

                chart_width = st.slider("Select chart width", min_value=1.0, max_value=30.0, value=10.0, step=0.1, key="correlation_matrix_chart_width")
                chart_height = st.slider("Select chart height", min_value=1.0, max_value=25.0, value=8.0, step=0.1, key="correlation_matrix_chart_height")
                markers = st.multiselect("Select markers for correlation", options=adata.var_names, key="markers_for_correlation")
                correlation_matrix_plot_file_format = st.radio("Select file format to save the figure", ('PNG', 'PDF', 'SVG'), key='correlation_matrix_plot_file_format')
                show_numbers = st.checkbox("Show correlation values on the heatmap", value=True, key="show_numbers")
                submit_correlation_matrix_plot = st.form_submit_button("Plot Correlation Matrix")

            if submit_correlation_matrix_plot:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                if len(markers) < 2:
                    st.error("Please select at least two markers.")
                else:
                    marker_data = pd.DataFrame(adata[:, markers].X, columns=markers)
                    correlation_matrix = marker_data.corr()

                    fig, ax = plt.subplots(figsize=(chart_width, chart_height))
                    sns.heatmap(correlation_matrix, annot=show_numbers, cmap="coolwarm", ax=ax)
                    st.pyplot(fig)

                    correlation_matrix_buffer = io.BytesIO()
                    fig.savefig(correlation_matrix_buffer, format=correlation_matrix_plot_file_format.lower(), dpi=300, bbox_inches="tight")
                    correlation_matrix_buffer.seek(0)
                    plt.close(fig)
                    plot_buffers.append(("CorrelationMatrix", correlation_matrix_buffer))

                    st.download_button(
                        label=f"Download Correlation Matrix Plot ({correlation_matrix_plot_file_format.upper()})",
                        data=correlation_matrix_buffer,
                        file_name=f"CorrelationMatrix_{random_number}.{correlation_matrix_plot_file_format.lower()}",
                        mime=f"image/{correlation_matrix_plot_file_format.lower()}",
                    )


        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab29:


            with st.form(key='Plot_29'):
                st.write("Plot Total Number of Cells in Each Cluster by Group:")

                cluster = st.selectbox(f"Select {cluster_col.capitalize()} Cluster", options=adata.obs[cluster_col].unique(), key="cluster_selector_ab")

                chart_width = st.slider("Select chart width", min_value=1.0, max_value=10.0, value=3.0, step=0.1, key="cellcount_chart_width_ab")
                chart_height = st.slider("Select chart height", min_value=1.0, max_value=10.0, value=5.0, step=0.1, key="cellcount_chart_height_ab")
                font_size = st.slider("Select the font size for labels", min_value=1, max_value=20, value=11, step=1, key="cellcount_font_size_ab")

                cellcount_plot_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG'), key='cellcount_plot_file_format_ab')

                submit_cell_count_plot = st.form_submit_button("Plot Cell Count")

            if submit_cell_count_plot:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                selected_data = adata.obs[adata.obs[cluster_col] == cluster]
                cell_counts = selected_data.groupby('Group').size().reset_index(name='Cell Count')

                fig_bar = px.bar(
                    cell_counts,
                    x='Group',
                    y='Cell Count',
                    title=f'Total Cells in {cluster_col.capitalize()} Cluster {cluster}',
                    labels={'Group': 'Group', 'Cell Count': 'Number of Cells'},
                    width=int(chart_width * 100),
                    height=int(chart_height * 100)
                )

                fig_bar.update_layout(
                    xaxis_title_font=dict(size=font_size),
                    yaxis_title_font=dict(size=font_size),
                    font=dict(size=font_size)
                )

                st.plotly_chart(fig_bar)

                bar_plot_buffer = io.BytesIO()
                fig_bar.write_image(bar_plot_buffer, format="png", engine="kaleido")
                bar_plot_buffer.seek(0)
                plot_buffers.append(("CellCount_BarPlot", bar_plot_buffer))

                st.download_button(
                    label=f"Download Bar Plot ({cluster_col.capitalize()} Cluster {cluster})",
                    data=bar_plot_buffer,
                    file_name=f"CellCount_BarPlot_{cluster_col}_{cluster}_{random_number}.png",
                    mime="image/png",
                )

                cell_counts['Frequency'] = cell_counts['Cell Count'] / cell_counts['Cell Count'].sum()

                fig_pie = px.pie(
                    cell_counts,
                    values='Cell Count',
                    names='Group',
                    title=f'Frequency of Cells in {cluster_col.capitalize()} Cluster {cluster} (as %)',
                    width=int(chart_width * 100),
                    height=int(chart_height * 100),
                    hole=0.4
                )

                fig_pie.update_layout(
                    font=dict(size=font_size),
                    annotations=[dict(text=f'Total: {cell_counts["Cell Count"].sum()}', x=0.5, y=0.5, font_size=font_size, showarrow=False)]
                )

                st.plotly_chart(fig_pie)

                pie_plot_buffer = io.BytesIO()
                fig_pie.write_image(pie_plot_buffer, format="png", engine="kaleido")
                pie_plot_buffer.seek(0)
                plot_buffers.append(("CellFrequency_PieChart", pie_plot_buffer))

                st.download_button(
                    label=f"Download Pie Chart ({cluster_col.capitalize()} Cluster {cluster})",
                    data=pie_plot_buffer,
                    file_name=f"CellFrequency_PieChart_{cluster_col}_{cluster}_{random_number}.png",
                    mime="image/png",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.png",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"CellCountPlots_{cluster_col}_{cluster}_{random_number}.zip",
                    mime="application/zip",
                )


        with tab30:

            with st.form(key='Plot_30'):
                st.write("Violin Plot and Bar Plot:")

                chart_width = st.slider("Select chart width", min_value=1.0, max_value=10.0, value=4.5, step=0.1, key="violin_chart_width")
                chart_height = st.slider("Select chart height", min_value=1.0, max_value=10.0, value=3.0, step=0.1, key="violin_chart_height")
                font_size = st.slider("Select the font size for labels", min_value=6, max_value=20, value=6, step=1, key="violin_font_size")
                marker = st.selectbox("Select a marker to plot", options=adata.var_names, key="violin_marker_selector")
                group_by = st.selectbox("Group by", options=adata.obs.columns, key="violin_groupby_selector")
                violin_plot_file_format = st.radio("Select file format to save the figures", ('PNG', 'PDF', 'SVG'), key='violin_plot_file_format')
                submit_violin_and_bar_plot = st.form_submit_button("Plot Violin and Bar Plots")

            if submit_violin_and_bar_plot:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                unique_categories = adata.obs[group_by].unique()
                palette = sns.color_palette("husl", n_colors=len(unique_categories))

                with plt.rc_context({'figure.figsize': (chart_width, chart_height)}):
                    ax = sc.pl.violin(adata, marker, groupby=group_by, show=False, palette=palette)

                ax.xaxis.label.set_size(font_size)
                ax.yaxis.label.set_size(font_size)
                ax.tick_params(labelsize=font_size)

                for tick in ax.get_xticklabels():
                    tick.set_rotation(90)

                fig = ax.figure
                st.pyplot(fig)

                violin_plot_buffer = io.BytesIO()
                fig.savefig(violin_plot_buffer, format=violin_plot_file_format.lower(), dpi=300, bbox_inches="tight")
                violin_plot_buffer.seek(0)
                plt.close(fig)
                plot_buffers.append(("ViolinPlot", violin_plot_buffer))

                st.download_button(
                    label=f"Download Violin Plot ({marker})",
                    data=violin_plot_buffer,
                    file_name=f"ViolinPlot_{marker}_{random_number}.{violin_plot_file_format.lower()}",
                    mime=f"image/{violin_plot_file_format.lower()}",
                )

                bar_data = pd.DataFrame({group_by: adata.obs[group_by], 'expression': adata[:, marker].X.toarray().flatten()})
                bar_data_mean = bar_data.groupby(group_by).mean().reset_index()

                bar_fig = px.bar(
                    bar_data_mean,
                    x=group_by,
                    y='expression',
                    title=f'Bar Plot of {marker}',
                    labels={group_by: group_by, 'expression': 'Mean Expression'}
                )

                st.plotly_chart(bar_fig)

                bar_plot_buffer = io.BytesIO()
                bar_fig.write_image(bar_plot_buffer, format="png", engine="kaleido")
                bar_plot_buffer.seek(0)
                plot_buffers.append(("BarPlot", bar_plot_buffer))

                st.download_button(
                    label=f"Download Bar Plot ({marker})",
                    data=bar_plot_buffer,
                    file_name=f"BarPlot_{marker}_{random_number}.png",
                    mime="image/png",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.{violin_plot_file_format.lower() if plot_name == 'ViolinPlot' else 'png'}",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"ViolinAndBarPlots_{marker}_{random_number}.zip",
                    mime="application/zip",
                )

        st.divider()

        st.subheader("Statistics")

        tab31, tab32, tab33, tab34, tab35, tab36, tab37 = st.tabs(["Plot31", "Plot32", "Plot33", "Plot34", "Plot35", "Plot36", "Plot37"])

        with tab31:

            with st.form(key="Plot_31"):
                if 'cell_type' in adata.obs.columns:
                    comparison_options = adata.obs['cell_type'].unique()
                    comparison_label = "Select Cell Types:"
                else:
                    comparison_options = adata.obs['leiden'].unique()
                    comparison_label = "Select Leiden Clusters:"

                selected_comparisons = st.multiselect(
                    comparison_label,
                    options=comparison_options
                )

                group_options = adata.obs['Group'].unique()
                group1 = st.selectbox("Select Group 1:", options=group_options)
                group2 = st.selectbox("Select Group 2:", options=group_options)

                group1_color = st.color_picker(f"Pick color for {group1}", value='#1f77b4')  # Default blue color
                group2_color = st.color_picker(f"Pick color for {group2}", value='#ff7f0e')  # Default orange color

                box_option = st.radio(
                    "Select what to show in plot:",
                    options=["Median (IQR)", "Mean (SD)", "Both"]
                )

                test_option = st.radio(
                    "Select statistical test:",
                    options=["T-test", "Welch's T-test", "Mann-Whitney U test", "Brunner-Munzel test"]
                )

                show_p_value = st.checkbox("Show p-value on plots", value=True)
                plot_width = st.slider("Plot Width", min_value=200, max_value=1200, value=600)
                plot_height = st.slider("Plot Height", min_value=200, max_value=1200, value=600)
                num_columns = st.slider("Number of Columns", min_value=1, max_value=5, value=2)
                font_size = st.slider("Font Size", min_value=8, max_value=24, value=14)

                submit_comparison = st.form_submit_button(label="Compare Groups")

            if submit_comparison:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                num_rows = (len(selected_comparisons) // num_columns) + 1
                fig = make_subplots(rows=num_rows, cols=num_columns, subplot_titles=selected_comparisons)

                for idx, comparison in enumerate(selected_comparisons):
                    filtered_data = adata.obs[
                        (adata.obs['Group'].isin([group1, group2])) &
                        ((adata.obs['cell_type'] == comparison) if 'cell_type' in adata.obs.columns else
                        (adata.obs['leiden'] == comparison))
                    ]

                    cell_counts_per_sample = filtered_data.groupby(['Group', 'SampleID']).size().reset_index(name='cell_count')

                    group1_data = cell_counts_per_sample[cell_counts_per_sample['Group'] == group1]['cell_count']
                    group2_data = cell_counts_per_sample[cell_counts_per_sample['Group'] == group2]['cell_count']

                    if group1_data.empty or group2_data.empty:
                        st.error(f"No data available for comparison in {comparison}. Skipping this cluster/cell type.")
                        continue

                    if test_option == "T-test":
                        stat, p_value = ttest_ind(group1_data, group2_data)
                    elif test_option == "Welch's T-test":
                        stat, p_value = ttest_ind(group1_data, group2_data, equal_var=False)
                    elif test_option == "Mann-Whitney U test":
                        stat, p_value = mannwhitneyu(group1_data, group2_data)
                    elif test_option == "Brunner-Munzel test":
                        stat, p_value = brunnermunzel(group1_data, group2_data)

                    row = (idx // num_columns) + 1
                    col = (idx % num_columns) + 1

                    if box_option == "Median (IQR)":
                        fig.add_trace(go.Box(
                            y=group1_data,
                            name=group1,
                            boxpoints='all',
                            jitter=0.5,
                            pointpos=0,
                            marker_color=group1_color,
                            boxmean=False
                        ), row=row, col=col)

                        fig.add_trace(go.Box(
                            y=group2_data,
                            name=group2,
                            boxpoints='all',
                            jitter=0.5,
                            pointpos=0,
                            marker_color=group2_color,
                            boxmean=False
                        ), row=row, col=col)

                    elif box_option == "Mean (SD)":
                        fig.add_trace(go.Scatter(
                            x=[group1] * len(group1_data),
                            y=group1_data,
                            mode='markers',
                            name=f'{group1} data',
                            marker_color=group1_color
                        ), row=row, col=col)

                        fig.add_trace(go.Scatter(
                            x=[group2] * len(group2_data),
                            y=group2_data,
                            mode='markers',
                            name=f'{group2} data',
                            marker_color=group2_color
                        ), row=row, col=col)

                    elif box_option == "Both":
                        fig.add_trace(go.Box(
                            y=group1_data,
                            name=group1,
                            boxpoints='all',
                            jitter=0.5,
                            pointpos=0,
                            marker_color=group1_color,
                            boxmean='sd'
                        ), row=row, col=col)

                        fig.add_trace(go.Box(
                            y=group2_data,
                            name=group2,
                            boxpoints='all',
                            jitter=0.5,
                            pointpos=0,
                            marker_color=group2_color,
                            boxmean='sd'
                        ), row=row, col=col)

                    if show_p_value:
                        fig.add_annotation(
                            text=f"p = {p_value:.4f}",
                            xref="x domain", yref="y domain",
                            x=0.5, y=1.05,
                            showarrow=False, font=dict(size=font_size),
                            row=row, col=col
                        )

                fig.update_layout(
                    title=f"Comparison of {len(selected_comparisons)} Comparisons",
                    width=plot_width,
                    height=plot_height,
                    font=dict(size=font_size),
                    showlegend=False,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    yaxis=dict(
                        gridcolor='lightgray'),
                    title_font_size=font_size
                )

                st.plotly_chart(fig)

                comparison_buffer = io.BytesIO()
                pio.write_image(fig, comparison_buffer, format="pdf")
                comparison_buffer.seek(0)
                plot_buffers.append(("SubplotComparison", comparison_buffer))

                st.download_button(
                    label=f"Download Comparison Plot (PDF)",
                    data=comparison_buffer,
                    file_name=f"SubplotComparison_{random_number}.pdf",
                    mime="application/pdf",
                )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(
                            f"{plot_name}_{random_number}.pdf",
                            plot_buffer.getvalue(),
                        )
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"ComparisonPlots_{random_number}.zip",
                    mime="application/zip",
                )

        cluster_label = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab32:

            with st.form(key="Plot_32"):
                st.write("Example: CD8+ expressing T cells compared between groups for Marker(hi) expression")
                selected_clusters = st.multiselect(
                    f"Select {cluster_label} to include in analysis:",
                    options=adata.obs[cluster_label].unique(),
                    default=adata.obs[cluster_label].unique()[0]
                )

                filtered_adata = adata[adata.obs[cluster_label].isin(selected_clusters)]

                selected_markers = st.multiselect(
                    "Select markers for analysis:",
                    options=filtered_adata.var_names
                )

                threshold = st.selectbox(
                    "Select a threshold value for marker expression:",
                    options=[100, 1000, 10000],
                    index=1
                )

                group_options = filtered_adata.obs['Group'].unique()
                group1 = st.selectbox("Select Group 1 for comparison:", options=group_options)
                group2 = st.selectbox("Select Group 2 for comparison:", options=group_options)

                color_group1 = st.color_picker(f"Select color for {group1}:", value='#1f77b4')
                color_group2 = st.color_picker(f"Select color for {group2}:", value='#ff7f0e')

                test_option = st.radio(
                    "Select statistical test:",
                    options=["T-test", "Welch's T-test", "Mann-Whitney U test"]
                )

                show_p_value = st.checkbox("Show p-value in plot", value=True)

                num_columns = st.slider("Select number of columns for subplots:", min_value=1, max_value=4, value=2)

                fig_width = st.slider("Select figure width:", min_value=400, max_value=2000, value=1200)
                fig_height = st.slider("Select figure height:", min_value=400, max_value=2000, value=900)

                submit_button = st.form_submit_button(label="Run Analysis")

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                num_rows = (len(selected_markers) + num_columns - 1) // num_columns

                fig = make_subplots(
                    rows=num_rows,
                    cols=num_columns,
                    subplot_titles=selected_markers
                )

                for idx, selected_marker in enumerate(selected_markers):
                    df_filtered = pd.DataFrame(
                        filtered_adata.X[:, filtered_adata.var_names.get_loc(selected_marker)],
                        columns=[selected_marker],
                        index=filtered_adata.obs.index
                    )

                    df_filtered = df_filtered[df_filtered[selected_marker] >= threshold]
                    df_filtered['SampleID'] = filtered_adata.obs['SampleID']
                    df_filtered['Group'] = filtered_adata.obs['Group']

                    df_grouped = df_filtered.groupby(['SampleID', 'Group'])[selected_marker].mean().reset_index()
                    df_grouped_cleaned = df_grouped.dropna(subset=[selected_marker])

                    group1_data = df_grouped_cleaned[df_grouped_cleaned['Group'] == group1][selected_marker]
                    group2_data = df_grouped_cleaned[df_grouped_cleaned['Group'] == group2][selected_marker]

                    if test_option == "T-test":
                        stat, p_value = ttest_ind(group1_data, group2_data)
                    elif test_option == "Welch's T-test":
                        stat, p_value = ttest_ind(group1_data, group2_data, equal_var=False)
                    else:
                        stat, p_value = mannwhitneyu(group1_data, group2_data)

                    row = (idx // num_columns) + 1
                    col = (idx % num_columns) + 1

                    fig.add_trace(go.Box(
                        y=group1_data,
                        name=group1,
                        marker_color=color_group1,
                        boxpoints='all',
                        jitter=0.5,
                        pointpos=0
                    ), row=row, col=col)

                    fig.add_trace(go.Box(
                        y=group2_data,
                        name=group2,
                        marker_color=color_group2,
                        boxpoints='all',
                        jitter=0.5,
                        pointpos=0
                    ), row=row, col=col)

                    if show_p_value:
                        fig.add_annotation(
                            text=f"p-value = {p_value:.4f}",
                            xref="x domain", yref="y domain",
                            x=0.5, y=1.0,
                            showarrow=False,
                            font=dict(size=12),
                            row=row, col=col
                        )

                fig.update_layout(
                    height=fig_height,
                    width=fig_width,
                    title_text="Marker Expression Comparisons",
                    showlegend=False,
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                )

                for i in range(1, num_rows * num_columns + 1):
                    fig.update_yaxes(tickformat=',', gridcolor='lightgray', row=(i - 1) // num_columns + 1, col=(i - 1) % num_columns + 1)

                st.plotly_chart(fig)

                subplot_comparison_buffer = io.BytesIO()
                fig.write_image(subplot_comparison_buffer, format="pdf")
                subplot_comparison_buffer.seek(0)
                plot_buffers.append(("SubplotComparison", subplot_comparison_buffer))

                st.download_button(
                    label=f"Download Marker Expression Comparisons (PDF)",
                    data=subplot_comparison_buffer,
                    file_name=f"MarkerExpressionComparison_{random_number}.pdf",
                    mime="application/pdf",
                )


        def cohen_d(group1_data, group2_data):
            mean1, mean2 = np.mean(group1_data), np.mean(group2_data)
            pooled_std = np.sqrt(((len(group1_data) - 1) * np.std(group1_data, ddof=1) ** 2 +
                                  (len(group2_data) - 1) * np.std(group2_data, ddof=1) ** 2) /
                                 (len(group1_data) + len(group2_data) - 2))
            return (mean1 - mean2) / pooled_std

        cluster_label = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'
        with tab33:
            with st.form(key='Plot_33'):
                group1 = st.selectbox('Select Control Group:', adata.obs['Group'].unique())
                group2 = st.selectbox('Select Treatment Group:', adata.obs['Group'].unique())
                effect_size_choice = st.radio('Choose Effect Size Calculation:', ("Cohen's d", "Cliff's Delta"))
                bar_color = st.color_picker('Pick a Bar Color', '#badcbe')
                submit_button = st.form_submit_button(label='Calculate Effect Size')

            if submit_button:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                cell_counts = adata.obs.groupby(['SampleID', 'Group', cluster_label]).size().reset_index(name='cell_count')
                filtered_data = cell_counts[cell_counts['Group'].isin([group1, group2])]

                if filtered_data.empty:
                    st.warning("No data found for the selected groups!")
                else:
                    results = []
                    for cell_type in filtered_data[cluster_label].unique():
                        control_data = filtered_data[(filtered_data['Group'] == group2) &
                                                     (filtered_data[cluster_label] == cell_type)]['cell_count']
                        treatment_data = filtered_data[(filtered_data['Group'] == group1) &
                                                       (filtered_data[cluster_label] == cell_type)]['cell_count']

                        if len(control_data) > 0 and len(treatment_data) > 0:
                            if effect_size_choice == "Cohen's d":
                                d_value = cohen_d(control_data, treatment_data)
                                results.append({
                                    'cell_type': cell_type,
                                    'effect_size_type': "Cohen's d",
                                    'value': d_value
                                })
                            elif effect_size_choice == "Cliff's Delta":
                                delta, _ = cliffs_delta(control_data, treatment_data)
                                results.append({
                                    'cell_type': cell_type,
                                    'effect_size_type': "Cliff's Delta",
                                    'value': delta
                                })

                    if not results:
                        st.warning("No valid effect sizes calculated!")
                    else:
                        results_df = pd.DataFrame(results)
                        fig = px.bar(
                            results_df,
                            x='cell_type',
                            y='value',
                            color='effect_size_type',
                            title=f"Effect Sizes for Control vs Treatment using {effect_size_choice}",
                            labels={'value': 'Effect Size', 'cell_type': 'Cell Type'},
                            color_discrete_sequence=[bar_color]
                        )
                        fig.update_layout(
                            showlegend=False,
                            autosize=True,
                            xaxis_tickangle=-90,
                            plot_bgcolor='white',
                            paper_bgcolor='white',
                            yaxis=dict(gridcolor='lightgray')
                        )
                        st.plotly_chart(fig)

                        effect_size_buffer = io.BytesIO()
                        fig.write_image(effect_size_buffer, format="pdf")
                        effect_size_buffer.seek(0)
                        plot_buffers.append(("EffectSizeComparison", effect_size_buffer))

                        st.download_button(
                            label=f"Download Effect Size Comparison Plot (PDF)",
                            data=effect_size_buffer,
                            file_name=f"EffectSizeComparison_{random_number}.pdf",
                            mime="application/pdf",
                        )

        def create_contingency_table(adata, cluster_col):
            contingency_table = pd.crosstab(adata.obs['Group'], adata.obs[cluster_col])
            return contingency_table

        def chi_square_test_with_residuals(contingency_table):
            chi2, p_value, dof, expected = chi2_contingency(contingency_table)
            residuals = (contingency_table - expected) / np.sqrt(expected)
            return chi2, p_value, dof, expected, residuals

        cluster_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'

        with tab34:
            with st.form("Plot_34"):
                st.write("**Is there a statistically significant association between clusters and groups?**")

                submitted = st.form_submit_button("Run Chi-Square Test")

            if submitted:
                random_number = random.randint(1000, 9999)
                file_buffers = []

                contingency_table = create_contingency_table(adata, cluster_col)
                st.write("Contingency Table (Clusters vs Groups):")
                st.table(contingency_table)

                chi2, p_value, dof, expected, residuals = chi_square_test_with_residuals(contingency_table)
                st.write(f"Chi-Square Test Results:")
                st.write(f"Chi-Square Statistic: {chi2}")
                st.write(f"p-value: {p_value}")
                st.write(f"Degrees of Freedom: {dof}")

                if p_value < 0.05:
                    st.write("There is a statistically significant association between clusters and groups (p < 0.05).")
                    st.write("Expected Counts:")
                    st.table(pd.DataFrame(expected, index=contingency_table.index, columns=contingency_table.columns))

                    st.subheader("Which clusters are more prevalent in their respective groups?")
                    st.write("Standardized Residuals (Positive residuals indicate higher prevalence than expected):")
                    st.table(residuals)
                    st.write("Clusters with positive standardized residuals are more prevalent in their respective groups.")

                    st.write("Heatmap of Standardized Residuals:")
                    plt.figure(figsize=(12, 3))
                    sns.heatmap(residuals, annot=True, cmap="viridis", center=0, cbar_kws={'label': 'Standardized Residuals'})
                    plt.title("Standardized Residuals Heatmap")
                    st.pyplot(plt.gcf())

                    contingency_buffer = io.StringIO()
                    contingency_table.to_csv(contingency_buffer)
                    contingency_buffer.seek(0)
                    file_buffers.append(("ContingencyTable", contingency_buffer))

                    residuals_buffer = io.StringIO()
                    residuals.to_csv(residuals_buffer)
                    residuals_buffer.seek(0)
                    file_buffers.append(("ChiSquareResiduals", residuals_buffer))

                    heatmap_buffer = io.BytesIO()
                    plt.savefig(heatmap_buffer, format="png", dpi=300, bbox_inches='tight')
                    heatmap_buffer.seek(0)
                    file_buffers.append(("ResidualsHeatmap", heatmap_buffer))

                    st.download_button(
                        label="Download Contingency Table (CSV)",
                        data=contingency_buffer.getvalue(),
                        file_name=f"ContingencyTable_{random_number}.csv",
                        mime="text/csv",
                    )

                    st.download_button(
                        label="Download Chi-Square Residuals (CSV)",
                        data=residuals_buffer.getvalue(),
                        file_name=f"ChiSquareResiduals_{random_number}.csv",
                        mime="text/csv",
                    )

                    st.download_button(
                        label="Download Residuals Heatmap (PNG)",
                        data=heatmap_buffer,
                        file_name=f"ResidualsHeatmap_{random_number}.png",
                        mime="image/png",
                    )

                    zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                        for file_name, file_buffer in file_buffers:
                            if isinstance(file_buffer, io.StringIO):
                                zip_file.writestr(f"{file_name}_{random_number}.csv", file_buffer.getvalue())
                            elif isinstance(file_buffer, io.BytesIO):
                                zip_file.writestr(f"{file_name}_{random_number}.png", file_buffer.getvalue())
                    zip_buffer.seek(0)

                    st.download_button(
                        label="Download All as ZIP",
                        data=zip_buffer,
                        file_name=f"ChiSquareAnalysis_{random_number}.zip",
                        mime="application/zip",
                    )
                else:
                    st.write("There is no statistically significant association between clusters and groups (p >= 0.05).")

        with tab35:
            with st.form("Plot_35"):
                st.subheader("Normality test: Choose Analysis")

                analysis_choice = st.selectbox("Choose analysis by:", ["Group", "Leiden Clusters"])
                submit_button = st.form_submit_button("Run Analysis")


                if submit_button:
                    if analysis_choice == "Group":
                        mean_by_group = adata.to_df().groupby(adata.obs['Group']).mean()
                        median_by_group = adata.to_df().groupby(adata.obs['Group']).median()

                        st.write("Mean expression levels by Group:")
                        st.table(mean_by_group)

                        st.write("Median expression levels by Group:")
                        st.table(median_by_group)

                    elif analysis_choice == "Leiden Clusters":
                        mean_by_cluster = adata.to_df().groupby(adata.obs['leiden']).mean()
                        median_by_cluster = adata.to_df().groupby(adata.obs['leiden']).median()

                        st.write("Mean expression levels by Leiden Cluster:")
                        st.table(mean_by_cluster)

                        st.write("Median expression levels by Leiden Cluster:")
                        st.table(median_by_cluster)

                    st.write('--------')
                    st.subheader("Is the expression of Marker X normally distributed within each group and each cluster?")
                    st.write("Pro tip: When you run the test, you’ll get a W-statistic and a p-value. W-Statistic value ranges from 0 to 1. The closer it is to 1, the more likely it is that the data is normally distributed")

                    def normality_test(df, label):
                        results = []
                        for column in df.columns:
                            stat, p_value = stats.shapiro(df[column])
                            results.append({'Marker': column, 'Test': 'Shapiro-Wilk', 'Label': label, 'W-Statistic': stat, 'p-value': p_value})
                        return pd.DataFrame(results)

                    markers = adata.var_names

                    grouped_by_group = adata.to_df().groupby(adata.obs['Group'])
                    grouped_by_cluster = adata.to_df().groupby(adata.obs['leiden'])

                    normality_results = []

                    for group, group_data in grouped_by_group:
                        result = normality_test(group_data[markers], f"Group: {group}")
                        normality_results.append(result)

                    for cluster, cluster_data in grouped_by_cluster:
                        result = normality_test(cluster_data[markers], f"Cluster {cluster}")
                        normality_results.append(result)

                    merged_results = pd.concat(normality_results, ignore_index=True)
                    group_results = merged_results[merged_results['Label'].str.contains('Group')]
                    cluster_results = merged_results[merged_results['Label'].str.contains('Cluster')]

                    st.write("Normality Test Results by Group:")
                    st.dataframe(group_results)

                    st.write("Normality Test Results by Cluster:")
                    st.dataframe(cluster_results)

                    def plot_marker_distributions(df_grouped, label_type):
                        markers = df_grouped.head(1).columns
                        num_markers = len(markers)
                        num_cols = 4
                        num_rows = (num_markers // num_cols) + 1

                        fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
                        axes = axes.flatten()

                        for i, marker in enumerate(markers):
                            sns.histplot(df_grouped[marker], kde=True, ax=axes[i])
                            axes[i].set_title(f'{marker} Distribution in {label_type}')

                        plt.tight_layout()
                        st.pyplot(fig)

                    for group, group_data in grouped_by_group:
                        st.write(f"Marker Distribution in Group {group}")
                        with st.spinner('Working...'):
                            plot_marker_distributions(group_data[markers], f"Group {group}")
                            st.success("Done!")

                    merged_results['Test_Recommendation'] = merged_results['p-value'].apply(lambda p: 'T-test' if p > 0.05 else 'Mann-Whitney U test')

                    st.write("Normality Test Results with Test Recommendations:")
                    st.dataframe(merged_results)

                    st.write('Pro tip: In case of confusion between visual inspection and statistical test result, prefer Shapiro-Wilk test result as it is more sensitive.')

        with tab36:
            with st.form("Plot_36"):
                color_map = st.selectbox(
                    "Select Color Map",
                    options=["coolwarm", "viridis", "plasma", "inferno", "magma", "cividis", "Blues", "Reds", "Greens"],
                    index=0
                )

                fig_size = st.slider("Set Figure Size (Height and Width)", min_value=5, max_value=20, value=10)
                plot_file_format = st.selectbox("Select Plot File Format", options=["PNG", "PDF", "SVG", "JPEG"], index=0)
                submitted = st.form_submit_button("Run Correlation and Dendrogram Analysis for Each Group")

            if submitted:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                cell_type_column = "cell_type" if "cell_type" in adata.obs.columns else "leiden"
                unique_groups = adata.obs['Group'].unique()

                for group in unique_groups:
                    st.write(f"### Correlation Plot for Group: {group}")
                    adata_group = adata[adata.obs['Group'] == group]
                    adata_group.obs[cell_type_column] = adata_group.obs[cell_type_column].astype("category")
                    clustered_data = adata_group.to_df().groupby(adata_group.obs[cell_type_column]).mean()
                    corr_matrix = clustered_data.T.corr()
                    st.write(f"Plotting correlation heatmap with dendrogram for Group: {group}")

                    g = sns.clustermap(
                        corr_matrix,
                        method='average',
                        cmap=color_map,
                        row_cluster=True,
                        col_cluster=True,
                        figsize=(fig_size, fig_size)
                    )

                    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

                    buffer = io.BytesIO()
                    plt.savefig(buffer, format=plot_file_format.lower(), dpi=300, bbox_inches='tight')
                    buffer.seek(0)
                    plot_buffers.append((f"Correlation_Group_{group}", buffer))

                    st.pyplot(g)
                    plt.close(g.fig)

                for plot_name, plot_buffer in plot_buffers:
                    st.download_button(
                        label=f"Download {plot_name} ({plot_file_format.upper()})",
                        data=plot_buffer,
                        file_name=f"{plot_name}_{random_number}.{plot_file_format.lower()}",
                        mime=f"image/{plot_file_format.lower()}",
                    )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(f"{plot_name}_{random_number}.{plot_file_format.lower()}", plot_buffer.getvalue())
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Correlation Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"CorrelationPlots_{random_number}.zip",
                    mime="application/zip",
                )

        with tab37:
            with st.form("Plot_37"):
                color_map = st.selectbox(
                    "Select Color Map",
                    options=["coolwarm", "viridis", "plasma", "inferno", "magma", "cividis", "Blues", "Reds", "Greens"],
                    index=0
                )

                fig_size = st.slider("Set Figure Size (Height and Width)", min_value=5, max_value=20, value=10)
                plot_file_format = st.selectbox("Select Plot File Format", options=["PNG", "PDF", "SVG", "JPEG"], index=0)
                submitted = st.form_submit_button("Run Correlation and Dendrogram Analysis for Each Group")

            if submitted:
                random_number = random.randint(1000, 9999)
                plot_buffers = []

                cell_type_column = "cell_type" if "cell_type" in adata.obs.columns else "leiden"
                unique_groups = adata.obs['Group'].unique()

                for group in unique_groups:
                    st.write(f"### Correlation Plot for Group: {group}")
                    adata_group = adata[adata.obs['Group'] == group]
                    adata_group.obs[cell_type_column] = adata_group.obs[cell_type_column].astype("category")
                    clustered_data = adata_group.to_df().groupby(adata_group.obs[cell_type_column]).mean()
                    corr_matrix = clustered_data.T.corr()

                    g = sns.clustermap(
                        corr_matrix,
                        method='average',
                        cmap=color_map,
                        row_cluster=True,
                        col_cluster=True,
                        figsize=(fig_size, fig_size)
                    )

                    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

                    buffer = io.BytesIO()
                    plt.savefig(buffer, format=plot_file_format.lower(), dpi=300, bbox_inches='tight')
                    buffer.seek(0)
                    plot_buffers.append((f"Correlation_Group_{group}", buffer))

                    st.pyplot(g)
                    plt.close(g.fig)

                for plot_name, plot_buffer in plot_buffers:
                    st.download_button(
                        label=f"Download {plot_name} ({plot_file_format.upper()})",
                        data=plot_buffer,
                        file_name=f"{plot_name}_{random_number}.{plot_file_format.lower()}",
                        mime=f"image/{plot_file_format.lower()}",
                    )

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for plot_name, plot_buffer in plot_buffers:
                        zip_file.writestr(f"{plot_name}_{random_number}.{plot_file_format.lower()}", plot_buffer.getvalue())
                zip_buffer.seek(0)

                st.download_button(
                    label="Download All Correlation Plots as ZIP",
                    data=zip_buffer,
                    file_name=f"CorrelationPlots_{random_number}.zip",
                    mime="application/zip",
                )
