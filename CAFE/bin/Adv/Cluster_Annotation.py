import streamlit as st
import scanpy as sc
import re
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.cluster.hierarchy as sch
from scipy.sparse import issparse
from io import BytesIO
import random
import io
import zipfile
import tempfile

st.set_page_config(layout="centered")

image_path = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path)

image_path = os.path.join('bin', 'img', 'logo_v2.png')
st.image(image_path, caption='', use_container_width=True)

st.title("Advanced Annotation")
st.write('*The module allows annotating leiden clusters into cell types and save the adata file.*')

with st.form(key='upload_form'):
    uploaded_file = st.file_uploader("Upload AnnData (.h5ad) file", type="h5ad")
    submit_upload = st.form_submit_button("Load AnnData")

if 'adata' not in st.session_state:
    st.session_state.adata = None

if uploaded_file and submit_upload:
    st.session_state.adata = sc.read_h5ad(uploaded_file)
    st.write(f"AnnData loaded with shape: {st.session_state.adata.shape}")

    if 'X_umap' not in st.session_state.adata.obsm.keys():
        sc.tl.umap(st.session_state.adata)

    st.write("UMAP visualization with 'tab20c' colormap:")
    fig, ax = plt.subplots()

    sc.pl.umap(
        st.session_state.adata,
        color='leiden',
        palette="tab20c",
        ax=ax,
        show=False,
        legend_loc='on data'
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    st.pyplot(fig)


if st.session_state.adata is not None:
    st.write("Enter cell type annotations in the format: CD4 T cell = 1, CD8 T cell = 2, ...")

    cell_type_input = st.text_area("Cell Type Annotations")

    if st.button("Apply Cell Type Annotations"):
        try:
            cell_type_mapping = {}
            annotations = cell_type_input.split(",")
            for annotation in annotations:
                cell_type, cluster = annotation.split("=")
                cell_type_mapping[cluster.strip()] = cell_type.strip()

            st.session_state.adata.obs['cell_type'] = st.session_state.adata.obs['leiden'].map(cell_type_mapping)
            st.success("Cell type annotations applied successfully!")

            random_number = random.randint(1000, 9999)
            zip_buffer = io.BytesIO()

            with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
                    temp_path = temp_file.name
                    st.session_state.adata.write_h5ad(temp_path)
                    zip_file.write(temp_path, arcname=f"annotated_adata_{random_number}.h5ad")
                    os.remove(temp_path)

            zip_buffer.seek(0)

            st.download_button(
                label="Download Annotated Data as ZIP",
                data=zip_buffer,
                file_name=f"annotated_data_{random_number}.zip",
                mime="application/zip"
            )

            st.write("UMAP visualization with 'cell_type' annotations:")
            fig, ax = plt.subplots()

            sc.pl.umap(
                st.session_state.adata,
                color='cell_type',
                palette="tab20c",
                ax=ax,
                show=False,
                legend_loc='on data'
            )

            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])

            st.pyplot(fig)

            umap_buffer = io.BytesIO()
            fig.savefig(umap_buffer, dpi=300, format='png', bbox_inches='tight')
            umap_buffer.seek(0)

            st.download_button(
                label="Download UMAP as PNG",
                data=umap_buffer,
                file_name=f"UMAP_cell_type_{random_number}.png",
                mime="image/png"
            )

        except Exception as e:
            st.error(f"Error: {e}")
