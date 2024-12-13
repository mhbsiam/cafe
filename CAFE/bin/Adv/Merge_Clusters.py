import streamlit as st
import scanpy as sc
import re
import matplotlib.pyplot as plt
import os
import zipfile
import tempfile
import random
import io

st.set_page_config(layout="centered")

image_path = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path)

# Display an image
image_path = os.path.join('bin', 'img', 'logo_v2.png')
st.image(image_path, caption='', use_container_width=True)

st.title("Merge Clusters")
st.write('*This module allows users to merge subclusters into metaclusters and save the new adata file.*')

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

if st.session_state.adata is not None:
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

    st.divider()
    image_path = os.path.join('bin', 'img', 'merging.png')
    st.image(image_path, caption='', use_container_width=True)
    st.divider()
    with st.form(key='merge_form'):
        merge_input = st.text_input(
            "Enter clusters to merge in the format {new_cluster_name:old_cluster,old_cluster}. For example: {1:7,3,8}, {2:5,9,4}",
            value="{1:7,3,8}, {2:5,9,4}"
        )
        submit_merge = st.form_submit_button("Merge Clusters")

    if merge_input and submit_merge:
        progress_bar = st.progress(0)
        matches = re.findall(r'{(\d+):([\d,]+)}', merge_input)
        leiden_copy = st.session_state.adata.obs['leiden'].astype(str).copy()

        cluster_mapping = {}
        for match in matches:
            new_cluster = match[0]
            clusters_to_merge = match[1].split(',')
            clusters_to_merge = [x.strip() for x in clusters_to_merge]
            for cluster in clusters_to_merge:
                cluster_mapping[cluster] = new_cluster

        leiden_copy.replace(cluster_mapping, inplace=True)
        st.session_state.adata.obs['leiden'] = leiden_copy

        st.write("Leiden clusters updated successfully.")
        progress_bar.progress(50)

        st.write("Updated UMAP visualization:")
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
        progress_bar.progress(75)

        random_number = random.randint(1000, 9999)
        merge_log = f"Cluster Merging Input:\n{merge_input}"

        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w") as zip_file:
            # Add the merged AnnData file to the zip
            with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
                temp_path = temp_file.name
                st.session_state.adata.write_h5ad(temp_path)
                zip_file.write(temp_path, arcname=f"merged_adata_{random_number}.h5ad")

            # Add the merge log to the zip
            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as log_temp_file:
                log_temp_path = log_temp_file.name
                with open(log_temp_path, "w") as log_file:
                    log_file.write(merge_log)
                zip_file.write(log_temp_path, arcname=f"merge_log_{random_number}.txt")

        zip_buffer.seek(0)
        progress_bar.progress(100)

        st.download_button(
            label=f"Download Merged Data and Log (merged_data_{random_number}.zip)",
            data=zip_buffer,
            file_name=f"merged_data_{random_number}.zip",
            mime="application/zip",
        )

        st.success("Merged data and log files are ready for download.")
