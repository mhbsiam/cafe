"""Shared landing-page content for the desktop and hosted CAFE builds.

Keeping the content in one place means the two entry points (``App/CAFE.py`` and
``web/CAFE.py``) stay in sync. The only differences are the workflow image and
whether the system-requirements section is shown.
"""
import os

import pandas as pd
import streamlit as st

from theme import info_card


_EXAMPLE_CSV_DATA = {
    "CD27": [27257.3, 16917.6, 47378.3, 1769.79, 40472.1, 487.425, 31111.4, -249.461, 18263.4],
    "CD38": [5442.59, 7401.0, 4412.47, 54090.1, 2639.16, 23893.5, 2010.18, 32396.4, 455.373],
    "CD127": [13068.2, 10038.9, 19117.0, 265.372, 20865.8, 2089.89, 3647.74, 1188.07, 2590.04],
    "HLA-DR": [-2432.2, -1056.87, -1147.41, 8136.79, -2138.81, 6163.82, -1939.37, 18938.1, -1000.35],
    "CD1c": [-3011.07, -743.682, -654.677, -2202.24, -1856.0, 555.821, 264.905, 455.37, -2214.69],
    "CD141": [2048.14, -3049.75, -5256.23, 3329.33, -3293.17, 3007.88, 1963.75, -3065.3, 2975.21],
    "CD45ra": [-516.212, -542.594, 68264.9, 49450.1, 81219.1, 69839.0, 12751.9, 17207.5, 421.45],
    "CD16": [10914.9, 3165.02, 17845.1, 3431.71, 10629.2, 7129.61, 8712.62, 1354.81, 5342.13],
    "CCR5": [4355.31, 4260.19, 8269.94, 1735.02, 4249.25, 1701.63, 4819.16, 1981.51, 2328.19],
    "CD4": [87660.1, 90922.7, 97427.5, 908.1, 91438.0, 164.528, 89793.9, -38.6853, 54503.7],
    "CD11c": [2256.65, 2309.36, 5624.15, 15106.8, 3298.7, 19372.5, 5114.97, 9176.99, 2726.56],
    "CD56": [1316.62, -686.368, 2726.43, 6331.52, 960.926, 9389.31, 1778.75, 17000.4, 1764.66],
}


def render_landing_page(
    workflow_image="workflow.png",
    how_to_title="How to Run the App",
    include_system_requirements=True,
):
    """Render the CAFE welcome page.

    Parameters
    ----------
    workflow_image : str
        Filename (under ``bin/img/``) for the workflow diagram.
    how_to_title : str
        Heading for the workflow section ("How to Run the App" for desktop,
        "How to use CAFE" for hosted).
    include_system_requirements : bool
        Whether to show the desktop-specific system-requirements image.
    """
    st.image(os.path.join("bin", "img", "logo.png"), caption="", width="stretch")
    st.logo(os.path.join("bin", "img", "s_logo.png"))

    st.title("Cell Analyzer for Flow Experiment")
    st.subheader("Explore and analyze flow cytometry data interactively")

    st.markdown(
        "This tool analyzes immune cells from spectral flow cytometry and generates "
        "publishable figures. For detailed documentation, visit the "
        "[GitHub page](https://github.com/mhbsiam/cafe)."
    )
    st.markdown("---")

    st.markdown("## Citation")
    st.markdown(
        """
        If you use CAFE in your research, please cite our papers:

        > Md Hasanul Banna Siam, Md Akkas Ali, Donald Vardaman, Mallikarjun Patil, Satwik Acharyya, and Daniel Tyrrell. “Cafe: An Integrated Web App for High-Dimensional Analysis and Visualization in Spectral Flow Cytometry,” *bioRxiv* (2024). https://doi.org/10.1101/2024.12.03.626714.

        > Donald Vardaman, Md Akkas Ali, Md Hasanul Banna Siam, Chase Bolding, Harrison Tidwell, Holly R. Stephens, Mallikarjun Patil, and Daniel J. Tyrrell. "Development of a Spectral Flow Cytometry Analysis Pipeline for High-Dimensional Immune Cell Characterization." *The Journal of Immunology* (2024). https://doi.org/10.4049/jimmunol.2400370.
        """
    )

    st.markdown("---")

    st.subheader(how_to_title)
    st.image(
        os.path.join("bin", "img", workflow_image),
        caption="CAFE workflow: upload CSVs, process into an AnnData object, then visualize and export results.",
        width="stretch",
    )

    st.subheader("1. Preparing the CSV files")
    st.markdown(
        "Open FlowJo or similar software and upload your FCS files. After manual "
        "gating to remove debris, dead cells, and doublets, export scaled CSV files. "
        "Gating on a cell type such as CD45+ will give more focused clustering results."
    )

    info_card(
        title="File naming convention",
        body="Name each CSV as <code>SampleID_Group.csv</code>, for example "
             "<code>ABC001_Aged.csv</code> and <code>ABC002_Young.csv</code>. "
             "CAFE expects two groups to run statistics properly.",
        kind="info",
    )

    st.markdown("Ensure that every CSV file:")
    st.markdown(
        "- Contains the same number of columns (markers)  \n"
        "- Uses identical marker names across files  \n"
        "- Uses unique sample names across files  \n"
        "- Example markers: CD45, TCRY, CD28, etc.  \n"
        "- Mismatched column names will cause errors."
    )

    st.markdown("**CSV file structure:**")
    st.markdown(
        "Marker names as column headers and numerical values representing expression levels."
    )
    st.dataframe(pd.DataFrame(_EXAMPLE_CSV_DATA), use_container_width=True)

    st.subheader("2. Load CSV files")
    st.markdown("Select all CSV files together and upload them on the Data Processing page.")

    st.subheader("3. Follow the on-screen steps")
    st.markdown(
        "Use the guided pipeline to run quality control, dimensionality reduction, "
        "clustering, and statistical comparisons."
    )

    st.subheader("4. Download the results")
    st.markdown(
        "Download the ZIP file containing publication-ready figures, tables, and the "
        "processed AnnData object."
    )

    st.subheader("5. Visualize the data")
    st.markdown(
        "Upload the saved <code>adata_.h5ad</code> file to the Visualization page for "
        "interactive exploration and additional plots.",
        unsafe_allow_html=True,
    )

    info_card(
        title="Quick restart",
        body="To restart the program, clear the cache from the top-right settings menu and reload the page.",
        kind="info",
    )

    st.divider()

    if include_system_requirements:
        st.subheader("System Requirements")
        st.image(
            os.path.join("bin", "img", "os_selection.png"),
            caption="Choosing the right CAFE build for your operating system.",
            width="stretch",
        )
        st.divider()

    st.image(os.path.join("bin", "img", "funding.png"), caption="", width="stretch")
    st.image(os.path.join("bin", "img", "uab.png"), caption="", width="stretch")
