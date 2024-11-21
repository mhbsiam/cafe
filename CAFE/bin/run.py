import streamlit as st

page1 = st.Page("App/CAFE.py", default=True)
page2 = st.Page("Tools/Data_Processing.py")
page3 = st.Page("Tools/Visualization.py")

page4 = st.Page("Adv/Merge_Clusters.py")
page5 = st.Page("Adv/Cluster_Annotation.py")

page6 = st.Page("Adv/Cluster_Evaluation.py")
page7 = st.Page("Adv/Selective_Clustering.py")
page8 = st.Page("Adv/Downsampling.py")

pg = st.navigation(
    {
        "Welcome": [page1],
        "Tools": [page2, page3],
        "Advanced": [page4, page5, page6, page7, page8],
    }
)

pg.run()
