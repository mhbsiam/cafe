import streamlit as st

from theme import hero_entrance, landing_ambient
from App.landing_content import render_landing_page

st.set_page_config(layout="centered")

landing_ambient()
hero_entrance()

render_landing_page(
    workflow_image="workflow_w.png",
    how_to_title="How to use CAFE",
    include_system_requirements=False,
)
