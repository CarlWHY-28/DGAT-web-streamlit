import streamlit as st
from style import page_style, footer
import pandas as pd
from persist import load_widget_state, persist
from views.utils import get_sample_dataframe
from style import define_layout
# st.cache_data.clear()
# st.cache_resource.clear()

import time, os, psutil
# --- PAGE SETUP ----
process = psutil.Process(os.getpid())



st.set_page_config(
        page_title='DGAT',
        page_icon= "./logo/gbm_ribbon.png",
        initial_sidebar_state="expanded",
)
if "has_upload" not in st.session_state:
    st.session_state["has_upload"] = False

define_layout(max_width='80%', padding_top='2rem', padding_right='0rem', padding_left='0rem', padding_bottom='0rem')

load_widget_state()

# df_sample = pd.read_csv('./data/dataset.csv')
# df_sample.index = df_sample.index + 1
# st.session_state['df_sample'] = df_sample
df_sample = get_sample_dataframe('./data/dataset.csv')
st.session_state['df_sample'] = df_sample
persist("sample_id")

emoji = "🔹"  # "🔸" "💠" ...

home_page = st.Page(
    page = "views/home.py",
    title = "Home",
    icon = emoji,
    default= True,
)

datasets_page = st.Page(
    page = "views/dataset.py",
    title = "Dataset Explorer",
    icon = emoji
)

gene_page = st.Page(
    page = "views/spatial_protein.py",
    title = "Protein Expression Maps",
    icon = emoji
)

s_tf_page = st.Page(
    page = "views/spatial_tf.py",
    title = "TF Activity Maps",
    icon = emoji
)

s_pathway_page = st.Page(
    page = "views/spatial_pathway.py",
    title = "Pathway Activity Maps",
    icon = emoji
)

contact_page = st.Page(
    page = "views/contact.py",
    title = "Contact us",
    icon = emoji
)

citation_page = st.Page(
    page = "views/citation.py",
    title = "Citation",
    icon = emoji
)

upload_page = st.Page(
    page = "views/upload.py",
    title = "Upload Data",
    icon = emoji
)

view_uploaded_page = st.Page(
    page = "views/view_uploaded.py",
    title = "View your data",
    icon = emoji
)

# -- NAVIGATION --
nav_groups = {
    "": [home_page, datasets_page],
    "Impute Your Data": [upload_page],  # Upload
   # "Analysis of Individual Samples": [gene_page, s_tf_page, s_pathway_page],
    #"Comparison Across Samples": [],
    "Resources": [citation_page, contact_page]

}
if st.session_state.get("has_upload", False):
    nav_groups["Impute Your Data"].append(view_uploaded_page)

pg = st.navigation(nav_groups)
pg.run()

st.divider()
st.markdown(footer, unsafe_allow_html=True)

while True:
    print(f"Memory: {process.memory_info().rss / 1024**2:.2f} MB", end="\r")
    time.sleep(1)