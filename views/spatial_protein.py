import streamlit as st
import urllib.request
from persist import persist


# IMG_REPO = 'https://raw.githubusercontent.com/matthewlu2/plot_spatial/main'
# IMG_REPO_2 = 'https://raw.githubusercontent.com/matthewlu2/plot_violin/main'

IMG_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main/spatial_protein'


st.markdown("<h2 style='text-align: center; color: black;'>Protein Expression Maps</h1>", unsafe_allow_html=True)
st.write("")

st.info("TBD")

file = open('text_files/common_protein_31.txt', 'r')
list = file.read().splitlines()

tabs_font_css = """
<style>
div[class*="stSelectbox"] label {
  color: purple;
}
</style>
"""
st.write(tabs_font_css, unsafe_allow_html=True)

a, b = st.columns(2)
df_sample = st.session_state.df_sample
sample_list = df_sample['Sample-ID'].values.tolist()
option = a.selectbox(
    label='Sample',
    options=sample_list,
    key = persist("sample_id")
    ) 

option2 = b.selectbox(
    'Protein',
    list)

def url_is_alive(url):
    """
    Checks that a given URL is reachable.
    :param url: A URL
    :rtype: bool
    """
    request = urllib.request.Request(url)
    request.get_method = lambda: 'HEAD'
    try:
        urllib.request.urlopen(request)
        return True
    except urllib.request.HTTPError:
        return False

image_na = "./logo/no_available_icon.png"

image_spatial_tissue = f"{IMG_REPO}/{option}/{option}.png"
# st.text(image_spatial)
if url_is_alive(image_spatial_tissue):
    a.image(image_spatial_tissue)
else:
    a.image(image_na)
# a.subheader('Spatial Plot')
image_spatial = f"{IMG_REPO}/{option}/{option2}.png"
# st.text(image_spatial)
if url_is_alive(image_spatial):
    b.image(image_spatial)
else:
    b.image(image_na)

