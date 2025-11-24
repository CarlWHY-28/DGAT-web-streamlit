# Code borrowed/modified from:
# (1) https://github.com/CarlWHY-28/DGAT-web-streamlit/blob/master/views/view_uploaded.py
# (2) https://github.com/osmanbeyoglulab/DGAT/blob/master/Reproduction/Figure_3AB_4AB_5AB_5DE_SupF2AB_SupF3AB.ipynb

import streamlit as st
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

adata_in = st.session_state.get("adata_in", None)
adata_out = st.session_state.get("adata_out", None)

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Leiden clustering</h2>", unsafe_allow_html=True)
st.write("")

def _plot_leiden_clustering(
    adata, 
    ax, 
    n_neighbors = 10, 
    resolution = 0.5, 
    title = None
    ):
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, use_rep = 'X')
    sc.tl.leiden(adata, resolution = resolution)
    sq.pl.spatial_scatter(adata, color = "leiden", title = title, ax = ax)

    ax.get_legend().set_title("Leiden cluster")

g1, g2 = st.columns(2)

with g1:
    fig_protein, ax_protein = plt.subplots()
    _plot_leiden_clustering(adata_out, title = "Protein", ax = ax_protein)
    st.pyplot(fig_protein, use_container_width = True)
    plt.close(fig_protein)

with g2:
    fig_mRNA, ax_mRNA = plt.subplots()
    _plot_leiden_clustering(adata_in, title = "mRNA", ax = ax_mRNA)
    st.pyplot(fig_mRNA, use_container_width = True)
    plt.close(fig_mRNA)