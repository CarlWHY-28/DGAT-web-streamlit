## Initialized by amazing ChatGPT and Gemini, modified by human.

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import tempfile
import os
from typing import Optional, Tuple, List, Dict, Any
import sys
import psutil
import squidpy as sq
import random

try:
    from scipy import sparse as sp
except Exception:
    sp = None



FIGSIZE = (4.8, 4.8)

IMAGE_NA_PATH = "./logo/no_available_icon.png"

st.markdown("<h2 style='text-align: center; color: black;'>View your imputed result</h2>", unsafe_allow_html=True)
st.write("")

adata_in = st.session_state.get("adata_in", None)
adata_out = st.session_state.get("adata_out", None)
protein_names = st.session_state.get("protein_names", [])

if (adata_in is None) or (adata_out is None):
    st.warning("No uploaded result found. Please go to **Upload Data** and run again.")
    st.stop()

c1, c2 = st.columns([3, 1])
with c1:
    st.subheader('Protein imputed data preview')
with c2:
    try:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            adata_out.write_h5ad(tmp.name)
            tmp_path = tmp.name
        with open(tmp_path, "rb") as f:
            st.download_button(
                label="⬇️ Download predicted data (.h5ad)",
                data=f.read(),
                file_name="adata_out.h5ad",
                mime="application/octet-stream"
            )
        os.remove(tmp_path)
    except Exception as e:
        st.error(f"Failed to prepare download: {e}")

def _varnames(adata) -> List[str]:
    return [str(x) for x in getattr(adata, "var_names", [])]

def _to_1d_vals(adata, gene) -> np.ndarray:
    X = adata[:, str(gene)].X
    if sp is not None and sp.issparse(X):
        return np.ravel(X.A)
    return np.asarray(X).ravel()

def _has_spatial_coords(adata) -> bool:
    try:
        return "spatial" in adata.obsm and getattr(adata.obsm["spatial"], "shape", (0, 0))[1] >= 2
    except Exception:
        return False

def _probe_spatial_meta(adata) -> Tuple[bool, Optional[str], Optional[str], Dict[str, Any]]:

    meta = {}
    try:
        if "spatial" not in adata.uns or not isinstance(adata.uns["spatial"], dict):
            return False, None, None, meta

        spatial_uns = adata.uns["spatial"]
        libs = list(spatial_uns.keys())
        if len(libs) == 0:
            return False, None, None, meta
        lib = libs[0]
        lib_dict = spatial_uns.get(lib, {})
        images_dict = lib_dict.get("images", {})

        candidates = ["hires", "image", "lowres"]
        img_key = None
        if isinstance(images_dict, dict) and len(images_dict) > 0:
            for k in candidates:
                if k in images_dict:
                    img_key = k
                    break
            if img_key is None:
                img_key = list(images_dict.keys())[0]
            return True, lib, img_key, {"libs": libs, "img_keys": list(images_dict.keys())}
        for k in candidates:
            if k in lib_dict:
                return True, lib, k, {"libs": libs, "img_keys": [k]}

        return False, lib, None, {"libs": libs, "img_keys": []}
    except Exception:
        return False, None, None, meta


def _plot_spatial_tissue_scanpy(adata, library_id: Optional[str], img_key: Optional[str]) -> Optional[plt.Figure]:

    if library_id is not None and img_key is not None:
        try:
            fig_obj = sc.pl.spatial(
                adata,
                color=None,
                library_id=library_id,
                img_key=img_key,
                show=False,
                return_fig=True,
                figsize=FIGSIZE,

            )
            return fig_obj if fig_obj is not None else plt.gcf()
        except Exception:
            pass

    try:
        fig_obj = sc.pl.spatial(
            adata,
            color=None,
            show=False,
            return_fig=True,
            figsize=FIGSIZE
        )
        return fig_obj if fig_obj is not None else plt.gcf()
    except Exception:
        return None


def _plot_spatial_expr_scanpy(adata, gene: str, library_id: Optional[str], img_key: Optional[str]) -> Optional[plt.Figure]:
    if library_id is None or img_key is None:
        return None
    try:
        fig_obj = sc.pl.spatial(
            adata,
            color=str(gene),
            library_id=library_id,
            img_key=img_key,
            show=False,
            return_fig=True,
            figsize=FIGSIZE
        )
        fig = fig_obj if fig_obj is not None else plt.gcf()
        return fig
    except Exception:
        return None

def _plot_scatter_expr(adata, gene: Optional[str]) -> plt.Figure:
    fig = plt.figure(figsize=FIGSIZE)
    if not _has_spatial_coords(adata) or gene is None or str(gene) not in _varnames(adata):
        plt.axis("off")
        plt.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16)
        return fig

    coords = adata.obsm["spatial"]
    vals = _to_1d_vals(adata, gene)
    sca = plt.scatter(coords[:, 0], coords[:, 1], s=20, c=vals)
    #plt.gca().invert_yaxis()
    plt.xticks([]); plt.yticks([])
    plt.colorbar(sca, shrink=0.75).set_label(str(gene))
    return fig

def _plot_image_placeholder(img_path: str) -> plt.Figure:
    fig = plt.figure(figsize=FIGSIZE)
    try:
        img = plt.imread(img_path)
        plt.imshow(img)
        plt.axis("off")
    except Exception:
        plt.axis("off")
        plt.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16)
    return fig


def _plot_tissue_only(adata, library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    try:
        if library_id and img_key:
            sc.pl.spatial(
                adata,
                color=None,
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax
            )
        else:
            sc.pl.spatial(adata, color=None, show=False, ax=ax)
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
        ax.set_aspect('equal', adjustable='box')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        return fig
    except Exception:
        plt.close(fig)

    # fallback: NA
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    ax.axis("off")
    ax.text(0.5, 0.5, "No tissue image", ha="center", va="center", fontsize=14, transform=ax.transAxes)
    return fig

def _plot_spatial_expr_mrna(adata, gene: Optional[str], library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    if gene is None or str(gene) not in _varnames(adata):
        ax.axis("off")
        ax.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16, transform=ax.transAxes)
        return fig

    if library_id and img_key:
        try:
            sc.pl.spatial(
                adata,
                color=str(gene),
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax,
                size=1.7,
                cmap="viridis"   # ⭐ 新增：mRNA 用 viridis
            )
            ax.set_xlim(ax.get_xlim())
            ax.set_ylim(ax.get_ylim())
            ax.set_aspect('equal', adjustable='box')
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            return fig
        except Exception:
            plt.close(fig)

    return _plot_scatter_expr(adata, gene)


def _plot_spatial_expr(adata, gene: Optional[str], library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    if gene is None or str(gene) not in _varnames(adata):
        ax.axis("off")
        ax.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16, transform=ax.transAxes)
        return fig

    if library_id and img_key:
        try:
            sc.pl.spatial(
                adata,
                color=str(gene),
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax,
                size=1.7,
                cmap="plasma"   # ⭐ 新增：protein 用 plasma
            )
            ax.set_xlim(ax.get_xlim())
            ax.set_ylim(ax.get_ylim())
            ax.set_aspect('equal', adjustable='box')
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            return fig
        except Exception:
            plt.close(fig)

    # fallback scatter
    return _plot_scatter_expr(adata, gene)

def _plot_leiden_clustering(
    adata, 
    ax, 
    n_neighbors = 10, 
    resolution = 0.5, 
    title = None, 
    seed = 0
    ):
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, use_rep = 'X', random_state = seed)
    sc.tl.leiden(adata, resolution = resolution, random_state = seed)
    if 'leiden_colors' in adata.uns.keys():
        adata.uns.pop('leiden_colors')
    sq.pl.spatial_scatter(adata, color = "leiden", title = title, ax = ax)

    ax.get_legend().set_title("Leiden cluster")

has_img_out, lib_id_out, img_key_out, spatial_meta_out = _probe_spatial_meta(adata_out)
has_img_in,  lib_id_in,  img_key_in,  spatial_meta_in  = _probe_spatial_meta(adata_in)

if len(protein_names) == 0:
    st.info("No proteins found in adata_out.var_names")
    gene = None
else:
    gene = st.selectbox("Select a protein", options=protein_names, index=0)

st.divider()


with st.expander("ℹ️ Data Preprocessing & Missing Value Handling"):
    st.markdown("""
    **mRNA Preprocessing Workflow:**
    The Spatial Transcriptomics (ST) data underwent the following preprocessing steps (`preprocess_ST`):
    1.  **Normalization:** Total counts per cell were normalized to a target sum of **10,000**.
    2.  **Log Transformation:** Data was transformed using `log1p` (natural logarithm of 1 + x).
    3.  **Scaling:** Data was scaled to unit variance, with values clipped to a maximum of **10**.

    **Handling Missing Genes:**
    * **Imputation of Zeros:** If a gene (mRNA) is present in the protein prediction target list but **missing** in the original ST sample, it is automatically filled with **0** for calculation and visualization purposes.
    """)



g1, g2, g3 = st.columns(3)

with g1:
    st.caption("Tissue image (adata_out)")
    fig1 = _plot_tissue_only(adata_out, lib_id_out if has_img_out else None, img_key_out if has_img_out else None)
    st.pyplot(fig1, use_container_width=True)
    plt.close(fig1)

with g2:
    st.caption("Spatial Protein Expression (Imputed)")
    fig2 = _plot_spatial_expr(adata_out, gene, lib_id_out if has_img_out else None, img_key_out if has_img_out else None)
    st.pyplot(fig2, use_container_width=True)
    plt.close(fig2)


with g3:
    st.caption("Spatial mRNA Expression (Your Data)")
    if gene is None or str(gene) not in _varnames(adata_in):
        fig3 = _plot_image_placeholder(IMAGE_NA_PATH)
    else:
        fig3 = _plot_spatial_expr_mrna(adata_in, gene, lib_id_in if has_img_in else None,
                                       img_key_in if has_img_in else None)
    st.pyplot(fig3, use_container_width=True)
    plt.close(fig3)
#
#
# process = psutil.Process(os.getpid())
# mem_info = process.memory_info()

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Leiden clustering</h2>", unsafe_allow_html=True)
st.write("")

g1_b, g2_b = st.columns(2)

with g1_b:
    #resolution_protein = st.slider("Select a resolution.", min_value = 0.2, max_value = 2.0, value = 0.5, step = 0.1)
    #n_neighbors_protein = st.slider("Select the number of neighbors.", min_value = 10, max_value = 100, value = 10, step = 10, key = "n_neighbors_protein_slider")
    resolution_protein = st.number_input("Resolution:", min_value = 0.2, max_value = 2.0, value = 0.5, step = 0.1)
    n_neighbors_protein = st.number_input("Number of neighbors:", min_value = 10, max_value = 100, value = 10, step = 10, key = "n_neighbors_protein_input")
    
    fig_protein, ax_protein = plt.subplots()
    _plot_leiden_clustering(adata_out, ax = ax_protein, n_neighbors = n_neighbors_protein, resolution = resolution_protein, title = "Protein")
    st.pyplot(fig_protein, use_container_width = True)
    plt.close(fig_protein)

with g2_b:

    resolution_mRNA = st.number_input("Resolution:", min_value = 0.2, max_value = 2.0, value = 1.0, step = 0.1)
    n_neighbors_mRNA = st.number_input("Number of neighbors:", min_value = 10, max_value = 100, value = 10, step = 10, key = "n_neighbors_mRNA_input")
    
    fig_mRNA, ax_mRNA = plt.subplots()
    _plot_leiden_clustering(adata_in, ax = ax_mRNA, n_neighbors = n_neighbors_mRNA, resolution = resolution_mRNA, title = "mRNA")
    st.pyplot(fig_mRNA, use_container_width = True)
    plt.close(fig_mRNA)
