## Initialized by amazing ChatGPT and Gemini, modified by human.

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import tempfile
import os
from typing import Optional, Tuple, List, Dict, Any

try:
    from scipy import sparse as sp
except Exception:
    sp = None

FIGSIZE = (4.8, 4.8)

IMAGE_NA_PATH = "./logo/no_available_icon.png"

st.markdown("<h2 style='text-align: center; color: black;'>View your uploaded result</h2>", unsafe_allow_html=True)
st.write("")

adata_in = st.session_state.get("adata_in", None)
adata_out = st.session_state.get("adata_out", None)
protein_names = st.session_state.get("protein_names", [])

if (adata_in is None) or (adata_out is None):
    st.warning("No uploaded result found. Please go to **Upload Data** and run again.")
    st.stop()

c1, c2 = st.columns([3, 1])
with c1:
    st.subheader('Sample name: "Uploaded data"')
with c2:
    try:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            adata_out.write_h5ad(tmp.name)
            tmp_path = tmp.name
        with open(tmp_path, "rb") as f:
            st.download_button(
                label="⬇️ Download adata_out (.h5ad)",
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
                # 用第一个 key 兜底
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
                figsize=FIGSIZE
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
    #
    # if library_id and img_key:
    #     try:
    #         img = adata.uns["spatial"][library_id]["images"][img_key]
    #         if img is not None and hasattr(img, 'shape'):
    #             fig, ax = plt.subplots(figsize=FIGSIZE)
    #             ax.imshow(img)
    #             ax.axis("off")
    #             plt.tight_layout()
    #             return fig
    #     except Exception as e:
    #         pass

    if library_id and img_key:
        try:
            fig, ax = plt.subplots(figsize=FIGSIZE)
            sc.pl.spatial(
                adata,
                color=None,
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax
            )
            return fig
        except Exception:
            plt.close(fig)

    try:
        fig, ax = plt.subplots(figsize=FIGSIZE)
        sc.pl.spatial(
            adata,
            color=None,
            show=False,
            ax=ax
        )
        return fig
    except Exception:
        plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axis("off")
    ax.text(0.5, 0.5, "No tissue image", ha="center", va="center", fontsize=14, transform=ax.transAxes)
    return fig


def _plot_spatial_expr(adata, gene: Optional[str], library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:

    if gene is None or str(gene) not in _varnames(adata):
        return _plot_scatter_expr(adata, gene)

    if library_id and img_key:
        try:
            fig, ax = plt.subplots(figsize=FIGSIZE)
            sc.pl.spatial(
                adata,
                color=str(gene),
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax,
                size = 1.7
            )
            return fig
        except Exception:
            plt.close(fig)

    return _plot_scatter_expr(adata, gene)


has_img_out, lib_id_out, img_key_out, spatial_meta_out = _probe_spatial_meta(adata_out)
has_img_in,  lib_id_in,  img_key_in,  spatial_meta_in  = _probe_spatial_meta(adata_in)
# #
# with st.expander("Debug (click to expand)"):
#     st.write("adata_out.uns keys:", list(getattr(adata_out, "uns_keys", lambda: adata_out.uns.keys())()))
#     st.write("adata_out.obsm keys:", list(getattr(adata_out, "obsm_keys", lambda: adata_out.obsm.keys())()))
#     st.write("Detected (out) library_id:", lib_id_out, "img_key:", img_key_out)
#     if _has_spatial_coords(adata_out):
#         st.write("obsm['spatial'] shape (adata_out):", adata_out.obsm["spatial"].shape)
#
#     st.write("---")
#     st.write("adata_in.uns keys:", list(getattr(adata_in, "uns_keys", lambda: adata_in.uns.keys())()))
#     st.write("adata_in.obsm keys:", list(getattr(adata_in, "obsm_keys", lambda: adata_in.obsm.keys())()))
#     st.write("Detected (in) library_id:", lib_id_in, "img_key:", img_key_in)
#     if _has_spatial_coords(adata_in):
#         st.write("obsm['spatial'] shape (adata_in):", adata_in.obsm["spatial"].shape)
#
#     st.write("---")
#     st.write("adata_out.var_names (first 10):", _varnames(adata_out)[:10])
#
# with st.expander("Debug - Image Data Check"):
#     st.write("Checking actual image data...")
#     try:
#         spatial_dict = adata_out.uns["spatial"]
#         lib_dict = spatial_dict["CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer"]
#         st.write("Library dict keys:", list(lib_dict.keys()))
#
#         if "images" in lib_dict:
#             images_dict = lib_dict["images"]
#             st.write("Images dict keys:", list(images_dict.keys()))
#
#             if "hires" in images_dict:
#                 hires_img = images_dict["hires"]
#                 st.write("Hires image type:", type(hires_img))
#                 st.write("Hires image shape:", getattr(hires_img, "shape", "No shape attribute"))
#                 st.write("Hires image dtype:", getattr(hires_img, "dtype", "No dtype attribute"))
#
#                 if hires_img is not None:
#                     fig_test = plt.figure(figsize=(5, 5))
#                     plt.imshow(hires_img)
#                     plt.axis("off")
#                     plt.title("Direct image display test")
#                     st.pyplot(fig_test)
#                     plt.close(fig_test)
#             else:
#                 st.error("'hires' key not found in images dict")
#         else:
#             st.error("'images' key not found in library dict")
#
#         if "scalefactors" in lib_dict:
#             st.write("Scalefactors:", lib_dict["scalefactors"])
#
#     except Exception as e:
#         st.error(f"Error checking image: {e}")
#         import traceback
#
#         st.code(traceback.format_exc())
#

if len(protein_names) == 0:
    st.info("No proteins found in adata_out.var_names")
    gene = None
else:
    gene = st.selectbox("Select a protein", options=protein_names, index=0)

st.divider()

g1, g2, g3 = st.columns(3)

with g1:
    st.caption("Tissue image (adata_out)")
    fig1 = _plot_tissue_only(adata_out, lib_id_out if has_img_out else None, img_key_out if has_img_out else None)
    st.pyplot(fig1, width="stretch")
    plt.close(fig1)

with g2:
    st.caption("Spatial protein expression (imputed)")
    fig2 = _plot_spatial_expr(adata_out, gene, lib_id_out if has_img_out else None, img_key_out if has_img_out else None)
    st.pyplot(fig2, width="stretch")
    plt.close(fig2)

with g3:
    st.caption("Spatial mRNA expression (your data)")
    if gene is None or str(gene) not in _varnames(adata_in):
        fig3 = _plot_image_placeholder(IMAGE_NA_PATH)
    else:
        fig3 = _plot_spatial_expr(adata_in, gene, lib_id_in if has_img_in else None, img_key_in if has_img_in else None)
    st.pyplot(fig3, width="stretch")
    plt.close(fig3)
