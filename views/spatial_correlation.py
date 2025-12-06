import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import tempfile
import os
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import sys
import psutil
from dgat_utils.downstream import compute_moran_single
from scipy import sparse as sp

FIGSIZE = (4.8, 4.8)
IMAGE_NA_PATH = "./logo/no_available_icon.png"

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Correlation Analysis</h2>", unsafe_allow_html=True)
st.write("")

adata_in = st.session_state.get("adata_in", None)
adata_out = st.session_state.get("adata_out", None)
protein_names = st.session_state.get("protein_names", [])

if (adata_in is None) or (adata_out is None):
    st.warning("No uploaded result found. Please go to **Upload Data** and run again.")
    st.stop()

if "moran_result_df" not in st.session_state:
    with st.spinner("Calculating Moran's I for Protein and mRNA (this may take a moment)..."):
        try:
            moran_df, moran_stats = compute_moran_single(
                adata_out,
                adata_in,
                tissue_name="Current Sample",
                coord_type="grid",
                n_perms=10
            )
            st.session_state["moran_result_df"] = moran_df
            st.session_state["moran_stats"] = moran_stats
        except Exception as e:
            st.error(f"Error computing Moran's I: {e}")
            st.session_state["moran_result_df"] = pd.DataFrame()
            st.session_state["moran_stats"] = {"stat": np.nan, "pvalue": np.nan}

moran_df = st.session_state["moran_result_df"]
moran_stats = st.session_state["moran_stats"]

st.subheader("1. Global Spatial Autocorrelation Analysis")

m1, m2 = st.columns(2)
with m1:
    st.metric("Wilcoxon Statistic", value=f"{moran_stats.get('stat', np.nan):.4f}")
with m2:
    pval = moran_stats.get('pvalue', np.nan)
    st.metric("P-value", value=f"{pval:.4e}" if pval is not None else "NaN")

st.markdown("**Moran's I per gene:**")
if not moran_df.empty:
    display_df = moran_df.copy()
    st.dataframe(display_df, use_container_width=True, height=300)

    csv = display_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="⬇️ Download Moran's I Table (.csv)",
        data=csv,
        file_name='moran_statistics.csv',
        mime='text/csv',
    )
else:
    st.warning("Moran's I calculation failed or returned empty result.")

st.divider()

c1, c2 = st.columns([3, 1])
with c1:
    st.subheader('2. Spatial Expression Visualization')
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
    try:
        X = adata[:, str(gene)].X
        if sp is not None and sp.issparse(X):
            return np.ravel(X.A)
        return np.asarray(X).ravel()
    except KeyError:
        return np.zeros(adata.n_obs)


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


def _plot_tissue_only(adata, library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    try:
        if library_id and img_key:
            sc.pl.spatial(adata, color=None, library_id=library_id, img_key=img_key, show=False, ax=ax)
        else:
            sc.pl.spatial(adata, color=None, show=False, ax=ax)
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
        ax.set_aspect('equal', adjustable='box')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        return fig
    except Exception:
        plt.close(fig)
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    ax.axis("off")
    ax.text(0.5, 0.5, "No tissue image", ha="center", va="center", fontsize=14, transform=ax.transAxes)
    return fig


def _plot_scatter_expr(adata, gene: Optional[str]) -> plt.Figure:
    fig = plt.figure(figsize=FIGSIZE)
    if not _has_spatial_coords(adata) or gene is None or str(gene) not in _varnames(adata):
        plt.axis("off")
        plt.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16)
        return fig
    coords = adata.obsm["spatial"]
    vals = _to_1d_vals(adata, gene)
    sca = plt.scatter(coords[:, 0], coords[:, 1], s=20, c=vals, cmap="viridis")
    plt.xticks([]);
    plt.yticks([])
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


def _plot_spatial_expr_mrna(adata, gene: Optional[str], library_id: Optional[str],
                            img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    if gene is None or str(gene) not in _varnames(adata):
        ax.axis("off")
        ax.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16, transform=ax.transAxes)
        return fig
    if library_id and img_key:
        try:
            sc.pl.spatial(
                adata, color=str(gene), library_id=library_id, img_key=img_key,
                show=False, ax=ax, size=1.7, cmap="viridis"
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
                adata, color=str(gene), library_id=library_id, img_key=img_key,
                show=False, ax=ax, size=1.7, cmap="plasma"
            )
            ax.set_xlim(ax.get_xlim())
            ax.set_ylim(ax.get_ylim())
            ax.set_aspect('equal', adjustable='box')
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            return fig
        except Exception:
            plt.close(fig)
    return _plot_scatter_expr(adata, gene)



has_img_out, lib_id_out, img_key_out, spatial_meta_out = _probe_spatial_meta(adata_out)
has_img_in, lib_id_in, img_key_in, spatial_meta_in = _probe_spatial_meta(adata_in)

if len(protein_names) == 0:
    st.info("No proteins found in adata_out.var_names")
    gene = None
else:
    def format_func(gene_name):
        return gene_name


    gene = st.selectbox("Select a protein to visualize:", options=protein_names, index=0, format_func=format_func)


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

# Columns for plots
g1, g2, g3 = st.columns(3)

with g1:
    st.caption("Tissue image")
    if has_img_out:
        fig1 = _plot_tissue_only(adata_out, lib_id_out, img_key_out)
    elif has_img_in:
        fig1 = _plot_tissue_only(adata_in, lib_id_in, img_key_in)
    else:
        fig1 = _plot_tissue_only(adata_out, None, None)
    st.pyplot(fig1, use_container_width=True)
    plt.close(fig1)

with g2:
    st.caption("Spatial Protein (Imputed)")
    fig2 = _plot_spatial_expr(adata_out, gene, lib_id_out if has_img_out else None,
                              img_key_out if has_img_out else None)
    st.pyplot(fig2, use_container_width=True)
    plt.close(fig2)

with g3:
    st.caption("Spatial mRNA (Normalized)")
    if gene is None or str(gene) not in _varnames(adata_in):
        fig3 = _plot_image_placeholder(IMAGE_NA_PATH)
    else:
        fig3 = _plot_spatial_expr_mrna(adata_in, gene, lib_id_in if has_img_in else None,
                                       img_key_in if has_img_in else None)
    st.pyplot(fig3, use_container_width=True)
    plt.close(fig3)