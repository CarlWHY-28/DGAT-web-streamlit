import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import tempfile
import os
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
from scipy import sparse as sp

from dgat_utils.downstream import compute_bivariate_moran_single

FIGSIZE = (4.8, 4.8)
IMAGE_NA_PATH = "./logo/no_available_icon.png"

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Colocalization Analysis</h2>",
            unsafe_allow_html=True)
st.write("")

adata_in = st.session_state.get("adata_in", None)
adata_out = st.session_state.get("adata_out", None)
protein_names = st.session_state.get("protein_names", [])

if (adata_in is None) or (adata_out is None):
    st.warning("No uploaded result found. Please go to **Upload Data** and run again.")
    st.stop()

if "bivariate_moran_df" not in st.session_state:
    with st.spinner("Calculating Bivariate Moran's I for all protein pairs (this may take a moment)..."):
        try:

            bivariate_df = compute_bivariate_moran_single(
                adata_out,
                tissue_name="Current Sample",
                output_dir=None,
                coord_type="grid"
            )
            st.session_state["bivariate_moran_df"] = bivariate_df
        except Exception as e:
            st.error(f"Error computing Bivariate Moran's I: {e}")
            st.session_state["bivariate_moran_df"] = pd.DataFrame()

bivariate_df = st.session_state["bivariate_moran_df"]

st.subheader("1. Colocalization Statistics (Bivariate Moran's I)")

search_col, download_col = st.columns([3, 1])
with search_col:
    search_term = st.text_input("🔍 Search for a protein pair (e.g. 'CD3', 'CD19')", "")

if not bivariate_df.empty:
    display_df = bivariate_df.copy()

    if search_term:
        mask = (
                display_df['Marker_A'].astype(str).str.contains(search_term, case=False) |
                display_df['Marker_B'].astype(str).str.contains(search_term, case=False)
        )
        display_df = display_df[mask]

    st.dataframe(
        display_df,
        use_container_width=True,
        height=300,
        column_config={
            "Bivariate_Moran_I": st.column_config.NumberColumn(
                "Moran's I",
                format="%.4f"
            )
        }
    )

    st.caption(f"Showing {len(display_df)} pairs.")

    with download_col:
        st.write("")  # Spacer
        st.write("")
        csv = display_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="⬇️ Download CSV",
            data=csv,
            file_name='bivariate_moran_colocalization.csv',
            mime='text/csv',
        )
else:
    st.warning("Calculation returned empty result.")

st.divider()


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
    sca = plt.scatter(coords[:, 0], coords[:, 1], s=20, c=vals, cmap="plasma")
    plt.xticks([]);
    plt.yticks([])
    plt.colorbar(sca, shrink=0.75).set_label(str(gene))
    return fig


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



st.subheader('2. Colocalization Visualization')

has_img_out, lib_id_out, img_key_out, spatial_meta_out = _probe_spatial_meta(adata_out)
has_img_in, lib_id_in, img_key_in, spatial_meta_in = _probe_spatial_meta(adata_in)

if len(protein_names) < 2:
    st.error("Not enough proteins to perform colocalization analysis (Need at least 2).")
else:
    sel_c1, sel_c2 = st.columns(2)
    with sel_c1:
        prot_a = st.selectbox("Select Protein A", options=protein_names, index=0, key="prot_a")
    with sel_c2:
        default_idx_b = 1 if len(protein_names) > 1 else 0
        prot_b = st.selectbox("Select Protein B", options=protein_names, index=default_idx_b, key="prot_b")

    if not bivariate_df.empty:
        pair_row = bivariate_df[
            ((bivariate_df['Marker_A'] == prot_a) & (bivariate_df['Marker_B'] == prot_b)) |
            ((bivariate_df['Marker_A'] == prot_b) & (bivariate_df['Marker_B'] == prot_a))
            ]
        if not pair_row.empty:
            val = pair_row.iloc[0]['Bivariate_Moran_I']
            st.info(f"🔗 **Colocalization Index (Bivariate Moran's I)** between {prot_a} and {prot_b}: **{val:.4f}**")
        elif prot_a == prot_b:
            st.info("Selected the same protein. Bivariate Moran's I is essentially Univariate Moran's I.")
        else:
            st.info("Pair not found in computed statistics.")

    with st.expander("ℹ️ About Bivariate Moran's I & Preprocessing"):
        st.markdown("""
        **What is Bivariate Moran's I?**
        It measures the spatial correlation between two different variables (proteins).
        * **Positive value (approx > 0.1):** High expression of Protein A is surrounded by high expression of Protein B (Colocalization).
        * **Negative value:** High expression of Protein A is surrounded by low expression of Protein B (Mutually exclusive).
        * **Near Zero:** Random spatial relationship.

        **Data Info:**
        The visualization uses the **imputed protein expression** values. 
        """)

    g1, g2, g3 = st.columns(3)

    # 1. Tissue Image
    with g1:
        st.caption("Tissue Image")
        if has_img_out:
            fig1 = _plot_tissue_only(adata_out, lib_id_out, img_key_out)
        elif has_img_in:
            fig1 = _plot_tissue_only(adata_in, lib_id_in, img_key_in)
        else:
            fig1 = _plot_tissue_only(adata_out, None, None)
        st.pyplot(fig1, use_container_width=True)
        plt.close(fig1)

    # 2. Protein A
    with g2:
        st.caption(f"Protein A: {prot_a}")
        fig2 = _plot_spatial_expr(adata_out, prot_a, lib_id_out if has_img_out else None,
                                  img_key_out if has_img_out else None)
        st.pyplot(fig2, use_container_width=True)
        plt.close(fig2)

    # 3. Protein B
    with g3:
        st.caption(f"Protein B: {prot_b}")
        fig3 = _plot_spatial_expr(adata_out, prot_b, lib_id_out if has_img_out else None,
                                  img_key_out if has_img_out else None)
        st.pyplot(fig3, use_container_width=True)
        plt.close(fig3)