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

# ... (Previous code) ...
# Paste this after your existing code

import plotly.graph_objects as go
from streamlit_agraph import agraph, Node, Edge, Config

st.subheader("3. Advanced Interaction: Heatmap & Network")

# 准备数据：将长表格转回矩阵形式用于热图
matrix_df = bivariate_df.pivot(index='Marker_A', columns='Marker_B', values='Bivariate_Moran_I')

# 确保矩阵行列顺序一致，为了美观
all_markers = sorted(list(set(bivariate_df['Marker_A']).union(set(bivariate_df['Marker_B']))))
matrix_df = matrix_df.reindex(index=all_markers, columns=all_markers)

# ==========================================
# 1. 交互式三角热图 (Plotly) - 单独一行
# ==========================================
st.markdown("#### A. Correlation Heatmap")
st.caption("Click on a cell to view details below.")

# 创建 Plotly 热图
fig_heatmap = go.Figure(data=go.Heatmap(
    z=matrix_df.values,
    x=matrix_df.columns,
    y=matrix_df.index,
    colorscale='RdBu_r',
    zmid=0,
    xgap=1, ygap=1,
    showscale=True,
    hoverongaps=False
))

fig_heatmap.update_layout(
    width=600, height=500,  # 适当增加尺寸使其更清晰
    xaxis_showgrid=False, yaxis_showgrid=False,
    xaxis_tickangle=-45,
    margin=dict(t=20, b=20, l=50, r=50),  # 增加边距防止标签被截断
    plot_bgcolor='rgba(0,0,0,0)'
)

# 使用 Streamlit 的 selection 事件
select_event = st.plotly_chart(
    fig_heatmap,
    use_container_width=False,  # 不强制占满，使用设定宽度
    on_select="rerun",
    selection_mode="points"
)

# 处理点击事件：输出信息
if select_event and len(select_event.selection["points"]) > 0:
    point = select_event.selection["points"][0]
    clicked_x = point["x"]
    clicked_y = point["y"]

    if clicked_x in matrix_df.columns and clicked_y in matrix_df.index:
        clicked_val = matrix_df.loc[clicked_y, clicked_x]
        if not pd.isna(clicked_val):
            st.info(f"**Selected Heatmap Pair:**\n\n"
                    f"🧬 {clicked_y} ↔ {clicked_x}\n\n"
                    f"🔗 Moran's I: **{clicked_val:.4f}**")
        else:
            st.write("No data for this pair.")

st.divider()

# ==========================================
# 2. 可拖拽网络图 (Streamlit-Agraph) - 带标注颜色 (修复图例显示)
# ==========================================
st.markdown("#### B. Interaction Network")

# --- 0. 定义标注数据 (Biomarker Annotations) ---
MARKER_INFO = {
    "ACTA2": "Stromal / endothelial", "BCL2": "Basal epithelial", "CCR7": "T cell",
    "CD14": "Myeloid / APC", "CD163": "Myeloid / APC", "CD19": "B cell / TLS",
    "CD27": "T cell", "CD274": "Immune checkpoint", "CD3E": "T cell",
    "CD4": "T cell", "CD40": "B cell / TLS", "CD68": "Myeloid / APC",
    "CD8A": "T cell", "CEACAM8": "Myeloid / APC", "CR2": "B cell / TLS",
    "CXCR5": "B cell / TLS", "EPCAM": "Basal epithelial", "FCGR3A": "Myeloid / APC",
    "ITGAM": "Myeloid / APC", "ITGAX": "Myeloid / APC", "KRT5": "Basal epithelial",
    "MS4A1": "B cell / TLS", "PAX5": "B cell / TLS", "PCNA": "Proliferation",
    "PDCD1": "T cell", "PECAM1": "Stromal / endothelial", "SDC1": "Basal epithelial",
    "VIM": "Stromal / endothelial"
}

TYPE_COLORS = {
    "B cell / TLS": "#33A02C",
    "Basal epithelial": "#6A3D9A",
    "Immune checkpoint": "#D73027",
    "Myeloid / APC": "#E66101",
    "Proliferation": "#FFD92F",
    "Stromal / endothelial": "#B15928",
    "T cell": "#1F78B4",
    "Other": "#B0B0B0"
}

# --- 1. 控制区 (垂直布局) ---
threshold = st.slider("Edge Threshold (Abs Moran's I)",
                      min_value=0.0, max_value=1.0, value=0.5, step=0.05)

# --- 2. 显示图例 (Legend) - 修复版 ---
# 修复说明：将 HTML 拼接改为单行字符串，避免 Python 代码缩进被 st.markdown 误识别为代码块。
legend_html = "<div style='margin-bottom: 10px; display: flex; flex-wrap: wrap; gap: 10px;'>"
for type_name, color in TYPE_COLORS.items():
    # 使用 f-string 构建单行 HTML
    item_html = f"<div style='display: flex; align-items: center;'><span style='width: 12px; height: 12px; background-color: {color}; border-radius: 50%; display: inline-block; margin-right: 5px;'></span><span style='font-size: 13px; color: #333;'>{type_name}</span></div>"
    legend_html += item_html
legend_html += "</div>"

# 渲染 HTML 图例
st.markdown(legend_html, unsafe_allow_html=True)

st.caption("Drag nodes to rearrange. Node color represents cell type/function.")

# --- 3. 准备数据与构建网络 ---
edges_df = bivariate_df[bivariate_df['Bivariate_Moran_I'].abs() >= threshold].reset_index(drop=True)
unique_nodes = set(edges_df['Marker_A']).union(set(edges_df['Marker_B']))

nodes = []
for marker in unique_nodes:
    m_type = MARKER_INFO.get(marker, "Other")
    m_color = TYPE_COLORS.get(m_type, TYPE_COLORS["Other"])
    nodes.append(Node(
        id=marker,
        label=marker,
        size=12,
        font={'size': 14, 'color': 'black'},
        color=m_color,
        title=f"Type: {m_type}"
    ))

edges = []
for _, row in edges_df.iterrows():
    edge_id = f"{row['Marker_A']}|{row['Marker_B']}|{row['Bivariate_Moran_I']}"
    edge_color = "#33FF57" if row['Bivariate_Moran_I'] > 0 else "#FF3357"
    edges.append(Edge(
        source=row['Marker_A'],
        target=row['Marker_B'],
        id=edge_id,
        label="",
        color=edge_color,
        width=2
    ))

config = Config(
    width=700,
    height=500,
    directed=False,
    physics=True,
    hierarchical=False,
    nodeHighlightBehavior=True,
    highlightColor="#F7A7A6",
    collapsible=False,
    # --- 新增修改：开启导航按钮 ---
    # 这会在画布右下角添加 +, -, 以及复位箭头的按钮
    interaction={'navigationButtons': True, 'dragView': True, 'zoomView': True}
)

# --- 4. 渲染图表 ---
# key参数有时候有助于强制刷新，防止缓存导致的配置不生效，这里作为可选项
selected_id = agraph(nodes=nodes, edges=edges, config=config)

# 处理点击事件
if selected_id:
    if "|" in selected_id:
        parts = selected_id.split("|")
        if len(parts) == 3:
            m_a, m_b, val = parts
            st.info(f"**Selected Network Edge:**\n\n"
                    f"🧬 {m_a} ↔ {m_b}\n\n"
                    f"🔗 Moran's I: **{float(val):.4f}**")
    else:
        m_type = MARKER_INFO.get(selected_id, "Other")
        st.info(f"**Selected Protein Node:** {selected_id}\n\n"
                f"🏷️ Type: **{m_type}**")



