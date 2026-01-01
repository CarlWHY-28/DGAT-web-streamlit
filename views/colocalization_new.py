import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go
from streamlit_agraph import agraph, Node, Edge, Config
from scipy import sparse as sp
from typing import Optional, Tuple, List, Dict, Any

# 假设这个函数存在于你的工具包中
# 如果没有，请确保你有正确的导入路径，或保留你原有的 compute_moran_single 逻辑
from dgat_utils.downstream import compute_bivariate_moran_single

# ==========================================
# 0. 配置与辅助函数
# ==========================================
FIGSIZE = (4.8, 4.8)
IMAGE_NA_PATH = "./logo/no_available_icon.png"


def _varnames(adata) -> List[str]:
    return [str(x) for x in getattr(adata, "var_names", [])]


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
        ax.set_aspect('equal', adjustable='box')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        return fig
    except Exception:
        plt.close(fig)
    return plt.figure()


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
            ax.set_aspect('equal', adjustable='box')
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            return fig
        except Exception:
            plt.close(fig)

    # Fallback scatter
    coords = adata.obsm["spatial"]
    X = adata[:, str(gene)].X
    vals = np.ravel(X.A) if sp.issparse(X) else np.asarray(X).ravel()
    sca = plt.scatter(coords[:, 0], coords[:, 1], s=20, c=vals, cmap="plasma")
    plt.colorbar(sca, shrink=0.75).set_label(str(gene))
    plt.axis("off")
    return fig


# ==========================================
# 1. 页面初始化与计算逻辑
# ==========================================
st.markdown("<h2 style='text-align: center; color: black;'>Spatial Colocalization Analysis</h2>",
            unsafe_allow_html=True)
st.write("")

adata_in = st.session_state.get("adata_in", None)
adata_out = st.session_state.get("adata_out", None)
protein_names = st.session_state.get("protein_names", [])

if (adata_in is None) or (adata_out is None):
    st.warning("No uploaded result found. Please go to **Upload Data** and run again.")
    st.stop()

# 自动计算 Bivariate Moran's I
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

# ==========================================
# 2. URL 参数处理 (Deep Linking)
# ==========================================
query_params = st.query_params

# 定义默认索引
default_idx_a = 0
default_idx_b = 1 if len(protein_names) > 1 else 0

# 如果 URL 中有参数，覆盖默认值
if "prot_a" in query_params and query_params["prot_a"] in protein_names:
    try:
        default_idx_a = protein_names.index(query_params["prot_a"])
    except ValueError:
        pass

if "prot_b" in query_params and query_params["prot_b"] in protein_names:
    try:
        default_idx_b = protein_names.index(query_params["prot_b"])
    except ValueError:
        pass


# 回调函数：当下拉框改变时，更新 URL
def update_url_params():
    st.query_params["prot_a"] = st.session_state.sel_prot_a
    st.query_params["prot_b"] = st.session_state.sel_prot_b


# ==========================================
# 3. 统计表格展示
# ==========================================
st.subheader("1. Colocalization Statistics")

if not bivariate_df.empty:
    with st.expander("Show Statistics Table", expanded=False):
        st.dataframe(
            bivariate_df,
            use_container_width=True,
            height=300,
            column_config={"Bivariate_Moran_I": st.column_config.NumberColumn("Moran's I", format="%.4f")}
        )
        csv = bivariate_df.to_csv(index=False).encode('utf-8')
        st.download_button("⬇️ Download CSV", csv, "bivariate_moran.csv", "text/csv")

st.divider()

# ==========================================
# 4. 空间可视化 (受 URL 控制)
# ==========================================
st.subheader('2. Colocalization Visualization')

if len(protein_names) < 2:
    st.error("Need at least 2 proteins.")
else:
    # 选择框：绑定到 session_state 并在变化时调用 update_url_params
    c_sel1, c_sel2 = st.columns(2)
    with c_sel1:
        prot_a = st.selectbox("Select Protein A", options=protein_names, index=default_idx_a,
                              key="sel_prot_a", on_change=update_url_params)
    with c_sel2:
        prot_b = st.selectbox("Select Protein B", options=protein_names, index=default_idx_b,
                              key="sel_prot_b", on_change=update_url_params)

    # 显示计算结果
    if not bivariate_df.empty:
        pair_row = bivariate_df[
            ((bivariate_df['Marker_A'] == prot_a) & (bivariate_df['Marker_B'] == prot_b)) |
            ((bivariate_df['Marker_A'] == prot_b) & (bivariate_df['Marker_B'] == prot_a))
            ]
        if not pair_row.empty:
            val = pair_row.iloc[0]['Bivariate_Moran_I']
            st.info(f"🔗 **Colocalization Index** ({prot_a} ↔ {prot_b}): **{val:.4f}**")

    # 绘图区域
    has_img_out, lib_id_out, img_key_out, _ = _probe_spatial_meta(adata_out)

    g1, g2, g3 = st.columns(3)

    with g1:
        st.caption("Tissue Image")
        fig1 = _plot_tissue_only(adata_out, lib_id_out, img_key_out)
        st.pyplot(fig1, use_container_width=True)
        plt.close(fig1)

    with g2:
        st.caption(f"Protein A: {prot_a}")
        fig2 = _plot_spatial_expr(adata_out, prot_a, lib_id_out, img_key_out)
        st.pyplot(fig2, use_container_width=True)
        plt.close(fig2)

    with g3:
        st.caption(f"Protein B: {prot_b}")
        fig3 = _plot_spatial_expr(adata_out, prot_b, lib_id_out, img_key_out)
        st.pyplot(fig3, use_container_width=True)
        plt.close(fig3)

st.divider()

# ==========================================
# 5. 高级交互 (Heatmap & Network)
# ==========================================
st.subheader("3. Advanced Interaction: Heatmap & Network")

# 准备矩阵数据
matrix_df = bivariate_df.pivot(index='Marker_A', columns='Marker_B', values='Bivariate_Moran_I')
all_markers = sorted(list(set(bivariate_df['Marker_A']).union(set(bivariate_df['Marker_B']))))
matrix_df = matrix_df.reindex(index=all_markers, columns=all_markers)

# --- A. Heatmap (热图) ---
st.markdown("#### A. Correlation Heatmap")
st.caption("Click on a cell to update the visualizations above (syncs to URL).")

fig_heatmap = go.Figure(data=go.Heatmap(
    z=matrix_df.values, x=matrix_df.columns, y=matrix_df.index,
    colorscale='RdBu_r', zmid=0, xgap=1, ygap=1,
    showscale=True, hoverongaps=False
))
fig_heatmap.update_layout(
    width=600, height=500, xaxis_showgrid=False, yaxis_showgrid=False,
    xaxis_tickangle=-45, margin=dict(t=20, b=20, l=50, r=50),
    plot_bgcolor='rgba(0,0,0,0)'
)

# 交互事件：点击热图
select_event = st.plotly_chart(fig_heatmap, use_container_width=False, on_select="rerun", selection_mode="points")

if select_event and len(select_event.selection["points"]) > 0:
    point = select_event.selection["points"][0]
    clicked_x = point["x"]  # Marker B
    clicked_y = point["y"]  # Marker A

    if clicked_x in matrix_df.columns and clicked_y in matrix_df.index:
        # 如果点击了热图，强制更新 URL 并重跑
        if (clicked_y != st.session_state.sel_prot_a) or (clicked_x != st.session_state.sel_prot_b):
            st.query_params["prot_a"] = clicked_y
            st.query_params["prot_b"] = clicked_x
            st.rerun()

st.divider()

# --- B. Network (网络图 - 带Marker颜色) ---
st.markdown("#### B. Interaction Network")

# 1. 定义 Biomarker 颜色信息 (恢复这部分)
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

# 2. 控制区 (垂直布局)
threshold = st.slider("Edge Threshold (Abs Moran's I)", 0.0, 1.0, 0.1, 0.05)

# 3. 渲染 HTML 图例 (带颜色的小圆点)
legend_html = "<div style='margin-bottom: 10px; display: flex; flex-wrap: wrap; gap: 10px;'>"
for type_name, color in TYPE_COLORS.items():
    item_html = f"<div style='display: flex; align-items: center;'><span style='width: 12px; height: 12px; background-color: {color}; border-radius: 50%; display: inline-block; margin-right: 5px;'></span><span style='font-size: 13px; color: #333;'>{type_name}</span></div>"
    legend_html += item_html
legend_html += "</div>"
st.markdown(legend_html, unsafe_allow_html=True)
st.caption("Drag nodes to rearrange. Thicker edges = stronger correlation.")

# 4. 构建网络数据
edges_df = bivariate_df[bivariate_df['Bivariate_Moran_I'].abs() >= threshold].reset_index(drop=True)
unique_nodes = set(edges_df['Marker_A']).union(set(edges_df['Marker_B']))

nodes = []
for marker in unique_nodes:
    # --- 关键修改：应用颜色 ---
    m_type = MARKER_INFO.get(marker, "Other")
    m_color = TYPE_COLORS.get(m_type, TYPE_COLORS["Other"])

    nodes.append(Node(
        id=marker,
        label=marker,
        size=12,
        font={'size': 14, 'color': 'black'},
        color=m_color,  # 使用分类颜色
        title=f"Type: {m_type}"
    ))

edges = []
for _, row in edges_df.iterrows():
    val = row['Bivariate_Moran_I']
    edge_id = f"{row['Marker_A']}|{row['Marker_B']}|{val}"
    edge_color = "#33FF57" if val > 0 else "#FF3357"
    # 动态边粗细
    calc_width = 1.0 + (abs(val) * 6.0)

    edges.append(Edge(
        source=row['Marker_A'],
        target=row['Marker_B'],
        id=edge_id,
        label="",
        color=edge_color,
        width=calc_width,
        title=f"Moran's I: {val:.4f}"
    ))

# 5. 配置与渲染 (带缩放按钮)
config = Config(
    width=700,
    height=500,
    directed=False,
    physics=True,
    hierarchical=False,
    nodeHighlightBehavior=True,
    highlightColor="#F7A7A6",
    collapsible=False,
    interaction={'navigationButtons': True, 'dragView': True, 'zoomView': True}
)

selected_id = agraph(nodes=nodes, edges=edges, config=config)

# 6. 网络图交互反馈
if selected_id:
    if "|" in selected_id:
        m_a, m_b, val_str = selected_id.split("|")
        st.info(f"**Selected Edge:** {m_a} ↔ {m_b} (I={float(val_str):.4f})")
    else:
        st.info(f"**Selected Node:** {selected_id} ({MARKER_INFO.get(selected_id, 'Other')})")