import streamlit as st
import pandas as pd
import io
import os
import plotly.graph_objects as go
from streamlit_agraph import agraph, Node, Edge, Config
from dgat_utils.task_manager import get_s3_client, get_image_url

# --- 1. 样式与布局设置 ---
st.markdown("""
    <style>
    /* 强制一行内的所有列垂直居中对齐 */
    [data-testid="stHorizontalBlock"] {
        align-items: center;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Colocalization Analysis</h2>",
            unsafe_allow_html=True)
st.write("")

# --- 2. 获取任务上下文 ---
feature_code = st.session_state.get("current_feature_code")
if not feature_code:
    st.warning("No Task ID found. Please run the analysis first.")
    st.stop()

s3 = get_s3_client()
bucket = os.getenv("BUCKET_NAME")
plot_prefix = f"task_{feature_code}/spatial_plots"


# --- 3. 加载计算结果 (从 S3 读取 CSV) ---
@st.cache_data(ttl=600)
def load_bivariate_data(f_code):
    """读取 Worker 计算好的双变量 Moran's I 表格"""
    try:
        obj = s3.get_object(Bucket=bucket, Key=f"task_{f_code}/bivariate_moran_colocalization.csv")
        return pd.read_csv(io.BytesIO(obj['Body'].read()))
    except Exception:
        return None


bivariate_df = load_bivariate_data(feature_code)

if bivariate_df is None:
    st.info("Colocalization statistics are being processed or not available. Please refresh later.")
    st.stop()

# --- 4. 统计表格展示 ---
st.subheader("1. Colocalization Statistics (Bivariate Moran's I)")

search_col, download_col = st.columns([3, 1])
with search_col:
    search_term = st.text_input("🔍 Search for a protein pair (e.g. 'CD3', 'CD19')", "")

display_df = bivariate_df.copy()

if search_term:
    mask = (
            display_df['Marker_A'].astype(str).str.contains(search_term, case=False) |
            display_df['Marker_B'].astype(str).str.contains(search_term, case=False)
    )
    display_df = display_df[mask]

st.dataframe(
    display_df,
    use_container_width= True,
    height=300,  # 移除 use_container_width=True
    column_config={
        "Bivariate_Moran_I": st.column_config.NumberColumn(
            "Moran's I",
            format="%.4f"
        )
    }
)
st.caption(f"Showing {len(display_df)} pairs.")

with download_col:
    st.write("")
    st.write("")
    csv = display_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="⬇️ Download CSV",
        data=csv,
        file_name='bivariate_moran_colocalization.csv',
        mime='text/csv',
    )

st.divider()

# --- 5. 空间可视化 (URL 图片模式) ---
st.subheader('2. Colocalization Visualization')

# 获取所有出现的蛋白名称
all_proteins = sorted(list(set(bivariate_df['Marker_A']).union(set(bivariate_df['Marker_B']))))

if len(all_proteins) < 2:
    st.error("Not enough proteins to perform colocalization analysis.")
else:
    sel_c1, sel_c2 = st.columns(2)
    with sel_c1:
        prot_a = st.selectbox("Select Protein A", options=all_proteins, index=0, key="prot_a")
    with sel_c2:
        default_idx_b = 1 if len(all_proteins) > 1 else 0
        prot_b = st.selectbox("Select Protein B", options=all_proteins, index=default_idx_b, key="prot_b")

    # 显示选定配对的数值
    pair_row = bivariate_df[
        ((bivariate_df['Marker_A'] == prot_a) & (bivariate_df['Marker_B'] == prot_b)) |
        ((bivariate_df['Marker_A'] == prot_b) & (bivariate_df['Marker_B'] == prot_a))
        ]
    if not pair_row.empty:
        val = pair_row.iloc[0]['Bivariate_Moran_I']
        st.info(f"🔗 **Colocalization Index** ({prot_a} ↔ {prot_b}): **{val:.4f}**")

    # --- 核心布局：Tissue(90%) | Protein A | Protein B ---
    g1, g2, g3 = st.columns([1, 1, 1])

    # 1. Tissue Image (特殊处理)
    with g1:
        st.caption("Tissue Image")
        # 嵌套列：左右留白 5%，中间 90%
        sub_l, sub_mid, sub_r = st.columns([0.05, 0.9, 0.05])
        with sub_mid:
            # 直接调用 URL，不使用 use_container_width
            tissue_url = get_image_url(s3, bucket, f"{plot_prefix}/tissue.png")
            st.image(tissue_url)

    # 2. Protein A Image (调用 URL)
    with g2:
        st.caption(f"Protein A: {prot_a}")
        pa_url = get_image_url(s3, bucket, f"{plot_prefix}/protein_{prot_a}.png")
        st.image(pa_url)

    # 3. Protein B Image (调用 URL)
    with g3:
        st.caption(f"Protein B: {prot_b}")
        pb_url = get_image_url(s3, bucket, f"{plot_prefix}/protein_{prot_b}.png")
        st.image(pb_url)

# --- 6. 高级交互分析 (基于 CSV 数据的轻量渲染) ---
st.divider()
st.subheader("3. Advanced Interaction: Heatmap & Network")

# 准备矩阵数据
matrix_df = bivariate_df.pivot(index='Marker_A', columns='Marker_B', values='Bivariate_Moran_I')
matrix_df = matrix_df.reindex(index=all_proteins, columns=all_proteins)

tab1, tab2 = st.tabs(["Correlation Heatmap", "Interaction Network"])

# === Tab 1: Heatmap ===
with tab1:
    st.caption("Click on a cell to view details.")
    fig_heatmap = go.Figure(data=go.Heatmap(
        z=matrix_df.values, x=matrix_df.columns, y=matrix_df.index,
        colorscale='RdBu_r', zmid=0, xgap=1, ygap=1
    ))
    fig_heatmap.update_layout(
        width=600, height=550,
        xaxis_showgrid=False, yaxis_showgrid=False,
        plot_bgcolor='rgba(0,0,0,0)',
        margin=dict(t=20, b=20, l=50, r=50)
    )

    select_event = st.plotly_chart(fig_heatmap, on_select="rerun", selection_mode="points")

    if select_event and len(select_event.selection["points"]) > 0:
        point = select_event.selection["points"][0]
        st.info(f"🧬 {point['y']} ↔ {point['x']} : **{matrix_df.loc[point['y'], point['x']]:.4f}**")

# === Tab 2: Network ===
with tab2:
    # 静态配置数据 (保留你原来的字典)
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
        "B cell / TLS": "#33A02C", "Basal epithelial": "#6A3D9A", "Immune checkpoint": "#D73027",
        "Myeloid / APC": "#E66101", "Proliferation": "#FFD92F", "Stromal / endothelial": "#B15928",
        "T cell": "#1F78B4", "Other": "#B0B0B0"
    }

    c_net1, c_net2 = st.columns([1, 3])

    with c_net1:
        threshold = st.slider("Edge Threshold", 0.0, 1.0, 0.5, 0.05)
        # 图例渲染
        legend_html = "<div style='display: flex; flex-direction: column; gap: 5px; font-size: 12px;'>"
        for t_name, color in TYPE_COLORS.items():
            legend_html += f"<div><span style='color:{color};'>●</span> {t_name}</div>"
        legend_html += "</div>"
        st.markdown(legend_html, unsafe_allow_html=True)

    with c_net2:
        edges_df = bivariate_df[bivariate_df['Bivariate_Moran_I'].abs() >= threshold]
        active_nodes = set(edges_df['Marker_A']).union(set(edges_df['Marker_B']))

        nodes = []
        for marker in active_nodes:
            m_type = MARKER_INFO.get(marker, "Other")
            nodes.append(Node(
                id=marker, label=marker, size=15,
                color=TYPE_COLORS.get(m_type, TYPE_COLORS["Other"]),
                title=f"Type: {m_type}"
            ))

        edges = []
        for _, row in edges_df.iterrows():
            edges.append(Edge(
                source=row['Marker_A'], target=row['Marker_B'],
                color="#33FF57" if row['Bivariate_Moran_I'] > 0 else "#FF3357"
            ))

        config = Config(width=600, height=500, physics=True, interaction={'navigationButtons': True},directed=False)
        agraph(nodes=nodes, edges=edges, config=config)