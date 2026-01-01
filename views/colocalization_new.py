import streamlit as st
import pandas as pd
import io
import os
import plotly.graph_objects as go
from streamlit_agraph import agraph, Node, Edge, Config
from task_manager import get_s3_client, get_image_url

# --- 1. 样式与配置 ---
st.markdown("""
    <style>
    /* 垂直居中所有列内容 */
    [data-testid="stHorizontalBlock"] {
        align-items: center;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown("<h2 style='text-align: center;'>Spatial Colocalization Analysis</h2>", unsafe_allow_html=True)

# 获取任务 ID
feature_code = st.session_state.get("current_feature_code") or st.query_params.get("task")
if not feature_code:
    st.warning("No Task ID found. Please run the analysis first.")
    st.stop()

s3 = get_s3_client()
bucket = os.getenv("BUCKET_NAME")
plot_prefix = f"task_{feature_code}/spatial_plots"


# --- 2. 数据加载 (从 S3 读取计算好的 Bivariate CSV) ---
@st.cache_data(ttl=600)
def load_bivariate_data(f_code):
    try:
        # 假设 Worker 计算完后存为该路径
        obj = s3.get_object(Bucket=bucket, Key=f"task_{f_code}/bivariate_moran_colocalization.csv")
        return pd.read_csv(io.BytesIO(obj['Body'].read()))
    except:
        return None


bivariate_df = load_bivariate_data(feature_code)

if bivariate_df is None:
    st.info("Colocalization data is being processed. Please refresh in a moment.")
    st.stop()

# --- 3. 统计表格展示 ---
st.subheader("1. Colocalization Statistics")
search_term = st.text_input("🔍 Search for a protein pair (e.g. 'CD3')", "")

display_df = bivariate_df.copy()
if search_term:
    mask = (display_df['Marker_A'].str.contains(search_term, case=False) |
            display_df['Marker_B'].str.contains(search_term, case=False))
    display_df = display_df[mask]

st.dataframe(display_df, height=250)

# --- 4. 空间可视化 (Gallery 模式) ---
st.divider()
st.subheader("2. Colocalization Visualization")

protein_names = sorted(list(set(bivariate_df['Marker_A']).union(set(bivariate_df['Marker_B']))))
sel_c1, sel_c2 = st.columns(2)
with sel_c1:
    prot_a = st.selectbox("Select Protein A", options=protein_names, index=0)
with sel_c2:
    prot_b = st.selectbox("Select Protein B", options=protein_names, index=min(1, len(protein_names) - 1))

# 布局：Tissue (90%) | Protein A | Protein B
g1, g2, g3 = st.columns([1, 1, 1])

with g1:
    st.caption("Tissue Image")
    # 嵌套列实现 90% 宽度并居中
    sub_l, sub_mid, sub_r = st.columns([0.05, 0.9, 0.05])
    with sub_mid:
        # 直接调用 Worker 预存的 Tissue 图
        st.image(get_image_url(s3, bucket, f"{plot_prefix}/tissue.png"))

with g2:
    st.caption(f"Protein A: {prot_a}")
    # 直接调用 Worker 预存的单蛋白图 (假设文件名格式为 protein_NAME.png)
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/protein_{prot_a}.png"))

with g3:
    st.caption(f"Protein B: {prot_b}")
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/protein_{prot_b}.png"))

# --- 5. 交互式分析 (热图 & 网络) ---
st.divider()
st.subheader("3. Interaction Analysis")

tab1, tab2 = st.tabs(["Correlation Heatmap", "Interaction Network"])

with tab1:
    # 热图数据准备
    matrix_df = bivariate_df.pivot(index='Marker_A', columns='Marker_B', values='Bivariate_Moran_I')
    fig_heatmap = go.Figure(data=go.Heatmap(
        z=matrix_df.values, x=matrix_df.columns, y=matrix_df.index,
        colorscale='RdBu_r', zmid=0
    ))
    fig_heatmap.update_layout(width=700, height=600, margin=dict(t=20, b=20, l=50, r=50))
    st.plotly_chart(fig_heatmap)

with tab2:
    threshold = st.slider("Edge Threshold (Abs Moran's I)", 0.0, 1.0, 0.5, 0.05)

    # 后面这段 agraph 逻辑保持不变，因为它本身就是基于 dataframe 数值实时生成的，不涉及复杂图像处理
    edges_df = bivariate_df[bivariate_df['Bivariate_Moran_I'].abs() >= threshold]
    nodes = [Node(id=m, label=m, size=15) for m in set(edges_df['Marker_A']).union(set(edges_df['Marker_B']))]
    edges = [Edge(source=r['Marker_A'], target=r['Marker_B'],
                  color="#33FF57" if r['Bivariate_Moran_I'] > 0 else "#FF3357")
             for _, r in edges_df.iterrows()]

    config = Config(width=700, height=500, physics=True, interaction={'navigationButtons': True})
    agraph(nodes=nodes, edges=edges, config=config)