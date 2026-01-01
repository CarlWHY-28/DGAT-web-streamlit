import os

import streamlit as st
import pandas as pd
import io
import requests
from dgat_utils.task_manager import get_s3_client, get_image_url

# --- 1. 页面配置与 CSS 注入 ---
st.markdown("""
    <style>
    /* 强制 columns 内容垂直居中 */
    [data-testid="stHorizontalBlock"] {
        align-items: center;
    }
    .centered-text {
        text-align: center;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Correlation Analysis</h2>", unsafe_allow_html=True)

# --- 2. 获取任务上下文 ---
# 假设你在 URL 参数或 SessionState 中存了 feature_code
feature_code = st.session_state.get("current_feature_code") or st.query_params.get("task")

if not feature_code:
    st.warning("No task ID found. Please go to **Upload Data** or check your link.")
    st.stop()

s3 = get_s3_client()
bucket = os.getenv("BUCKET_NAME")
plot_prefix = f"task_{feature_code}/spatial_plots"


# --- 3. 加载 Moran's I 数据 (从 S3 读取 CSV) ---
# 注意：你的 Worker 需要在计算完后把 moran_df 存为 task_{feature_code}/moran_statistics.csv
@st.cache_data(ttl=3600)
def load_moran_data(f_code):
    try:
        obj = s3.get_object(Bucket=bucket, Key=f"task_{f_code}/moran_statistics.csv")
        df = pd.read_csv(io.BytesIO(obj['Body'].read()))
        # 假设统计数据存放在 S3 的另一个小 json 里，或者直接从 df 计算
        return df
    except Exception:
        return None


moran_df = load_moran_data(feature_code)

# --- 4. 渲染统计指标 ---
st.subheader("1. Global Spatial Autocorrelation Analysis")

if moran_df is not None:
    # 简单的统计展示
    m1, m2 = st.columns(2)
    with m1:
        # 这里假设你 CSV 里存了这些汇总信息，或者从 session_state 获取
        st.metric("Analyzed Markers", value=len(moran_df))
    with m2:
        st.write("Moran's I reflects the spatial consistency between Protein and mRNA.")

    st.dataframe(moran_df, height=250)

    # 下载按钮
    csv_data = moran_df.to_csv(index=False).encode('utf-8')
    st.download_button("⬇️ Download Moran's I Table", data=csv_data, file_name='moran_stats.csv')
else:
    st.info("Moran's I table is still being processed or not available.")

st.divider()

# --- 5. 空间表达可视化 (Gallery 模式) ---
st.subheader("2. Spatial Expression Visualization")

# 选择蛋白质 (从 CSV 的 marker 列获取)
if moran_df is not None:
    protein_list = moran_df['marker'].unique().tolist()
    selected_p = st.selectbox("Select a protein to visualize:", options=protein_list)
else:
    # 备选方案：手动定义或从 session 获取
    selected_p = st.selectbox("Select a protein:", options=st.session_state.get("protein_names", []))

# 布局：Tissue(90%) | Protein | mRNA
col1, col2, col3 = st.columns([1, 1, 1])

with col1:
    st.caption("Tissue Image")
    # 嵌套子列实现 90% 宽度并居中
    sub_l, sub_mid, sub_r = st.columns([0.05, 0.9, 0.05])
    with sub_mid:
        tissue_url = get_image_url(s3, bucket, f"{plot_prefix}/tissue.png")
        st.image(tissue_url)

with col2:
    st.caption(f"Protein: {selected_p} (Imputed)")
    prot_url = get_image_url(s3, bucket, f"{plot_prefix}/protein_{selected_p}.png")
    st.image(prot_url)

with col3:
    st.caption(f"mRNA: {selected_p} (Original)")
    mrna_url = get_image_url(s3, bucket, f"{plot_prefix}/mrna_{selected_p}.png")
    st.image(mrna_url)

# --- 6. 数据下载区域 ---
st.divider()
st.subheader("3. Data Download")
d_col1, d_col2 = st.columns(2)

with d_col1:
    # 生成预签名的 H5AD 下载链接
    h5ad_url = get_image_url(s3, bucket, f"task_{feature_code}/output.h5ad", expires_in=3600)
    st.link_button("⬇️ Download Predicted Data (.h5ad)", h5ad_url)

with d_col2:
    st.info("The .h5ad file contains imputed protein layers and spatial coordinates.")