import streamlit as st
import os
from dgat_utils.task_manager import get_s3_client

# 核心：生成 Presigned URL 的辅助函数
def get_image_url(s3_client, bucket, key):
    return s3_client.generate_presigned_url(
        'get_object',
        Params={'Bucket': bucket, 'Key': key},
        ExpiresIn=1800 # 30分钟有效期
    )

st.title("View Result Gallery")

# 从 Session State 获取当前查看的特征码
feature_code = st.session_state.get("current_feature_code")
if not feature_code:
    st.warning("Please query a feature code first.")
    st.stop()

s3 = get_s3_client()
bucket = os.getenv("BUCKET_NAME")
plot_prefix = f"task_{feature_code}/plots"

# --- 1. 蛋白质选择与展示 ---
protein_names = st.session_state.get("protein_names", [])
selected_p = st.selectbox("Select Protein", protein_names)

col1, col2, col3 = st.columns(3)
with col1:
    st.caption("Tissue Image")
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/tissue.png"))

with col2:
    st.caption(f"Protein: {selected_p} (Imputed)")
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/protein_{selected_p}.png"))

with col3:
    st.caption(f"mRNA: {selected_p} (Original)")
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/mrna_{selected_p}.png"))

st.divider()

# --- 2. Leiden 聚类图片浏览器 ---
st.subheader("Spatial Leiden Clustering")
c1, c2 = st.columns(2)

with c1:
    res_p = st.select_slider("Protein Resolution", options=[0.3, 0.4, 0.5, 0.6], value=0.5)
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/leiden_prot_{res_p}.png"))

with c2:
    res_m = st.select_slider("mRNA Resolution", options=[0.8, 0.9, 1.0, 1.1], value=1.0)
    st.image(get_image_url(s3, bucket, f"{plot_prefix}/leiden_mrna_{res_m}.png"))