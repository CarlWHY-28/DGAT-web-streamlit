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

# --- 诊断代码段 ---
with st.expander("🔍 Debug: Check Bucket Files"):
    try:
        s3 = get_s3_client()
        bucket = os.getenv("BUCKET_NAME")
        prefix = f"task_{feature_code}/spatial_plots/"
        response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)

        if 'Contents' in response:
            st.write("✅ Files found in S3:")
            for obj in response['Contents']:
                st.write(f"- {obj['Key']}")
        else:
            st.error(f"❌ No files found in prefix: {prefix}. Did the Worker finish drawing?")
    except Exception as e:
        st.error(f"Error connecting to S3: {e}")
# --- 诊断结束 ---


if not feature_code:
    st.warning("Please query a feature code first.")
    st.stop()

s3 = get_s3_client()
bucket = os.getenv("BUCKET_NAME")
plot_prefix = f"task_{feature_code}/spatial_plots"


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



# --- 1. 蛋白质选择与展示 ---
protein_names = st.session_state.get("protein_names", [])
selected_p = st.selectbox("Select Protein", protein_names)



# 1. 注入 CSS 样式
st.markdown("""
    <style>
    /* 定位到 columns 的父容器，强制其子元素垂直居中 */
    [data-testid="stHorizontalBlock"] {
        align-items: center;
    }
    </style>
    """, unsafe_allow_html=True)

# 2. 正常的列布局
col1, col2, col3 = st.columns(3)

with col1:
    # 标题建议放在 sub_mid 内部或图片上方
    sub_l, sub_mid, sub_r = st.columns([0.02, 0.96, 0.02])
    with sub_mid:
        st.image(
            get_image_url(s3, bucket, f"{plot_prefix}/tissue.png"),
            channels="RGB"
        )




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