import os

import streamlit as st
import pandas as pd
import io
import requests
from dgat_utils.task_manager import get_s3_client, get_image_url

st.markdown("""
    <style>
    [data-testid="stHorizontalBlock"] {
        align-items: center;
    }
    .centered-text {
        text-align: center;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown("<h2 style='text-align: center; color: black;'>Spatial Correlation Analysis</h2>", unsafe_allow_html=True)


feature_code = st.session_state.get("current_feature_code") or st.query_params.get("task")

if not feature_code:
    st.warning("No task ID found. Please go to **Upload Data** or check your link.")
    st.stop()

s3 = get_s3_client()
bucket = os.getenv("BUCKET_NAME")
plot_prefix = f"task_{feature_code}/spatial_plots"


@st.cache_data(ttl=3600)
def load_moran_data(f_code):
    try:
        obj = s3.get_object(Bucket=bucket, Key=f"task_{f_code}/moran_statistics.csv")
        df = pd.read_csv(io.BytesIO(obj['Body'].read()))
        return df
    except Exception:
        return None


moran_df = load_moran_data(feature_code)

st.subheader("1. Global Spatial Autocorrelation Analysis")

if moran_df is not None:
    m1, m2 = st.columns(2)
    with m1:
        st.write(f"Analyzed Markers: {len(moran_df)}")
    with m2:
        st.write("Moran's I reflects the spatial consistency")

    st.dataframe(moran_df, use_container_width=True, height=300)
    csv_data = moran_df.to_csv(index=False).encode('utf-8')
    st.download_button("⬇️ Download Moran's I Table", data=csv_data, file_name='moran_stats.csv')
else:
    st.info("Moran's I table is still being processed or not available.")

st.divider()

st.subheader("2. Spatial Expression Visualization")

if moran_df is not None:
    protein_list = moran_df['marker'].unique().tolist()
    selected_p = st.selectbox("Select a protein to visualize:", options=protein_list)
else:
    selected_p = st.selectbox("Select a protein:", options=st.session_state.get("protein_names", []))

col1, col2, col3 = st.columns([1, 1, 1])

with col1:
    st.caption("Tissue Image")
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
