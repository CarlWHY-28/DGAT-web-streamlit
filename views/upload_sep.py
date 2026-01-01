import streamlit as st
import uuid
import os
from dgat_utils.task_manager import Session, ProteinTask, get_s3_client
import gc
st.header("Upload and Submit Your Task")
st.info("Upload your `.h5ad` file and we will notify you via email once computation is done.")

import requests

URL_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'

common_protein = requests.get(f"{URL_REPO}/common_protein_31.txt").text.strip().splitlines()
st.session_state["protein_names"] = common_protein

if "uploading" not in st.session_state:
    st.session_state.uploading = False


with st.form("upload_form"):
    email = st.text_input("Email Address (for notification)")
    uploaded_file = st.file_uploader("Choose an h5ad file", type="h5ad")
    submit = st.form_submit_button(
        "Submit Task",
        type="primary",
        disabled=st.session_state.uploading  # ⛔上传时禁用
    )



if submit:
    if uploaded_file and email:
        st.session_state.uploading = True
        st.session_state["current_feature_code"] = None
        feature_code = str(uuid.uuid4())[:8].upper()
        bucket_name = os.getenv("BUCKET_NAME")
        input_key = f"task_{feature_code}/input.h5ad"

        try:
            with st.spinner("Uploading to secure storage..."):
                s3 = get_s3_client()
                s3.upload_fileobj(uploaded_file, bucket_name, input_key)

                session = Session()
                new_task = ProteinTask(
                    feature_code=feature_code,
                    email=email,
                    input_path=input_key,
                    status='pending'
                )
                session.add(new_task)
                session.commit()
                session.close()



            st.success(f"✅ Submitted! Your Feature Code: **{feature_code}**")
            st.balloons()
        except Exception as e:
            st.error(f"Submission failed: {e}")
        del uploaded_file

        # 4. 强制进行垃圾回收
        gc.collect()
    else:
        st.warning("Please provide both file and email.")
    st.session_state.uploading = False


st.write("")
st.divider()
tab1, tab2 = st.tabs(["Frequently Asked Questions (FAQ)", "Data Requirements"])

with tab1:
    st.subheader("Frequently Asked Questions (FAQ)")
    st.markdown("""
        **Q1: What is this upload page for?**  
        This page is used to upload your spatial transcriptomics data (in `h5ad` format) for protein imputation.

        **Q2: What information do I need to provide?**  
        Only an `.h5ad` file is required. No personal information is needed. We will not store your data beyond the session so feel free to upload sensitive data.

        **Q3: How long does the upload and processing take?**  
        It depends on your file size and Internet speed. Please keep this tab open during processing. Normally, it takes 5~10 minutes for files under 300MB. For a faster imputation, consider using our DGAT on local machines, view our GitHub Repo [here](https://github.com/osmanbeyoglulab/DGAT).
    """)

with tab2:
    st.subheader("Data Requirements")
    # st.markdown("""
    #     To ensure the analysis tools can process your data correctly, please make sure your `.h5ad` file meets the following requirements:

    #     **1. File Format:** `h5ad` (Anndata)

    #     **2. Required Data Slots:**
    #     * `adata.X`: Raw count matrix (prefer non-normalized counts; preprocessing will be handled by our pipeline).
    #     * `adata.obs`: Cell/Spot metadata.
    #     * `adata.var`: Gene metadata (e.g., gene names).
    #     * `adata.obsm['spatial']`: **(CRITICAL)** Spatial coordinates (N×2 array for x/y). **Spatial analysis is impossible without this.**

    #     **3. File Size:**
    #     Keep the file size under ~200MB for smooth upload performance.
    # """)

    st.markdown("""


    To ensure successful processing and visualization, please prepare your dataset in the following format and structure **before uploading** to DGATviz.

    ---

    ##### **1. File Format**

    - The input file must be in **`.h5ad` (AnnData)** format.

    ---

    ##### **2. Required Data Components**

    - **`adata.X`** — Raw count matrix  
      *(Preferably non-normalized; preprocessing and normalization are handled automatically by the DGAT pipeline.)*

    - **`adata.obs`** — Cell or spot metadata  
      *(e.g., barcodes, sample identifiers, or annotations.)*

    - **`adata.var`** — Gene names  


    - **`adata.obsm['spatial']`** — Spatial coordinates as an **N×2 array** (x/y positions).  

      ⚠️ **This field is essential.** Spatial analysis cannot be performed without valid coordinates.

    ---

    ##### **3. File Size**

    - **Recommended maximum file size:** ≤ **300 MB**  
    """)

