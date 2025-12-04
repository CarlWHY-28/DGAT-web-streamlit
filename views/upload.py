import streamlit as st
from datetime import datetime
from dgat_utils.predict_util import web_predict,test_Preprocessing
import anndata as ad
import traceback
import tempfile
import os
import scanpy as sc

st.session_state["adata_out"] = sc.read_h5ad('./data/CID44971.h5ad')

url_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'

st.markdown("<h2 style='text-align: center; color: black;'>Upload and Impute Your Data</h2>", unsafe_allow_html=True)
st.write("")
st.info("For testing purpose, we do not predict data. However, we still run the preprocessing steps, so please upload the data CID44971.")

if "results_ready" not in st.session_state:
    st.session_state["results_ready"] = False

with st.form(key="upload_form"):
    uploaded_file = st.file_uploader("Choose an h5ad file", type="h5ad", key="upload_file")
    submit_button = st.form_submit_button(label="Upload and Submit")

if submit_button:
    if uploaded_file is not None:
        now = datetime.now()
        timestamp_str = now.strftime("%Y%m%d-%H%M%S")
        st.session_state['current_time'] = timestamp_str

        with st.spinner(f"Uploading '{uploaded_file.name}' ..."):
            upload_success = True
            if upload_success:
                st.session_state['data_uploaded'] = True
                st.balloons()

                try:
                    uploaded_file.seek(0)
                    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_in:
                        tmp_in.write(uploaded_file.getbuffer())
                        tmp_in_path = tmp_in.name

                    st.warning("⚠️ Imputation is running. **Do NOT close or refresh this tab** until it finishes. Normally it takes 3~5 minutes for files under 300MB.", icon="⚠️")

                    with st.status("Preparing to run imputation…", state="running") as status:
                        status.update(label="Loading input file… (Do **not** close this page)", state="running")
                        adata_in = ad.read_h5ad(tmp_in_path)

                        status.update(label="Predicting proteins… (Do **not** close this page)", state="running")
                        ### Here we do not need adata_out
                        adata_in = test_Preprocessing(url_REPO, adata_in)
                        adata_out = st.session_state.get("adata_out", None)
                        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_out:
                            adata_out.write_h5ad(tmp_out.name)
                            final_out_path = tmp_out.name

                        status.update(label="Finalizing results…", state="complete")

                    with open(final_out_path, "rb") as f_out:
                        result_bytes = f_out.read()

                    download_name = f"{os.path.splitext(uploaded_file.name)[0]}_DGAT_pred_{timestamp_str}.h5ad"
                    st.success("✅ Prediction finished. You can download the result below:")
                    st.download_button(
                        label="⬇️ Download Predicted Data",
                        data=result_bytes,
                        file_name=download_name,
                        mime="application/octet-stream"
                    )

                    st.session_state["adata_in"] = adata_in

                    try:
                        st.session_state["protein_names"] = [str(x) for x in getattr(adata_out, "var_names", [])]
                    except Exception:
                        st.session_state["protein_names"] = []
                    st.session_state["has_upload"] = True
                    st.session_state["results_ready"] = True

                except Exception as e:
                    st.error("Prediction failed:")
                    st.exception(e)
                    st.code(traceback.format_exc(), language="python")
                finally:
                    try:
                        if 'tmp_in_path' in locals() and os.path.exists(tmp_in_path):
                            os.remove(tmp_in_path)
                        if 'final_out_path' in locals() and os.path.exists(final_out_path):
                            os.remove(final_out_path)
                    except Exception:
                        pass
            else:
                st.session_state['data_uploaded'] = False
    else:
        st.warning("Please select an `.h5ad` file before submitting.")

if st.session_state.get("results_ready", False):
    st.info("Your results are available under **Upload Your Data → View your data** in the left sidebar. You might download the imputed data there. Uploading a new file will overwrite the previous results. Please do not refresh this page to avoid losing the imputed data.", icon="ℹ️")
    col1, col2 = st.columns([1, 3])
    with col1:
        open_now = st.button("Open results page", type="primary", use_container_width=True)
    with col2:
        st.caption("Click to jump to the results view. If the jump fails, please open it from the left sidebar.")

    if open_now:
        try:
            st.switch_page("views/view_uploaded.py")
        except Exception:
            st.rerun()

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
        It depends on your file size and Internet speed. Please keep this tab open during processing. Normally, it takes 5~10 minutes for files under 300MB. For a faster imputation, consider using our DGAT on local machines, view our GitHub Repo <a href = 'https://github.com/osmanbeyoglulab/DGAT' >here/</a>.
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
    
    - **`adata.var`** — Gene metadata  
      *(e.g., gene names or Ensembl IDs.)*
    
    - **`adata.obsm['spatial']`** — Spatial coordinates as an **N×2 array** (x/y positions).  
    
      ⚠️ **This field is essential.** Spatial analysis cannot be performed without valid coordinates.
    
    ---
    
    ##### **3. File Size Limit**
    
    - **Recommended maximum file size:** ≤ **200 MB**  
      *(to ensure smooth and reliable upload performance.)*
    """)

