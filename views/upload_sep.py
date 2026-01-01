import streamlit as st
import uuid
import os
from dgat_utils.task_manager import Session, ProteinTask, get_s3_client

st.set_page_config(page_title="DGAT Protein Imputation", layout="wide")

tab1, tab2 = st.tabs(["📤 Upload & Submit", "🔍 Check Status & Result"])

# --- Tab 1: 上传逻辑 ---
with tab1:
    st.header("Upload Your Data")
    with st.form("upload_form"):
        email = st.text_input("Email Address (for notification)")
        uploaded_file = st.file_uploader("Choose an h5ad file", type="h5ad")
        submit = st.form_submit_button("Submit Task")

    if submit:
        if uploaded_file and email:
            feature_code = str(uuid.uuid4())[:8].upper()  # 简短特征码
            bucket_name = os.getenv("BUCKET_NAME")
            # 建立独立文件夹结构: task_{code}/input.h5ad
            input_key = f"task_{feature_code}/input.h5ad"

            try:
                # 1. 上传 Bucket
                s3 = get_s3_client()
                s3.upload_fileobj(uploaded_file, bucket_name, input_key)

                # 2. 写入数据库
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

                st.success(f"Successfully submitted! Your Feature Code: **{feature_code}**")
                st.info("We will notify you via email when results are ready.")
            except Exception as e:
                st.error(f"Submission failed: {e}")
        else:
            st.warning("Please provide both file and email.")

# --- Tab 2: 查询逻辑 ---
with tab2:
    st.header("Retrieve Your Results")
    search_code = st.text_input("Enter your Feature Code:")
    if st.button("Query"):
        session = Session()
        task = session.query(ProteinTask).filter_by(feature_code=search_code).first()
        session.close()

        if task:
            st.write(f"**Status:** {task.status.capitalize()}")
            if task.status == 'completed':
                st.success("Result is ready!")
                # 这里可以扩展展示代码或下载代码
                st.write(f"Result file path in Bucket: `{task.output_path}`")
            elif task.status == 'running':
                st.info("The worker is currently processing your data...")
            elif task.status == 'failed':
                st.error("Prediction failed. Please check your data.")
        else:
            st.warning("Code not found.")