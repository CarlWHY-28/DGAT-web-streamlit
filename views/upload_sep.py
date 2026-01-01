import streamlit as st
import uuid
import os
from dgat_utils.task_manager import Session, ProteinTask, get_s3_client

st.header("Upload and Submit Your Task")
st.info("Upload your `.h5ad` file and we will notify you via email once computation is done.")

with st.form("upload_form"):
    email = st.text_input("Email Address (for notification)")
    uploaded_file = st.file_uploader("Choose an h5ad file", type="h5ad")
    submit = st.form_submit_button("Submit Task", type="primary")

if submit:
    if uploaded_file and email:
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