import streamlit as st
import os
from dgat_utils.task_manager import Session, ProteinTask, get_s3_client

st.header("Retrieve Your Prediction Results")
st.write("Enter the 8-digit feature code sent to your email or shown after upload.")

search_code = st.text_input("Feature Code (e.g., A1B2C3D4):").strip().upper()

if st.button("Check Status", type="primary"):
    if search_code:
        session = Session()
        task = session.query(ProteinTask).filter_by(feature_code=search_code).first()
        session.close()

        if task:
            st.divider()
            # 状态展示
            status_colors = {"pending": "gray", "running": "blue", "completed": "green", "failed": "red"}
            color = status_colors.get(task.status, "gray")
            st.markdown(f"**Current Status:** :{color}[{task.status.upper()}]")

            if task.status == 'completed':
                st.success("✨ Your prediction is ready for download!")

                # 生成安全的临时下载链接 (有效期 1 小时)
                try:
                    s3 = get_s3_client()
                    download_url = s3.generate_presigned_url(
                        'get_object',
                        Params={'Bucket': os.getenv("BUCKET_NAME"), 'Key': task.output_path},
                        ExpiresIn= 900  # 15 minutes
                    )

                    st.link_button("⬇️ Download Result (.h5ad)", download_url, use_container_width=True)
                    st.caption("Note: This link is temporary and will expire in 15 min.")

                except Exception as e:
                    st.error(f"Error generating download link: {e}")

            elif task.status == 'running':
                st.info("⏳ Your data is currently being processed by our backend worker. Please check back later.")

            elif task.status == 'failed':
                st.error("❌ The computation failed. A notification with error details has been sent to your email.")
        else:
            st.error("🚫 Invalid Code. Please double-check your entry.")
    else:
        st.warning("Please enter a code.")