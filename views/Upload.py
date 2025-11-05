import streamlit as st
import boto3
from botocore.exceptions import NoCredentialsError
import io
from datetime import datetime
import smtplib
from email.mime.text import MIMEText
from email.utils import formataddr
from dgat_utils.predict_util import web_predict
import anndata as ad
import traceback
import tempfile
import os

url_REPO = 'https://raw.githubusercontent.com/CarlWHY-28/DGAT-web-resource/main'

def send_email_smtp(to_addr: str, subject: str, body: str):
    """
    Gmail SMTP + App Password
    """
    gmail_cfg = st.secrets.get("gmail", {})
    gmail_user = gmail_cfg.get("gmail_user")
    gmail_app_password = gmail_cfg.get("gmail_app_password")

    if not gmail_user or not gmail_app_password:
        st.error("Lost Gmail config: Check secrets.toml.")
        return False

    msg = MIMEText(body, "plain", "utf-8")
    msg["From"] = formataddr(("DGAT Portal", gmail_user))
    msg["To"] = to_addr
    msg["Subject"] = subject

    try:
        with smtplib.SMTP("smtp.gmail.com", 587) as server:
            server.ehlo()
            server.starttls()
            server.login(gmail_user, gmail_app_password)
            server.sendmail(gmail_user, [to_addr], msg.as_string())
        return True
    except Exception as e:
        st.error(f"Error: {e}")
        return False


# Config
try:
    s3_client = boto3.client(
        "s3",
        aws_access_key_id=st.secrets["aws"]["aws_access_key_id"],
        aws_secret_access_key=st.secrets["aws"]["aws_secret_access_key"],
        region_name=st.secrets["aws"]["aws_region"]
    )
    BUCKET_NAME = st.secrets["aws"]["s3_bucket_name"]
    s3_setup_success = True

except KeyError:
    st.error("S3 credentials failed! Check .streamlit/secrets.toml")
    s3_setup_success = False
    st.stop()
except NoCredentialsError:
    st.error("AWS credentials not found")
    s3_setup_success = False
    st.stop()
except Exception as e:
    st.error(f"Unexpected Error {e}")
    s3_setup_success = False
    st.stop()


def upload_to_s3(file_buffer, file_name, first_name, email):
    if not s3_setup_success:
        st.error("S3 configuration is not set up correctly.")
        return False, None

    try:
        current_time = st.session_state.get('current_time', 'notime')
        object_name = f"{first_name.replace(' ', '')}_{email}_{current_time}_{file_name}"

        file_buffer.seek(0)
        s3_client.upload_fileobj(file_buffer, BUCKET_NAME, object_name)

        st.success(f"File '{object_name}' uploaded！")
        return True, object_name

    except Exception as e:
        st.error(f"Upload failed: {e}")
        return False, None


st.markdown("<h2 style='text-align: center; color: black;'>Upload Your Data</h1>", unsafe_allow_html=True)
st.write("")
st.info(
    "Upload your spatial transcriptomics data in h5ad format. Ensure your data is correctly formatted for analysis.")

with st.form(key="upload_form"):
    first_name = st.text_input("First Name", key="first_name")
    email = st.text_input("Email", key="email")
    uploaded_file = st.file_uploader("Choose an h5ad file", type="h5ad", key="upload_file")
    submit_button = st.form_submit_button(label="Upload and Submit")

if submit_button:
    if uploaded_file is not None and first_name and email:
        now = datetime.now()
        timestamp_str = now.strftime("%Y%m%d-%H%M%S")
        st.session_state['current_time'] = timestamp_str

        with st.spinner(f"Uploading '{uploaded_file.name}' ..."):
            upload_success, object_name = upload_to_s3(
                file_buffer=uploaded_file,
                file_name=uploaded_file.name,
                first_name=first_name,
                email=email
            )

            if upload_success:
                st.session_state['data_uploaded'] = True
                st.balloons()

                admin_email = st.secrets["gmail"].get("admin_email", None)
                user_email = email

                subject_user = "DGAT Portal: Your data has been received"
                body_user = (
                    f"Hi {first_name},\n\n"
                    f"Your file has been successfully uploaded.\n"
                    # f"- S3 bucket: {BUCKET_NAME}\n"
                    # f"- Object key: {object_name}\n"
                    # f"- Upload time: {timestamp_str}\n\n"
                    f"We'll start processing your data shortly. Thank you!"
                )

                subject_admin = "DGAT Portal: New upload received"
                body_admin = (
                    f"A new file has been uploaded.\n\n"
                    f"- Uploader: {first_name} <{user_email}>\n"
                    f"- S3 bucket: {BUCKET_NAME}\n"
                    f"- Object key: {object_name}\n"
                    f"- Upload time: {timestamp_str}\n"
                    f"- Original filename: {uploaded_file.name}"
                )

                ok_user = send_email_smtp(user_email, subject_user, body_user)
                if ok_user:
                    st.info(f"Confirmation sent: {user_email}")

                if admin_email:
                    ok_admin = send_email_smtp(admin_email, subject_admin, body_admin)
                    if ok_admin:
                        st.info(f"Admin informed: {admin_email}")

                # Prediction
                try:
                    uploaded_file.seek(0)
                    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_in:
                        tmp_in.write(uploaded_file.getbuffer())
                        tmp_in_path = tmp_in.name

                    with st.status("Running Imputation (3~5 min)", state="running") as status:
                        status.update(label="Loading", state="running")
                        adata_in = ad.read_h5ad(tmp_in_path)

                        status.update(label="Predicting (Don't close the page)", state="running")
                        adata_out = web_predict(url_REPO,adata_in)
                        print("Prediction done.")
                        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_out:
                            adata_out.write_h5ad(tmp_out.name)
                            final_out_path = tmp_out.name

                        status.update(label="Getting ready for download...", state="complete")

                    # 回到内存，生成下载按钮
                    with open(final_out_path, "rb") as f_out:
                        result_bytes = f_out.read()

                    download_name = f"{os.path.splitext(uploaded_file.name)[0]}_DGAT_pred_{timestamp_str}.h5ad"
                    st.success("Prediction finished, you can download the result below:")
                    st.download_button(
                        label="⬇️ Download Predicted Data",
                        data=result_bytes,
                        file_name=download_name,
                        mime="application/octet-stream"
                    )

                    #Visualization link

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
        st.warning("Please fill in all fields and select a file before submitting.")

st.write("")
#st.markdown("")

st.divider()
tab1, tab2 = st.tabs(["Frequently Asked Questions (FAQ)", "Data Requirements"])

with tab1:
    st.subheader("Frequently Asked Questions (FAQ)")
    st.markdown("""
        **Q1: What is this upload page for?**
        A: This page is used to upload your spatial transcriptomics data (in `h5ad` format) for protein imputation.

        **Q2: What information do I need to provide?**
        A: Please provide your First Name and Email, and select an `.h5ad` file. All fields are mandatory.

        **Q3: Why do you need my name and email?**
        A: This information is used to uniquely identify your data file. The uploaded filename will include this information (and a timestamp) to prevent conflicts and help track data provenance.

        **Q4: Is my data secure?**
        A: Yes. Your data is uploaded to a private Amazon S3 bucket. While the filename is processed, the file content itself is kept confidential.

        **Q5: How long does the upload take?**
        A: Upload time depends on your file size and internet connection speed. Please keep this browser window open during the upload process.

        **Q6: What if the upload fails?**
        A: Please check your internet connection and ensure your file is a valid `.h5ad` format. If the problem persists, please contact the administrator.
    """)

with tab2:
    st.subheader("Data Requirements")
    st.markdown("""
        To ensure the analysis tools can process your data correctly, please make sure your `.h5ad` file meets the following requirements:

        **1. File Format:**
        * Must be in `h5ad` (Anndata) format.

        **2. Required Data Slots:**
        * `adata.X`: Raw count matrix. We recommend using non-normalized raw counts. Preprocessing will be handled by our pipeline. 
        * `adata.obs`: Cell/Spot metadata.
        * `adata.var`: Gene metadata (e.g., gene names).
        * `adata.obsm['spatial']`: **(CRITICAL)** Spatial coordinates. This must be a NumPy array where each row corresponds to a spot/cell and the two columns represent (x, y) coordinates. **Spatial analysis is impossible without this.**

        **3. File Size:**
        * Please keep the file size under 200MB for optimal upload performance.
    """)
