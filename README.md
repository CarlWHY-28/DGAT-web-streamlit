# 🧬 DGAT Protein Inference Platform

### Introduction
This is a [**Protein Sequence Inference Platform**](https://dgat-protein-imputation.up.railway.app) based on Deep Graph Attention Networks ([DGAT_Repo](https://github.com/CarlWHY-28/DGAT)).

---

### How It Works

The platform separates the **User Interface** from the **Calculation** to ensure stability and speed:

1.  **Submit Task**: You upload your `.h5ad` data file on the web page. The system saves the file to a storage bucket and gives you a unique **Tracking Code**. You can then close the webpage.
2.  **Background Processing**: A dedicated AI Worker automatically downloads your data and runs the complex DGAT deep learning model in the background.
3.  **Email Notification**: Once the inference is complete, the system sends you an email notification with the **Tracking Code**.
4.  **Retrieve Results**: You return to the website, enter your **Tracking Code**, and instantly view the visualization charts or download the result files.

---

### Tech Stack

* **Frontend**: Streamlit
* **Backend**: Python (Custom Inference Worker)
* **Database**: MySQL (Task Status Management)
* **Storage**: S3 Compatible Bucket (File Storage)
* **Deployment**: Railway Cloud

---
>https://dgat-protein-imputation.up.railway.app
> 
> *Haoyu Wang - 2026*
---
### 01/06/2025 Update
