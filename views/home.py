import streamlit as st


#import streamlit_analytics

#streamlit_analytics.start_tracking()
#streamlit_analytics.stop_tracking()

img = "./images/DGAT_Training.png"


st.markdown("<h2 style='text-align: center;  color: black;'>DGAT: A Dual-Graph Attention Network for Inferring Spatial Protein Landscapes from Transcriptomics", unsafe_allow_html=True)

st.write("")

st.image(img)

st.write("")
st.write("")

# st.markdown('➡️ Read our paper here: https://www.google.com')   

desc1 = "Spatial Transcriptomics (ST) technologies have revolutionized our understanding of the tumor microenvironment (TME) by enabling the mapping of gene expression within the spatial context of tissues. However, mRNA levels do not always correlate with protein abundance, which is crucial for understanding cellular functions and interactions. To address this gap, we developed DGAT (Dual-Graph Attention Network), a novel computational framework that integrates gene profiles and spatial information to predict protein expression from ST data."
desc2 = "In this study, we applied DGAT to publicly available ST datasets and additionally generated analysis based on predicted protein profiles."

st.markdown(f"<p style='text-align: justify; color: black; font-size:18px'>{desc1}</p>", unsafe_allow_html=True) 
st.markdown(f"<p style='text-align: justify; color: black; font-size:18px'>{desc2}</p>", unsafe_allow_html=True) 

st.markdown('📚 Read our preprint here! [DGAT: A Dual-Graph Attention Network for Inferring Spatial Protein Landscapes from Transcriptomics](https://www.biorxiv.org/content/10.1101/2025.07.05.662121v1)')

# st.markdown(' Read our paper from here 👉 https://aacrjournals.org/cancerrescommun/article/4/8/2133/747011/Spatial-Landscape-of-Malignant-Pleural-and')    
