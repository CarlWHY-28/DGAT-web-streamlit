# Adapted from a template provided by Microsoft Copilot

import streamlit as st
from pyvis.network import Network
import streamlit.components.v1 as components
import pandas as pd

adata_out = st.session_state.get("adata_out", None)

st.markdown("<h2 style='text-align: center; color: black;'>Protein co-expression graph</h2>", unsafe_allow_html=True)
st.write("")

adata_pandas = pd.DataFrame(adata_out.X, columns = adata_out.var_names)

# Calculate Spearman's rank correlation coefficient (srcc) between the expression of all protein pairs
srcc_matrix = adata_pandas.corr(method = "spearman")

threshold = st.slider("Show edges with Spearman's rank correlation coefficients in this range:", min_value = -1.0, max_value = 1.0, value = (0.8, 1.0), step = 0.05)

st.write("Selected range:", threshold)
st.write("Click on a protein to see its drug target ID from the [Therapeutic Target Database](https://ttd.idrblab.cn/) and the corresponding link.")

net = Network(height = "600px", width = "100%", bgcolor = "#222222", font_color = "white")

n_nodes = srcc_matrix.shape[0]

# Links to drug target IDs from TTD
title_string_dict = {'BCL2': "<a href='https://idrblab.net/ttd/data/target/details/t31309' target='_blank'>T31309</a>",
 'CCR7': "<a href='https://idrblab.net/ttd/data/target/details/t84981' target='_blank'>T84981</a>",
 'CD14': "<a href='https://idrblab.net/ttd/data/target/details/t23212' target='_blank'>T23212</a>",
 'CD163': "<a href='https://idrblab.net/ttd/data/target/details/t36576' target='_blank'>T36576</a>",
 'CD19': "<a href='https://idrblab.net/ttd/data/target/details/t56365' target='_blank'>T56365</a>",
 'CD27': "<a href='https://idrblab.net/ttd/data/target/details/t85554' target='_blank'>T85554</a>",
 'CD274': "<a href='https://idrblab.net/ttd/data/target/details/t99948' target='_blank'>T99948</a>",
 'CD3E': "<a href='https://idrblab.net/ttd/data/target/details/t87075' target='_blank'>T87075</a>",
 'CD4': "<a href='https://idrblab.net/ttd/data/target/details/t10191' target='_blank'>T10191</a>",
 'CD40': "<a href='https://idrblab.net/ttd/data/target/details/t45758' target='_blank'>T45758</a>",
 'CD8A': "<a href='https://idrblab.net/ttd/data/target/details/t20978' target='_blank'>T20978</a>",
 'CR2': "<a href='https://idrblab.net/ttd/data/target/details/t18059' target='_blank'>T18059</a>",
 'CXCR5': "<a href='https://idrblab.net/ttd/data/target/details/t09528' target='_blank'>T09528</a>",
 'EPCAM': "<a href='https://idrblab.net/ttd/data/target/details/t47863' target='_blank'>T47863</a>",
 'FCGR3A': "<a href='https://idrblab.net/ttd/data/target/details/t82240' target='_blank'>T82240</a>",
 'HLA_DRA': "<a href='https://idrblab.net/ttd/data/target/details/t49146' target='_blank'>T49146</a>",
 'ITGAM': "<a href='https://idrblab.net/ttd/data/target/details/t95616' target='_blank'>T95616</a>",
 'MS4A1': "<a href='https://idrblab.net/ttd/data/target/details/t73215' target='_blank'>T73215</a>",
 'PAX5': "<a href='https://idrblab.net/ttd/data/target/details/t29818' target='_blank'>T29818</a>",
 'PCNA': "<a href='https://idrblab.net/ttd/data/target/details/t21782' target='_blank'>T21782</a>",
 'PDCD1': "<a href='https://idrblab.net/ttd/data/target/details/t59631' target='_blank'>T59631</a>",
 'PECAM1': "<a href='https://idrblab.net/ttd/data/target/details/t89056' target='_blank'>T89056</a>",
 'PTPRC_1': "<a href='https://idrblab.net/ttd/data/target/details/t43115' target='_blank'>T43115</a>",
 'PTPRC_2': "<a href='https://idrblab.net/ttd/data/target/details/t43115' target='_blank'>T43115</a>",
 'SDC1': "<a href='https://idrblab.net/ttd/data/target/details/t13017' target='_blank'>T13017</a>"}

for node in srcc_matrix.columns:
    if node in title_string_dict:
        title_string = title_string_dict[node]
    else:
        title_string = "Not a known drug target"
    
    net.add_node(node, size = 5, title = title_string)

for i in range(n_nodes):
    row_id = srcc_matrix.columns[i]
    for j in range(i + 1, n_nodes):
        col_id = srcc_matrix.columns[j]
        weight = round(srcc_matrix.loc[row_id, col_id], 2)
        title = f"SRCC: {weight}"
        if weight >= threshold[0] and weight <= threshold[1]:
            net.add_edge(row_id, col_id, title = title, width = 1, smooth = False)

net.toggle_physics(False)

# Save and render in Streamlit
net.save_graph("graph.html")
with open("graph.html", "r", encoding="utf-8") as f:
    components.html(f.read(), height=650)