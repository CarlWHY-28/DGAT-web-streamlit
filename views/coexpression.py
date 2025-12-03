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
st.write("""Edge weights represent the Spearman's rank correlation coefficients (SRCCs) between different pairs of proteins across all spots.
    Hover over an edge to view a specific SRCC.""")

net = Network(height = "600px", width = "100%", bgcolor = "#222222", font_color = "white")

n_nodes = srcc_matrix.shape[0]

for node in srcc_matrix.columns:
    net.add_node(node, size = 5)

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
net.save_graph("coexpression_graph.html")
with open("coexpression_graph.html", "r", encoding="utf-8") as f:
    components.html(f.read(), height=650)

dt_df = pd.DataFrame({
    'Protein': ['ACTA2', 'BCL2', 'CCR7', 'CD14', 'CD163', 'CD19', 'CD27', 'CD274', 'CD3E', 'CD4', 'CD40', 'CD68', 
                'CD8A', 'CEACAM8', 'CR2', 'CXCR5', 'EPCAM', 'FCGR3A', 'HLA_DRA', 'ITGAM', 'ITGAX', 'KRT5', 'MS4A1', 
                'PAX5', 'PCNA', 'PDCD1', 'PECAM1', 'PTPRC_1', 'PTPRC_2', 'SDC1', 'VIM'], 
    'UniProt ID': ['P62736', 'P10415', 'P32248', 'P08571', 'Q86VB7', 'P15391', 'P26842', 'Q9NZQ7', 'P07766', 'P01730', 
                   'P25942', 'P34810', 'P01732', 'P31997', 'P20023', 'P32302', 'P16422', 'P08637', 'P01903', 'P11215', 
                   'P20702', 'P13647', 'P11836', 'Q02548', 'P12004', 'Q15116', 'P16284', 'P08575', 'P08575', 'P18827', 
                   'P08670'],
    'TTD Drug Target ID': ['Not a known drug target',
 "<a href='https://idrblab.net/ttd/data/target/details/t31309' target='_blank'>T31309</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t84981' target='_blank'>T84981</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t23212' target='_blank'>T23212</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t36576' target='_blank'>T36576</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t56365' target='_blank'>T56365</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t85554' target='_blank'>T85554</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t99948' target='_blank'>T99948</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t87075' target='_blank'>T87075</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t10191' target='_blank'>T10191</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t45758' target='_blank'>T45758</a>",
 'Not a known drug target',
 "<a href='https://idrblab.net/ttd/data/target/details/t20978' target='_blank'>T20978</a>",
 'Not a known drug target',
 "<a href='https://idrblab.net/ttd/data/target/details/t18059' target='_blank'>T18059</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t09528' target='_blank'>T09528</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t47863' target='_blank'>T47863</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t82240' target='_blank'>T82240</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t49146' target='_blank'>T49146</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t95616' target='_blank'>T95616</a>",
 'Not a known drug target',
 'Not a known drug target',
 "<a href='https://idrblab.net/ttd/data/target/details/t73215' target='_blank'>T73215</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t29818' target='_blank'>T29818</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t21782' target='_blank'>T21782</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t59631' target='_blank'>T59631</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t89056' target='_blank'>T89056</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t43115' target='_blank'>T43115</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t43115' target='_blank'>T43115</a>",
 "<a href='https://idrblab.net/ttd/data/target/details/t13017' target='_blank'>T13017</a>",
 'Not a known drug target']
})

st.markdown("<h2 style='text-align: center; color: black;'>TTD Drug Target IDs</h2>", unsafe_allow_html=True)
st.write("Drug target IDs are from the [Therapeutic Target Database](https://ttd.idrblab.cn/) (TTD).")

html_table = dt_df.to_html(escape = False)  # escape = False keeps HTML tags
st.markdown(html_table, unsafe_allow_html = True)