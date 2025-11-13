import streamlit as st

a, b, c, d = st.columns([3, 1, 1, 2 ])

a.write("")
a.write("")


a.subheader("Please contact:")
a.write("")
c.write("")
c.write("")
c.image('logo/dbmi.png', use_column_width=True)

a2, b2, c2 = st.columns([3,2,2])
a2.markdown("""<span style="font-size:16px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
a2.markdown("""<span style="font-size:16px;">Haoyu Wang<br>PhD Student<br>✉️ haw309@pitt.edu</span>""", unsafe_allow_html=True)

b2.markdown("<p style='text-align: center;'>University of Pittsburgh<br>Department of Biomedical Informatics</p>", unsafe_allow_html=True)
# b.markdown("<p style='text-align: center;'>Department of Biomedical Informatics</p>", unsafe_allow_html=True)
1
# b.markdown("""<div style="position:relative;padding-bottom:100%;">
#                 <iframe style="width:100%;height:100%;position:absolute;left:0px;top:0px;"
#                 frameborder="0" width="100%" height="100%" allowfullscreen allow="autoplay" 
#                src="https://maps.google.com/maps?width=100%25&amp;height=600&amp;hl=en&amp;q=5067%20CBaum%20Blvd%20Pittsburgh,%20PA%2015206+(My%20Business%20Name)&amp;t=p&amp;z=14&amp;ie=UTF8&amp;iwloc=B&amp;output=embed">
#                 <a href="https://www.maps.ie/distance-area-calculator.html">measure distance on map</a></iframe></div>""", unsafe_allow_html=True)
    

