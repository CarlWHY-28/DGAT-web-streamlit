import streamlit as st


a, b = st.columns([1, 1.6 ])

a.write("")
a.write("")
a.write("")
a.write("")
a.write("")
a.write("")
a.write("")
a.subheader("Please contact:")

a.markdown("""<span style="font-size:16px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
a.markdown("""<span style="font-size:16px;">Haoyu Wang<br>PhD Student<br>✉️ haw309@pitt.edu</span>""", unsafe_allow_html=True)
# b.image('https://wexfordscitech.com/wp-content/uploads/2021/03/Assembly-web-5.png')
# b.markdown("""<span style="font-size:18px;">University of Pittsburgh, UPMC Hillman Cancer Center, Assembly Building</span>""", unsafe_allow_html=True)
b.markdown("""<div style="position:relative;padding-bottom:100%;">
                <iframe style="width:100%;height:100%;position:absolute;left:0px;top:0px;"
                frameborder="0" width="100%" height="100%" allowfullscreen allow="autoplay" 
               src="https://maps.google.com/maps?width=100%25&amp;height=600&amp;hl=en&amp;q=5051%20Centre%20Avenue%20Pittsburgh,%20PA%2015213+(My%20Business%20Name)&amp;t=p&amp;z=14&amp;ie=UTF8&amp;iwloc=B&amp;output=embed">
                <a href="https://www.maps.ie/distance-area-calculator.html">measure distance on map</a></iframe></div>""", unsafe_allow_html=True)
    
st.markdown("#")



# import streamlit as st

# page config
# st.set_page_config(page_title="Contact Us", layout="centered")

# some custom CSS for styling
st.markdown(
    """
    <style>
      .contact-form {
        max-width: 600px;
        margin: auto;
        padding: 30px;
        background-color: #f9f9f9;
        border-radius: 8px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
      }
      .contact-form h1 {
        font-size: 32px;
        margin-bottom: 10px;
      }
      .contact-form p {
        font-size: 16px;
        color: #555;
        margin-bottom: 20px;
      }
      .contact-form input, .contact-form textarea {
        width: 100%;
        padding: 12px;
        margin-bottom: 15px;
        border: 1px solid #ccc;
        border-radius: 4px;
        font-size: 14px;
      }
      .contact-form textarea {
        height: 150px;
        resize: vertical;
      }
      .contact-form button {
        background-color: #0066cc;
        color: white;
        border: none;
        padding: 14px 20px;
        font-size: 16px;
        border-radius: 4px;
        cursor: pointer;
      }
      .contact-form button:hover {
        background-color: #005bb5;
      }
    </style>
    """,
    unsafe_allow_html=True
)

# content
st.markdown('<div class="contact-form">', unsafe_allow_html=True)
st.markdown("<h1>Contact Us</h1>", unsafe_allow_html=True)
st.markdown("<p>We’re here to help. Send us a message and we'll get back to you asap.</p>", unsafe_allow_html=True)

name = st.text_input("Name *", placeholder="Your full name")
email = st.text_input("Email *", placeholder="you@example.com")
subject = st.text_input("Subject", placeholder="What’s this about?")
message = st.text_area("Message *", placeholder="Write your message here…")

if st.button("Send Message"):
    # handle submission (e.g., send email or store in DB)
    st.info("Thanks for reaching out! We’ll be in touch soon.")
    # optionally clear inputs or redirect

st.markdown("</div>", unsafe_allow_html=True)

# Additional contact info
st.markdown("---")
st.markdown("**Address:** 5607 Baum Blvd, Pittsburgh PA, 15215")  
# st.markdown("**Phone:** (123) 456-7890")  
# st.markdown("**Email:** support@yourdomain.com")
st.markdown("""<span style="font-size:16px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
st.markdown("""<span style="font-size:16px;">Haoyu Wang<br>PhD Student<br>✉️ haw309@pitt.edu</span>""", unsafe_allow_html=True)
