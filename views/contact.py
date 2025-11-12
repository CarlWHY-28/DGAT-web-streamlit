# import streamlit as st
######
## option 1
# a, b = st.columns([1, 1.6 ])

# a.write("")
# a.write("")
# a.write("")
# a.write("")
# a.write("")
# a.write("")
# a.subheader("Please contact:")
# a.write("")
# a.markdown("""<span style="font-size:16px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
# a.markdown("""<span style="font-size:16px;">Haoyu Wang<br>PhD Student<br>✉️ haw309@pitt.edu</span>""", unsafe_allow_html=True)
# # b.image('https://wexfordscitech.com/wp-content/uploads/2021/03/Assembly-web-5.png')
# # b.markdown("""<span style="font-size:18px;">University of Pittsburgh, UPMC Hillman Cancer Center, Assembly Building</span>""", unsafe_allow_html=True)
# b.markdown("""<div style="position:relative;padding-bottom:100%;">
#                 <iframe style="width:100%;height:100%;position:absolute;left:0px;top:0px;"
#                 frameborder="0" width="100%" height="100%" allowfullscreen allow="autoplay" 
#                src="https://maps.google.com/maps?width=100%25&amp;height=600&amp;hl=en&amp;q=5067%20CBaum%20Blvd%20Pittsburgh,%20PA%2015206+(My%20Business%20Name)&amp;t=p&amp;z=14&amp;ie=UTF8&amp;iwloc=B&amp;output=embed">
#                 <a href="https://www.maps.ie/distance-area-calculator.html">measure distance on map</a></iframe></div>""", unsafe_allow_html=True)
    

##########
# option 2, needs to Add the Secret based on the hosting service, AI tool can give the soultion


# st.markdown(
#     """
#     <style>
#       .contact-form {
#         max-width: 600px;
#         margin: auto;
#         padding: 30px;
#         background-color: #f9f9f9;
#         border-radius: 8px;
#         box-shadow: 0 2px 10px rgba(0,0,0,0.05);
#       }
#       .contact-form h1 {
#         font-size: 32px;
#         margin-bottom: 10px;
#       }
#       .contact-form p {
#         font-size: 16px;
#         color: #555;
#         margin-bottom: 20px;
#       }
#       .contact-form input, .contact-form textarea {
#         width: 100%;
#         padding: 12px;
#         margin-bottom: 15px;
#         border: 1px solid #ccc;
#         border-radius: 4px;
#         font-size: 14px;
#       }
#       .contact-form textarea {
#         height: 150px;
#         resize: vertical;
#       }
#       .contact-form button {
#         background-color: #0066cc;
#         color: white;
#         border: none;
#         padding: 14px 20px;
#         font-size: 16px;
#         border-radius: 4px;
#         cursor: pointer;
#       }
#       .contact-form button:hover {
#         background-color: #005bb5;
#       }
#     </style>
#     """,
#     unsafe_allow_html=True
# )

# # content
# st.markdown('<div class="contact-form">', unsafe_allow_html=True)
# st.markdown("<h1>Contact Us</h1>", unsafe_allow_html=True)
# st.markdown("<p>We’re here to help. Send us a message and we'll get back to you asap.</p>", unsafe_allow_html=True)

# name = st.text_input("Name *", placeholder="Your full name")
# email = st.text_input("Email *", placeholder="haw309@pitt.edu")
# subject = st.text_input("Subject", placeholder="What’s this about?")
# message = st.text_area("Message *", placeholder="Write your message here…")

# if st.button("Send Message"):
#     # handle submission (e.g., send email or store in DB)
#     st.info("Thanks for reaching out! We’ll be in touch soon.")
#     # optionally clear inputs or redirect

# st.markdown("</div>", unsafe_allow_html=True)

# # Additional contact info
# st.markdown("---")
# st.markdown("**Address:** 5607 Baum Blvd, Pittsburgh PA, 15215")  
# # st.markdown("**Phone:** (123) 456-7890")  
# # st.markdown("**Email:** support@yourdomain.com")
# st.markdown("""<span style="font-size:16px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
# st.markdown("""<span style="font-size:16px;">Haoyu Wang<br>PhD Student<br>✉️ haw309@pitt.edu</span>""", unsafe_allow_html=True)


# import requests
# # import streamlit as st

# SENDGRID_API_KEY = st.secrets["SENDGRID_API_KEY"]

# def send_email_via_sendgrid(name, email, message):
#     url = "https://api.sendgrid.com/v3/mail/send"
#     headers = {
#         "Authorization": f"Bearer {SENDGRID_API_KEY}",
#         "Content-Type": "application/json"
#     }
#     data = {
#         "personalizations": [{
#             "to": [{"email": "your_email@domain.com"}],
#             "subject": f"New contact form message from {name}"
#         }],
#         "from": {"email": "no-reply@yourdomain.com"},
#         "content": [{
#             "type": "text/plain",
#             "value": f"Name: {name}\nEmail: {email}\n\nMessage:\n{message}"
#         }]
#     }
#     r = requests.post(url, headers=headers, json=data)
#     return r.status_code

# st.title("Contact Us")

# name = st.text_input("Name *")
# email = st.text_input("Email *")
# message = st.text_area("Message *")

# if st.button("Send Message"):
#     if not name or not email or not message:
#         st.warning("Please fill in all fields.")
#     else:
#         status = send_email_via_sendgrid(name, email, message)
#         if status == 202:
#             st.success("✅ Message sent successfully!")
#         else:
#             st.error(f"❌ Error sending email (code {status})")


# ##########
# # option 3
# # use https://formspree.io/ free plan
# import streamlit as st
# import requests

# # Custom CSS styling
# st.markdown(
#     """
#     <style>
#       .contact-form {
#         max-width: 600px;
#         margin: auto;
#         padding: 30px;
#         background-color: #f9f9f9;
#         border-radius: 8px;
#         box-shadow: 0 2px 10px rgba(0,0,0,0.05);
#       }
#       .contact-form h1 {
#         font-size: 32px;
#         margin-bottom: 10px;
#       }
#       .contact-form p {
#         font-size: 16px;
#         color: #555;
#         margin-bottom: 20px;
#       }
#     </style>
#     """,
#     unsafe_allow_html=True
# )

# # Contact form
# st.markdown('<div class="contact-form">', unsafe_allow_html=True)
# st.markdown("<h1>Contact Us</h1>", unsafe_allow_html=True)
# st.markdown("<p>We're here to help. Send us a message and we'll get back to you asap.</p>", unsafe_allow_html=True)

# name = st.text_input("Name *", placeholder="Your full name")
# email = st.text_input("Email *", placeholder="your.email@example.com")
# subject = st.text_input("Subject", placeholder="What's this about?")
# message = st.text_area("Message *", placeholder="Write your message here…")

# if st.button("Send Message"):
#     # Validate inputs
#     if not name or not email or not message:
#         st.warning("⚠️ Please fill in all required fields (Name, Email, and Message).")
#     else:
#         # Send via Formspree
#         formspree_url = "https://formspree.io/f/xblqljeo"  # Get this from formspree.io
        
#         data = {
#             "name": name,
#             "email": email,
#             "subject": subject,
#             "message": message
#         }
        
#         try:
#             response = requests.post(formspree_url, data=data)
#             if response.status_code == 200:
#                 st.success("✅ Message sent successfully! We'll get back to you soon.")
#             else:
#                 st.error("❌ Something went wrong. Please try again or email us directly.")
#         except Exception as e:
#             st.error("❌ Unable to send message. Please contact us directly at the email below.")

# st.markdown("</div>", unsafe_allow_html=True)

# # Contact information
# st.markdown("---")
# st.markdown("### Direct Contact Information")
# st.markdown("**Address:** 5607 Baum Blvd, Pittsburgh PA, 15215")
# st.markdown("")
# st.markdown(
#     """<span style="font-size:16px;">
#     <strong>Hatice Osmanbeyoglu</strong><br>
#     Principal Investigator<br>
#     ✉️ osmanbeyogluhu@pitt.edu
#     </span>""", 
#     unsafe_allow_html=True
# )
# st.markdown("")
# st.markdown(
#     """<span style="font-size:16px;">
#     <strong>Haoyu Wang</strong><br>
#     PhD Student<br>
#     ✉️ haw309@pitt.edu
#     </span>""", 
#     unsafe_allow_html=True
# )

##########
# option 4
import streamlit as st
from urllib.parse import quote

# Custom CSS styling
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
      .email-button {
        display: inline-block;
        background-color: #0066cc;
        color: white;
        padding: 14px 28px;
        text-decoration: none;
        border-radius: 4px;
        font-size: 16px;
        margin-top: 10px;
        transition: background-color 0.3s;
      }
      .email-button:hover {
        background-color: #005bb5;
      }
    </style>
    """,
    unsafe_allow_html=True
)

# Contact form
st.markdown('<div class="contact-form">', unsafe_allow_html=True)
st.markdown("<h1>Contact Us</h1>", unsafe_allow_html=True)
st.markdown("<p>We're here to help. Send us a message and we'll get back to you asap.</p>", unsafe_allow_html=True)

name = st.text_input("Name *", placeholder="Your full name")
email = st.text_input("Your Email *", placeholder="your.email@example.com")
subject = st.text_input("Subject", placeholder="What's this about?")
message = st.text_area("Message *", placeholder="Write your message here…")

if st.button("Send Message"):
    if not name or not email or not message:
        st.warning("⚠️ Please fill in all required fields.")
    else:
        # Create mailto link
        to_email = "xim33@pitt.edu"
        email_subject = quote(subject or "Contact Form Submission")
        email_body = quote(f"Name: {name}\nEmail: {email}\n\n{message}")
        mailto_link = f"mailto:{to_email}?subject={email_subject}&body={email_body}"
        
        st.success("✅ Opening your email client...")
        st.markdown(
            f'<a href="{mailto_link}" class="email-button" target="_blank">Click here if email didn\'t open</a>',
            unsafe_allow_html=True
        )
        
        # Auto-redirect with JavaScript
        st.markdown(
            f'<script>window.location.href = "{mailto_link}";</script>',
            unsafe_allow_html=True
        )

st.markdown("</div>", unsafe_allow_html=True)

# Contact information
st.markdown("---")
st.markdown("### Direct Contact Information")
st.markdown("**Address:** 5607 Baum Blvd, Pittsburgh PA, 15215")
st.markdown("")
st.markdown(
    """<span style="font-size:16px;">
    <strong>Hatice Osmanbeyoglu</strong><br>
    Principal Investigator<br>
    ✉️ <a href="mailto:osmanbeyogluhu@pitt.edu">osmanbeyogluhu@pitt.edu</a>
    </span>""", 
    unsafe_allow_html=True
)
st.markdown("")
st.markdown(
    """<span style="font-size:16px;">
    <strong>Haoyu Wang</strong><br>
    PhD Student<br>
    ✉️ <a href="mailto:haw309@pitt.edu">haw309@pitt.edu</a>
    </span>""", 
    unsafe_allow_html=True
)
