# main.py  —— DGAT Streamlit App (with per-session resource GC)

import os, psutil, time, uuid, threading
from typing import Callable, Any, Dict

import streamlit as st
from style import page_style, footer
import pandas as pd
from persist import load_widget_state, persist
from views.utils import get_sample_dataframe
from style import define_layout

st.set_page_config(
    page_title='DGAT',
    page_icon="./logo/gbm_ribbon.png",
    initial_sidebar_state="expanded",
)


INACTIVE_SECS = 1200
SWEEP_INTERVAL_SECS = 120

@st.cache_resource
def _resource_store() -> Dict[str, Dict[str, Dict[str, Any]]]:

    return {}

@st.cache_resource
def _lock() -> threading.RLock:
    return threading.RLock()

def _touch_heartbeat(sid: str) -> None:
    with _lock():
        store = _resource_store()
        bucket = store.get(sid, {})
        bucket["_t"] = time.time()
        store[sid] = bucket

def register_resource(sid: str, name: str, obj: Any, cleanup: Callable[[], None] | None = None) -> None:
    with _lock():
        store = _resource_store()
        bucket = store.get(sid, {"_t": time.time()})
        bucket[name] = {"obj": obj, "cleanup": cleanup}
        bucket["_t"] = time.time()
        store[sid] = bucket

def get_resource(sid: str, name: str, default=None):
    with _lock():
        store = _resource_store()
        return store.get(sid, {}).get(name, {}).get("obj", default)

def dispose_sid(sid: str) -> None:
    with _lock():
        store = _resource_store()
        bucket = store.get(sid, {})
        for k, v in list(bucket.items()):
            if k == "_t":
                continue
            cleanup = v.get("cleanup", None)
            try:
                if callable(cleanup):
                    cleanup()
            except Exception:
                pass
        if sid in store:
            del store[sid]

def _sweep_inactive() -> None:
    now = time.time()
    with _lock():
        store = _resource_store()
        dead = []
        for sid, bucket in list(store.items()):
            t = bucket.get("_t", 0)
            if now - t > INACTIVE_SECS:
                dead.append(sid)
        for sid in dead:
            bucket = store.get(sid, {})
            for k, v in list(bucket.items()):
                if k == "_t":
                    continue
                cleanup = v.get("cleanup", None)
                try:
                    if callable(cleanup):
                        cleanup()
                except Exception:
                    pass
            if sid in store:
                del store[sid]

@st.cache_resource
def _ensure_sweeper_thread() -> bool:
    def loop():
        while True:
            _sweep_inactive()
            time.sleep(SWEEP_INTERVAL_SECS)
    t = threading.Thread(target=loop, daemon=True)
    t.start()
    return True

# session id
if "_sid" not in st.session_state:
    st.session_state["_sid"] = str(uuid.uuid4())
SID = st.session_state["_sid"]

# try to start sweeper thread
_ensure_sweeper_thread()
_touch_heartbeat(SID)


process = psutil.Process(os.getpid())

if "has_upload" not in st.session_state:
    st.session_state["has_upload"] = False

define_layout(
    max_width='80%',
    padding_top='2rem',
    padding_right='0rem',
    padding_left='0rem',
    padding_bottom='0rem'
)

load_widget_state()

df_sample = get_sample_dataframe('./data/dataset.csv')
st.session_state['df_sample'] = df_sample
persist("sample_id")

emoji = "🔹"  # "🔸" "💠" ...

# page definitions
home_page = st.Page(page="views/home.py",              title="Home",               icon=emoji, default=True)
datasets_page = st.Page(page="views/dataset.py",       title="Dataset Explorer",   icon=emoji)
gene_page = st.Page(page="views/spatial_protein.py",   title="Protein Expression Maps", icon=emoji)
s_tf_page = st.Page(page="views/spatial_tf.py",        title="TF Activity Maps",   icon=emoji)
s_pathway_page = st.Page(page="views/spatial_pathway.py", title="Pathway Activity Maps", icon=emoji)
contact_page = st.Page(page="views/contact.py",        title="Contact us",         icon=emoji)
citation_page = st.Page(page="views/citation.py",      title="Citation",           icon=emoji)
termofuse_page = st.Page(page="views/termofuse.py",      title="Term of Use",           icon=emoji)
upload_page = st.Page(page="views/upload.py",          title="Upload Data",        icon=emoji)
view_uploaded_page = st.Page(page="views/view_uploaded.py", title="View your data", icon=emoji)
spatial_corr_page = st.Page(page="views/spatial_correlation.py", title="Spatial Correlation Analysis", icon=emoji)
colocalization_page = st.Page(page="views/colocalization.py", title="Colocalization Analysis", icon=emoji)
nav_groups = {
    "": [home_page],
    "Impute Your Data": [upload_page,datasets_page],
    "Resources": [termofuse_page, citation_page, contact_page],
}
if st.session_state.get("has_upload", False):
    nav_groups["Impute Your Data"].append(view_uploaded_page)
    nav_groups["Impute Your Data"].append(spatial_corr_page)
    nav_groups["Impute Your Data"].append(colocalization_page)

pg = st.navigation(nav_groups)
pg.run()

st.divider()
st.markdown(footer, unsafe_allow_html=True)
