# DEV
import os

import streamlit as st

import control
import tests.test_control
import view

st.set_page_config("Summarisation for biomedical texts", page_icon="🖐️")
control.session_states()

view.title()
view.pmid()
view.summarisation()

# DEV
with open(os.path.join("tests", f"{tests.test_control.DEFAULT_PMID}.json"), "r") as f:
    st.session_state.json = f.read()
st.session_state.pmid = tests.test_control.DEFAULT_PMID

view.upload()
