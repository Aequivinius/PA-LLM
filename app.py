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
view.upload()
