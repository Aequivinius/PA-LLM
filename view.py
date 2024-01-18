import streamlit as st

import constants
import control

st.set_page_config("Summarisation for biomedical texts", page_icon="🖐️")


def block(f):
    def wrapper():
        with st.spinner(f"Executing {f.__name__}."):
            f()
        st.divider()

    return wrapper


@block
def title():
    st.title("🖐️ PA-LLM")


@block
def pmid():
    """Allows you to enter a number as PMID and issue fetching the article"""
    st.markdown("*Enter a PMID to fetch the abstract:*")

    left, right = st.columns(2)
    with left:
        st.text_input(
            "Enter a PMID:",
            label_visibility="collapsed",
            key="pmid",
            on_change=control.validate_pmid,
            placeholder="for example: 33495752",
        )

        if st.session_state.pmid_validation_error:
            st.warning(st.session_state.pmid_validation_error)

    with right:
        st.button(
            "Fetch abstract",
            on_click=control.fetch_abstract,
            use_container_width=True,
            type="primary",
            disabled=st.session_state.pmid is None,
        )

    if st.session_state.abstract:
        st.success(st.session_state.abstract)

    if st.session_state.fetch_abstract_error:
        st.error(st.session_state.fetch_abstract_error)


@block
def summarisation() -> None:
    """Select LLM and request summary"""

    st.markdown(
        "*Select with LLM you want to use for summarisation and simplification:*"
    )

    left, right = st.columns(2)
    with left:
        st.selectbox(
            "Which LLM to use for summarisation?",
            constants.LLMS_,
            label_visibility="collapsed",
            key="llm",
            disabled=st.session_state.abstract is None,
            on_change=control.clear_summary,
        )

    with right:
        st.button(
            "Get summary",
            on_click=control.fetch_summary,
            use_container_width=True,
            type="secondary",
            disabled=st.session_state.abstract is None,
        )

    if st.session_state.summary:
        st.success(st.session_state.summary)
