import streamlit as st

import constants
import control


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
            disabled=st.session_state.pmid is constants.SESSION_STATES["pmid"],
        )

    if st.session_state.abstract:
        st.success(st.session_state.abstract)

    if st.session_state.fetch_abstract_error:
        st.error(st.session_state.fetch_abstract_error)


@block
def summarisation():
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
            disabled=st.session_state.abstract is constants.SESSION_STATES["abstract"],
            on_change=control.clear_summary,
        )

    with right:
        st.button(
            "Get summary",
            on_click=control.fetch_summary,
            use_container_width=True,
            type="secondary",
            disabled=st.session_state.abstract is constants.SESSION_STATES["abstract"],
        )

    if st.session_state.summary:
        st.success(st.session_state.summary)
        control.jsonify()


@block
def upload():
    st.markdown("*Use the summary:*")

    left, right = st.columns(2)
    with left:
        st.button(
            "⬆️ Let's upload summaries to PubAnnotation",
            on_click=control.upload,
            type="primary",
            use_container_width=True,
            disabled=st.session_state.json is constants.SESSION_STATES["json"],
        )

    with right:
        st.download_button(
            "⬇️ Let's download the `.json`",
            data=st.session_state.json,
            key="jsonify",
            file_name=f"{st.session_state.pmid}_summary.json",
            mime="application/json",
            use_container_width=True,
            disabled=st.session_state.json is constants.SESSION_STATES["json"],
        )

    if (
        st.session_state.upload
        and st.session_state.summary in st.session_state.upload["blocks"][0]["obj"]
    ):
        st.success(
            f"Uploaded! Go check it out at **[pubannotation.org](https://pubannotation.org/projects/PA-LLM/docs/sourcedb/PubMed/sourceid/{st.session_state.pmid})**! 👏"
        )
