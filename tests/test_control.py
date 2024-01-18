import os

from streamlit.testing.v1 import AppTest

DEFAULT_TIMEOUT = 10


def test_fetch_abstract():
    """Get PubMed abstract for a given PMID"""

    at = AppTest.from_file("app.py", default_timeout=DEFAULT_TIMEOUT).run()
    at.text_input(key="pmid").set_value("33495752")
    at.text_input(key="pmid").run()
    at.button[0].click().run()

    with open(os.path.join("tests", "33495752.txt"), "r") as f:
        abstract = f.read()

    assert abstract in at.session_state.abstract


def test_fetch():
    at = AppTest.from_file("app.py", default_timeout=DEFAULT_TIMEOUT).run()

    with open(os.path.join("tests", "33495752.txt"), "r") as f:
        at.session_state.abstract = f.read()

    at.selectbox("llm").set_value("ChatGPT")
    at.button[1].click().run()
    assert len(at.session_state.summary) > 0

    at.selectbox("llm").set_value("BARD")
    at.button[1].click().run()
    assert len(at.session_state.summary) > 0
