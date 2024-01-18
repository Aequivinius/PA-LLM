from streamlit.testing.v1 import AppTest


def test_fetch_abstract_button():
    """User can only click fetch button when they enter a number"""
    at = AppTest.from_file("app.py").run()
    at.text_input(key="pmid").set_value("not digits")
    at.text_input(key="pmid").run()
    assert at.button[0].disabled
    at.text_input(key="pmid").set_value("33495752")
    at.text_input(key="pmid").run()
    assert not at.button[0].disabled


def test_fetch_summary_button():
    """Can only be clicked if abstract is there"""
    at = AppTest.from_file("app.py").run()
    assert at.button[1].disabled
    at.session_state.abstract = "Lorem ipsum"
    at.run()
    assert not at.button[1].disabled
