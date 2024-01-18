import json
import os

from streamlit.testing.v1 import AppTest

import control

DEFAULT_TIMEOUT = 10
DEFAULT_PMID = "33495752"


def test_fetch_abstract():
    """Get PubMed abstract for a given PMID"""

    at = AppTest.from_file("app.py", default_timeout=DEFAULT_TIMEOUT).run()
    at.text_input(key="pmid").set_value(DEFAULT_PMID)
    at.text_input(key="pmid").run()
    at.button[0].click().run()

    with open(os.path.join("tests", f"{DEFAULT_PMID}.txt"), "r") as f:
        abstract = f.read()

    assert abstract in at.session_state.abstract


def test_fetch():
    at = AppTest.from_file("app.py", default_timeout=DEFAULT_TIMEOUT).run()

    with open(os.path.join("tests", f"{DEFAULT_PMID}.txt"), "r") as f:
        at.session_state.abstract = f.read()

    at.selectbox("llm").set_value("ChatGPT")
    at.button[1].click().run()
    assert len(at.session_state.summary) > 0

    at.selectbox("llm").set_value("BARD")
    at.button[1].click().run()
    assert len(at.session_state.summary) > 0


def test_jsonify():
    with open(os.path.join("tests", f"{DEFAULT_PMID}.txt"), "r") as f:
        text = f.read()

    with open(os.path.join("tests", f"{DEFAULT_PMID}_summary.txt"), "r") as f:
        summary = f.read()

    jsonified = control.jsonify_(text, sourceid=DEFAULT_PMID, summary=summary)

    with open(os.path.join("tests", "33495752.json"), "r") as f:
        jsonified_ = json.dumps(json.loads(f.read()))

    assert jsonified == jsonified_
