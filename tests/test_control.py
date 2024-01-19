import json
import os
import random
import string

import requests
from dotenv import load_dotenv
from requests.auth import HTTPBasicAuth
from streamlit.testing.v1 import AppTest

import control

DEFAULT_TIMEOUT = 10
DEFAULT_PMID = "33495752"
PA_URL = "https://pubannotation.org/projects/PA-LLM-test"


def test_fetch_abstract():
    """Get PubMed abstract for a given PMID"""

    abstract = control.fetch_abstract_(DEFAULT_PMID)

    with open(os.path.join("tests", f"{DEFAULT_PMID}.txt"), "r") as f:
        abstract_ = f.read()

    assert abstract in abstract_


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


def random_word(length):
    letters = string.ascii_lowercase
    return "".join(random.choice(letters) for i in range(length))


def test_upload():
    r = requests.get(f"{PA_URL}/docs.json")
    assert DEFAULT_PMID in [doc["sourceid"] for doc in r.json()]

    with open(os.path.join("tests", f"{DEFAULT_PMID}.json"), "r") as f:
        jsonified = json.loads(f.read())

    summary_ = random_word(25)
    jsonified["blocks"][0]["obj"] = summary_

    r = control.upload_(DEFAULT_PMID, json.dumps(jsonified))

    assert summary_ in r["blocks"][0]["obj"]
