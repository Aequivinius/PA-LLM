from typing import Union

import google.generativeai as genai
from Bio import Entrez
from openai import OpenAI as openai
from streamlit import cache_data, secrets, session_state

import constants

Entrez.email = "nicola.colic@supsi.ch"
genai.configure(api_key=secrets["GEMINI_API_KEY"])
gpt_client = openai(api_key=secrets["OPENAI_API_KEY"])


def validate_pmid() -> None:
    """PMID can only be numbers"""
    pmid = session_state.pmid
    if pmid and not pmid.isdigit():
        session_state.pmid = None
        session_state.pmid_validation_error = "Only digits allowed here"


def clear_summary() -> None:
    session_state.summary = None


def session_states() -> None:
    """Initialises session variables"""
    for state, value in constants.SESSION_STATES.items():
        if state not in session_state:
            session_state[state] = value


def fetch_abstract() -> None:
    """Wraps around fetch_abstract so we can cache results"""
    abstract = fetch_abstract_(session_state.pmid)
    session_state.abstract = abstract


@cache_data
def fetch_abstract_(pmid: int) -> Union[str, None]:
    """Gets abstract for PMID from PM via Entrez"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        xml_data = Entrez.read(handle)
        article = xml_data["PubmedArticle"][0]["MedlineCitation"]["Article"]
        abstract = article["Abstract"]["AbstractText"]
        return abstract[0]
    except (IndexError, IOError) as e:
        session_state.abstract = None
        session_state.fetch_abstract_error = e
        return None


def fetch_summary():
    """Wraps around fetch_summary so we can cache results"""
    session_state.summary = fetch_summary_(session_state.abstract, session_state.llm)


@cache_data
def fetch_summary_(abstract: str, llm: str) -> str:
    """Asks LLM to give summary"""
    query = f"Summarise the following text and make it understandable to a teenager: {abstract}"

    control_module = __import__(__name__)
    query_function = getattr(control_module, constants.LLMS[llm])
    return query_function(query)


def query_bard(query: str) -> str:
    model = genai.GenerativeModel("gemini-pro")
    response = model.generate_content(query)
    return response.text


def query_chatgpt(query: str) -> str:
    messages = [
        {
            "role": "system",
            "content": "You are a helpful assistant for text summarization.",
        },
        {
            "role": "user",
            "content": query,
        },
    ]

    res = gpt_client.chat.completions.create(
        model="gpt-3.5-turbo",
        max_tokens=constants.MAX_TOKENS,
        temperature=constants.TEMPERATURE,
        messages=messages,
    )

    # finish_reason = res.choices[0].finish_reason
    return res.choices[0].message.content
