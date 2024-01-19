ARTICLE_OPTIONS = {"abstract": "Abstract only", "article": "Full article"}

LLMS_ = ["BARD", "ChatGPT"]
LLMS = {llm: f"query_{llm.lower()}" for llm in LLMS_}

"""token limit for gpt queries"""
MAX_TOKENS = 500

PA_URL = "https://pubannotation.org/projects/PA-LLM"

"""session variables and their defaults"""
SESSION_STATES = {
    "abstract": None,
    "fetch_abstract_error": None,
    "json": "",
    "llm": "BARD",
    "pmid": None,
    "pmid_validation_error": None,
    "summary": None,
    "upload": None,
}

"""temperature default for gpt queries"""
TEMPERATURE = 0.2
