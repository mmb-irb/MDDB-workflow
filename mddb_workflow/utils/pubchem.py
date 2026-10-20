import time
import requests
from mddb_workflow.utils.auxiliar import RequestError, RateLimitError, retry_request, get_retry_after
from mddb_workflow.utils.type_hints import *
from functools import lru_cache


def _make_get_request(url: str) -> requests.Response:
    return requests.get(url, timeout=30)


def _make_post_request(url, payload):
    return requests.post(url, data=payload, timeout=30)


@retry_request
def make_pubchem_request(
    url: str,
    payload: Optional[dict] = None,
    error_context: str = 'PubChem',
    no_result_codes: tuple[int, ...] = (),
) -> Optional[dict]:
    """Make a request to PubChem API, with retries and error handling.

    Args:
        url (str): The URL to make the request to.
        payload (Optional[dict]): The payload to send in a POST request. If None, a GET request will be made.
        error_context (str): Context for error messages.
        no_result_codes (tuple[int, ...]): Status codes which mean PubChem has no result for this request.
            In these cases return None. Any other error status raises an error.

    """
    if payload is None:
        r = _make_get_request(url)
    else:
        r = _make_post_request(url, payload)
    # Add a small sleep if we are in the green zone to avoid hitting the rate limit
    # https://pubchem.ncbi.nlm.nih.gov/docs/dynamic-request-throttling
    throttling = r.headers.get('X-Throttling-Control', 'No throttling').split(',')[0]

    if 'Yellow' in throttling:
        time.sleep(0.2)
    elif 'Red' in throttling:
        time.sleep(0.5)
    elif 'Black' in throttling:
        time.sleep(1)
    if r.status_code == 200:
        return r.json()
    if r.status_code in no_result_codes:
        return None
    if r.status_code == 429:
        raise RateLimitError(f"Too many requests in the {error_context} request. Throttling: {throttling}",
            retry_after=get_retry_after(r.headers))
    if 'temporarily unable' in r.text:
        raise requests.exceptions.ConnectionError("PubChem is temporarily unable to process requests. Please try again later.")
    raise RequestError(f"Something went wrong with the {error_context} request (error " + str(r.status_code) + "): " + r.text + "\nThrottling: " + throttling)


@lru_cache(maxsize=None)
def pubchem_standardization(inchi: str) -> Optional[list[dict]]:
    """Try PubChem standardization service to write things like the bonds from aromatic rings in the same order.
    You can see examples on the paper https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0293-8.
    """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/standardize/inchi/JSON"
    # PubChem answers 400 when it is not able to standardize the structure
    data = make_pubchem_request(url, payload={'inchi': inchi},
        error_context='PubChem standardization', no_result_codes=(400,))
    if data:
        results = []
        for compound in data['PC_Compounds']:
            if compound['id']['id']['cid']:
                cid = str(compound['id']['id']['cid'])
            else:
                continue
            standard_inchi = None
            for prop in compound['props']:
                if prop['urn']['label'] == 'InChI':
                    standard_inchi = prop['value']['sval']
            results.append({'pubchem': cid, 'inchi': standard_inchi})
        return results
    return None
