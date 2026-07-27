#!/usr/bin/env python3
"""
ncbi_http.py
PT-BR: Utilitarios HTTP compartilhados: retry com backoff e rate limit do NCBI E-utilities.
EN-US: Shared HTTP utilities: retry with backoff and NCBI E-utilities rate limiting.
"""

import os
import time
import json
import ssl
import urllib.request
import urllib.error
import urllib.parse

USER_AGENT = "bafes-urease-mining/0.3 (https://github.com/cbafes-unb/bafes-urease)"

# PT-BR: Codigos HTTP transitorios que merecem nova tentativa.
# EN-US: Transient HTTP codes worth retrying.
RETRYABLE_STATUS = {403, 429, 500, 502, 503, 504}

# PT-BR: Limites do NCBI E-utilities: 3 req/s sem chave, 10 req/s com chave.
# EN-US: NCBI E-utilities limits: 3 req/s without a key, 10 req/s with one.
NCBI_DELAY_NO_KEY = 0.40
NCBI_DELAY_WITH_KEY = 0.11

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get_ssl_context():
    """
    PT-BR: Contexto SSL com verificacao normal; relaxa so se BAFES_INSECURE_SSL=1.
    EN-US: SSL context with normal verification; relaxed only if BAFES_INSECURE_SSL=1.
    """
    ctx = ssl.create_default_context()
    if os.environ.get("BAFES_INSECURE_SSL") == "1":
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    return ctx


def http_request(url, method="GET", attempts=5, timeout=30, read_body=True, quiet=False):
    """
    PT-BR: Requisicao HTTP com retry e backoff exponencial, honrando Retry-After.
    EN-US: HTTP request with retry and exponential backoff, honouring Retry-After.

    PT-BR: Retorna (status, corpo, erro). Erro None indica sucesso.
    EN-US: Returns (status, body, error). A None error means success.
    """
    ctx = get_ssl_context()
    last_error = None
    backoff = 1.0

    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method=method)
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
                body = resp.read().decode("utf-8", errors="replace") if read_body else ""
                return resp.status, body, None
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code} {exc.reason}"
            if exc.code not in RETRYABLE_STATUS or attempt == attempts:
                return exc.code, "", last_error
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            wait = backoff
            if retry_after:
                try:
                    wait = max(backoff, float(retry_after))
                except ValueError:
                    pass
        except Exception as exc:  # URLError, timeout, DNS, TLS
            last_error = str(exc)
            if attempt == attempts:
                return None, "", last_error
            wait = backoff

        if not quiet:
            print(f"      [RETRY {attempt}/{attempts - 1}] {last_error} — aguardando / waiting {wait:.0f}s...")
        time.sleep(wait)
        backoff *= 2

    return None, "", last_error


def ncbi_api_key():
    return os.environ.get("NCBI_API_KEY", "").strip()


def ncbi_delay():
    """
    PT-BR: Intervalo minimo entre requisicoes ao E-utilities.
    EN-US: Minimum interval between E-utilities requests.
    """
    return NCBI_DELAY_WITH_KEY if ncbi_api_key() else NCBI_DELAY_NO_KEY


def eutils_url(endpoint, **params):
    """
    PT-BR: Monta uma URL do E-utilities com tool/email/api_key exigidos pela politica do NCBI.
    EN-US: Builds an E-utilities URL with the tool/email/api_key required by NCBI policy.
    """
    params.setdefault("tool", "bafes-urease-mining")
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if email:
        params["email"] = email
    key = ncbi_api_key()
    if key:
        params["api_key"] = key
    return f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"


def eutils_json(endpoint, **params):
    """
    PT-BR: Chama o E-utilities e devolve (dados, erro). Dados None em caso de falha.
    EN-US: Calls E-utilities and returns (data, error). Data is None on failure.
    """
    params.setdefault("retmode", "json")
    _, body, error = http_request(eutils_url(endpoint, **params))
    if error is not None:
        return None, error
    try:
        return json.loads(body), None
    except json.JSONDecodeError as exc:
        return None, f"resposta invalida do NCBI / invalid NCBI response ({exc})"
