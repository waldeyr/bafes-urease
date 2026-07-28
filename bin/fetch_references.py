#!/usr/bin/env python3
"""
fetch_references.py
PT-BR: Baixa do UniProt as sequencias de referencia de urease definidas por uma consulta
       REST: proteinas de Firmicutes (taxonomy_id:1239) revisadas (Swiss-Prot) cujo nome
       de proteina (protein_name) e urease.
EN-US: Downloads the urease reference sequences from UniProt as defined by a REST query:
       Firmicutes proteins (taxonomy_id:1239), reviewed (Swiss-Prot), whose protein_name
       is urease.

PT-BR: Toda sequencia vem do endpoint /stream do UniProt. Nada e sintetizado localmente:
       falha de rede, resposta truncada ou conjunto vazio abortam sem escrever o FASTA.
EN-US: Every sequence comes from UniProt's /stream endpoint. Nothing is synthesised locally:
       a network failure, a truncated response or an empty result set aborts without
       writing the FASTA.
"""

import sys
import os
import gzip
import time
import hashlib
import argparse
import urllib.parse
import urllib.request
import urllib.error

# PT-BR: Consulta oficial das referencias de urease deste estudo.
# EN-US: Official urease reference query for this study.
UNIPROT_QUERY = "(protein_name:urease) AND (taxonomy_id:1239) AND (reviewed:true)"
UNIPROT_STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?compressed=true&format=fasta"
    "&query=%28protein_name%3Aurease%29+AND+%28taxonomy_id%3A1239%29+AND+%28reviewed%3Atrue%29"
)
USER_AGENT = "bafes-urease-mining/0.4 (https://github.com/cbafes-unb/bafes-urease)"

# PT-BR: Piso de sanidade — a consulta devolve milhares de registros; um numero muito
#        baixo indica resposta truncada ou consulta alterada, e nao um resultado real.
# EN-US: Sanity floor — the query returns thousands of records; a very low number means a
#        truncated response or a changed query, not a real result.
MIN_SEQUENCES = 100


def build_url(query):
    """
    PT-BR: Monta a URL do endpoint /stream para uma consulta arbitraria.
    EN-US: Builds the /stream endpoint URL for an arbitrary query.
    """
    params = urllib.parse.urlencode({
        'compressed': 'true',
        'format': 'fasta',
        'query': query,
    })
    return f"https://rest.uniprot.org/uniprotkb/stream?{params}"


def fetch_stream(url, attempts=4, timeout=300):
    """
    PT-BR: Baixa o FASTA completo da consulta, com retry, backoff e gunzip.
    EN-US: Downloads the full FASTA for the query, with retry, backoff and gunzip.

    PT-BR: Devolve (texto, cabecalhos_http). O gzip so descomprime por inteiro se o corpo
           chegou completo — resposta cortada no meio levanta erro em vez de virar FASTA parcial.
    EN-US: Returns (text, http_headers). gzip only inflates if the body arrived complete —
           a cut-off response raises instead of becoming a partial FASTA.
    """
    backoff = 2.0
    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = resp.read()
                headers = dict(resp.headers.items())
            if payload[:2] == b'\x1f\x8b':
                payload = gzip.decompress(payload)
            return payload.decode('utf-8'), headers
        except Exception as exc:  # HTTPError, URLError, timeout, BadGzipFile, EOFError
            last_error = exc
            if attempt == attempts:
                break
            print(f"    [RETRY {attempt}/{attempts - 1}] {exc} — aguardando / waiting {backoff:.0f}s")
            time.sleep(backoff)
            backoff *= 2

    raise RuntimeError(f"Falha ao baixar do UniProt / Failed to download from UniProt: {last_error}")


def parse_uniprot_header(header):
    """
    PT-BR: Extrai (accession, gene, organismo, status) de um cabecalho FASTA do UniProt,
           no formato db|accession|entry Nome OS=... OX=... GN=... PE=... SV=...
    EN-US: Extracts (accession, gene, organism, status) from a UniProt FASTA header in the
           db|accession|entry Name OS=... OX=... GN=... PE=... SV=... format.
    """
    ident, _, description = header.partition(' ')
    parts = ident.split('|')
    if len(parts) != 3:
        raise ValueError(f"Cabecalho UniProt inesperado / Unexpected UniProt header: {header!r}")

    db, accession, entry_name = parts
    status = 'reviewed' if db == 'sp' else 'unreviewed'

    fields = {}
    key = None
    for token in description.split():
        if len(token) > 3 and token[2] == '=' and token[:2].isalpha() and token[:2].isupper():
            key = token[:2]
            fields[key] = token[3:]
        elif key is not None:
            fields[key] += '_' + token

    gene = fields.get('GN') or entry_name
    organism = fields.get('OS') or 'unknown_organism'
    return accession, gene, organism, status


def wrap(sequence, width=60):
    return '\n'.join(sequence[i:i + width] for i in range(0, len(sequence), width))


def parse_fasta(text):
    """
    PT-BR: Percorre o FASTA multi-registro devolvendo pares (cabecalho, sequencia).
    EN-US: Walks the multi-record FASTA yielding (header, sequence) pairs.
    """
    header = None
    chunks = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(chunks)
            header = line[1:]
            chunks = []
        else:
            if header is None:
                raise ValueError("FASTA invalido: sequencia antes do cabecalho / "
                                 "Invalid FASTA: sequence before header")
            chunks.append(line)
    if header is not None:
        yield header, ''.join(chunks)


def write_provenance(path, url, headers, n_seqs, fasta_path, digest):
    """
    PT-BR: Registra a procedencia do download. A consulta e viva: o conjunto muda a cada
           release do UniProt, entao o que garante reprodutibilidade e este registro
           (release, data, contagem e hash) junto do FASTA versionado.
    EN-US: Records the download provenance. The query is live: the set changes with every
           UniProt release, so reproducibility rests on this record (release, date, count
           and hash) alongside the versioned FASTA.
    """
    lines = [
        "# PT-BR: Procedencia das referencias de urease / EN-US: Urease reference provenance",
        f"query_url\t{url}",
        f"uniprot_release\t{headers.get('X-UniProt-Release', 'unknown')}",
        f"uniprot_release_date\t{headers.get('X-UniProt-Release-Date', 'unknown')}",
        f"downloaded_at\t{headers.get('Date', 'unknown')}",
        f"fasta\t{os.path.basename(fasta_path)}",
        f"sequences\t{n_seqs}",
        f"sha256\t{digest}",
    ]
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Fetch urease reference sequences from UniProt (reviewed, Firmicutes)")
    parser.add_argument("--url", default=None,
                        help="URL completa do endpoint /stream / Full /stream endpoint URL")
    parser.add_argument("--query", default=None,
                        help=f"Consulta UniProt / UniProt query (default: {UNIPROT_QUERY!r})")
    parser.add_argument("--out", default="data/urease_references.fasta",
                        help="FASTA de saida / Output FASTA")
    parser.add_argument("--provenance", default=None,
                        help="Arquivo de procedencia / Provenance file "
                             "(default: <out>.provenance.txt)")
    parser.add_argument("--min-sequences", type=int, default=MIN_SEQUENCES,
                        help="Minimo aceitavel de sequencias / Minimum acceptable sequence count")
    args = parser.parse_args()

    if args.url and args.query:
        print("[ERRO/ERROR] Use --url ou --query, nao os dois / Use --url or --query, not both.")
        return 2

    if args.url:
        url = args.url
    elif args.query:
        url = build_url(args.query)
    else:
        url = UNIPROT_STREAM_URL

    print(">>> Baixando referencias de urease do UniProt / Downloading urease references from UniProt...")
    print(f"    {url}")

    try:
        text, headers = fetch_stream(url)
    except RuntimeError as exc:
        print(f"[ERRO/ERROR] {exc}")
        print("             Nenhuma sequencia sera inventada / No sequence will be invented.")
        return 1

    records = []
    seen = set()
    unreviewed = 0

    try:
        for header, sequence in parse_fasta(text):
            accession, gene, organism, status = parse_uniprot_header(header)
            if not sequence:
                raise ValueError(f"Sequencia vazia / Empty sequence: {accession}")
            if accession in seen:
                continue
            seen.add(accession)
            if status != 'reviewed':
                unreviewed += 1
            organism = '_'.join(organism.split())
            records.append(
                f">ref|{gene}|UniProt:{accession}|{organism}|{status}|{len(sequence)}aa\n{wrap(sequence)}"
            )
    except ValueError as exc:
        print(f"[ERRO/ERROR] Resposta do UniProt invalida / Invalid UniProt response: {exc}")
        print("             O arquivo de saida NAO foi escrito / The output file was NOT written.")
        return 1

    if len(records) < args.min_sequences:
        print(f"[ERRO/ERROR] Apenas {len(records)} sequencia(s) obtida(s), minimo {args.min_sequences} "
              f"/ only {len(records)} sequence(s) retrieved, minimum {args.min_sequences}.")
        print("             Resposta provavelmente truncada / Response likely truncated.")
        print("             O arquivo de saida NAO foi escrito / The output file was NOT written.")
        return 1

    if unreviewed:
        print(f"  [AVISO/WARNING] {unreviewed} registro(s) nao revisado(s) na resposta "
              f"/ unreviewed record(s) in the response.")

    body = '\n'.join(records) + '\n'
    digest = hashlib.sha256(body.encode('utf-8')).hexdigest()

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    with open(args.out, 'w', encoding='utf-8') as handle:
        handle.write(body)

    provenance = args.provenance or (os.path.splitext(args.out)[0] + '.provenance.txt')
    write_provenance(provenance, url, headers, len(records), args.out, digest)

    print(f"\n[OK] {len(records)} referencias reais escritas em / real references written to: {args.out}")
    print(f"[OK] Procedencia registrada em / Provenance recorded in: {provenance}")
    print(f"     UniProt release {headers.get('X-UniProt-Release', 'unknown')} | sha256 {digest[:16]}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
