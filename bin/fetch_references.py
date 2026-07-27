#!/usr/bin/env python3
"""
fetch_references.py
PT-BR: Baixa do UniProt as sequencias de referencia curadas a partir de accessions fixadas.
EN-US: Downloads the curated reference sequences from UniProt using pinned accessions.

PT-BR: Toda sequencia vem do UniProt e tem o comprimento conferido contra o valor
       fixado na tabela. Nada e sintetizado localmente: divergencia aborta com erro.
EN-US: Every sequence comes from UniProt and has its length checked against the value
       pinned in the table. Nothing is synthesised locally: a mismatch aborts with an error.
"""

import sys
import os
import time
import argparse
import urllib.request
import urllib.error

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
USER_AGENT = "bafes-urease-mining/0.3 (https://github.com/cbafes-unb/bafes-urease)"
REQUEST_DELAY = 0.35


def read_accession_table(path):
    """
    PT-BR: Le a tabela de accessions, ignorando comentarios e o cabecalho.
    EN-US: Reads the accession table, skipping comments and the header.
    """
    entries = []
    with open(path, 'r', encoding='utf-8') as handle:
        for raw in handle:
            line = raw.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split('\t')
            if parts[0] == 'target_gene':
                continue
            if len(parts) < 5:
                raise ValueError(f"Linha malformada / Malformed line: {line!r}")
            entries.append({
                'gene': parts[0],
                'accession': parts[1],
                'expected_length': int(parts[2]),
                'status': parts[3],
                'organism': parts[4],
                'note': parts[5] if len(parts) > 5 else '',
            })
    return entries


def fetch_fasta(accession, attempts=4, timeout=30):
    """
    PT-BR: Busca o FASTA de um accession no UniProt, com retry e backoff.
    EN-US: Fetches an accession's FASTA from UniProt, with retry and backoff.
    """
    url = UNIPROT_FASTA.format(acc=accession)
    backoff = 1.0
    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode('utf-8')
        except Exception as exc:  # HTTPError, URLError, timeout
            last_error = exc
            if attempt == attempts:
                break
            print(f"    [RETRY {attempt}/{attempts - 1}] {accession}: {exc} — aguardando / waiting {backoff:.0f}s")
            time.sleep(backoff)
            backoff *= 2

    raise RuntimeError(f"Falha ao baixar / Failed to download {accession}: {last_error}")


def parse_single_fasta(text):
    """
    PT-BR: Extrai (cabecalho, sequencia) de um FASTA de registro unico.
    EN-US: Extracts (header, sequence) from a single-record FASTA.
    """
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith('>'):
        raise ValueError("Resposta do UniProt nao e um FASTA valido / UniProt response is not a valid FASTA")
    return lines[0][1:], ''.join(lines[1:])


def wrap(sequence, width=60):
    return '\n'.join(sequence[i:i + width] for i in range(0, len(sequence), width))


def main():
    parser = argparse.ArgumentParser(description="Fetch curated urease reference sequences from UniProt")
    parser.add_argument("--accessions", default="data/urease_reference_accessions.tsv",
                        help="Tabela de accessions fixadas / Pinned accession table")
    parser.add_argument("--out", default="data/urease_references.fasta",
                        help="FASTA de saida / Output FASTA")
    args = parser.parse_args()

    if not os.path.exists(args.accessions):
        print(f"[ERRO/ERROR] Tabela nao encontrada / Table not found: {args.accessions}")
        return 1

    entries = read_accession_table(args.accessions)
    print(f">>> Baixando {len(entries)} referencias do UniProt / Downloading {len(entries)} references from UniProt...")

    records = []
    failures = []

    for index, entry in enumerate(entries):
        if index > 0:
            time.sleep(REQUEST_DELAY)

        acc = entry['accession']
        try:
            header, sequence = parse_single_fasta(fetch_fasta(acc))
        except (RuntimeError, ValueError) as exc:
            print(f"  [ERRO/ERROR] {entry['gene']} ({acc}): {exc}")
            failures.append(acc)
            continue

        if len(sequence) != entry['expected_length']:
            print(f"  [ERRO/ERROR] {entry['gene']} ({acc}): comprimento {len(sequence)} != esperado "
                  f"{entry['expected_length']} / length mismatch. O registro no UniProt mudou; "
                  f"revise a tabela / the UniProt record changed, review the table.")
            failures.append(acc)
            continue

        organism = entry['organism'].replace(' ', '_')
        records.append(
            f">ref|{entry['gene']}|UniProt:{acc}|{organism}|{entry['status']}|{len(sequence)}aa\n{wrap(sequence)}"
        )
        print(f"  [OK] {entry['gene']:<5} {acc} {len(sequence):>5} aa  {entry['organism']}")

    if failures:
        print(f"\n[FALHA/FAILURE] {len(failures)} referencia(s) nao obtida(s) / reference(s) not retrieved: "
              f"{', '.join(failures)}")
        print("                O arquivo de saida NAO foi escrito / The output file was NOT written.")
        return 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(records) + '\n')

    print(f"\n[OK] {len(records)} referencias reais escritas em / real references written to: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
