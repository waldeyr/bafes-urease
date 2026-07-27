#!/usr/bin/env python3
"""
download_genome.py
PT-BR: Baixa o assembly genomico real correspondente a um accession GenBank.
EN-US: Downloads the real genome assembly matching a GenBank accession.

PT-BR: Os accessions do estudo sao WGS master (ex.: VKHW00000000.1). O registro WGS
       master NAO contem sequencia — um efetch nele devolve zero bytes. Por isso o
       accession e primeiro resolvido para o assembly (GCA_/GCF_) e so entao baixado.
EN-US: The study's accessions are WGS master records (e.g. VKHW00000000.1). A WGS
       master record carries NO sequence — an efetch on it returns zero bytes. The
       accession is therefore resolved to its assembly (GCA_/GCF_) before download.

PT-BR: Este script NUNCA fabrica sequencia. Se o genoma real nao for obtido, sai com erro.
EN-US: This script NEVER fabricates sequence. If the real genome cannot be obtained, it exits with an error.
"""

import sys
import os
import time
import glob
import gzip
import shutil
import argparse
import subprocess

from ncbi_http import http_request, eutils_json, ncbi_delay

# PT-BR: Um assembly bacteriano legitimo nao tem menos que isso.
# EN-US: A legitimate bacterial assembly is not smaller than this.
MIN_GENOME_BASES = 500_000
VALID_BASES = set("ACGTURYKMSWBDHVN")


def resolve_assembly(accession):
    """
    PT-BR: Resolve um accession (WGS master ou GCA/GCF) para o accession do assembly.
    EN-US: Resolves an accession (WGS master or GCA/GCF) to its assembly accession.

    PT-BR: Retorna (accession_assembly, ftp_path, erro).
    EN-US: Returns (assembly_accession, ftp_path, error).
    """
    if accession.upper().startswith(("GCA_", "GCF_")):
        term = accession
    else:
        # PT-BR: A busca no db=assembly usa o prefixo WGS sem a versao.
        # EN-US: The db=assembly search uses the WGS prefix without the version.
        term = accession.split('.')[0]

    data, error = eutils_json("esearch.fcgi", db="assembly", term=term)
    if error is not None:
        return None, None, f"esearch falhou / failed: {error}"

    idlist = data.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return None, None, f"nenhum assembly ligado ao accession / no assembly linked to accession {accession}"

    time.sleep(ncbi_delay())
    data, error = eutils_json("esummary.fcgi", db="assembly", id=idlist[0])
    if error is not None:
        return None, None, f"esummary falhou / failed: {error}"

    result = data.get("result", {})
    uids = result.get("uids", [])
    if not uids:
        return None, None, "esummary sem resultados / returned no results"

    record = result[uids[0]]
    asm_acc = record.get("assemblyaccession", "")
    # PT-BR: Prefere RefSeq; cai para GenBank quando ausente.
    # EN-US: Prefers RefSeq; falls back to GenBank when absent.
    ftp_path = record.get("ftppath_refseq") or record.get("ftppath_genbank") or ""

    if not asm_acc:
        return None, None, "esummary sem assemblyaccession / missing assemblyaccession"

    return asm_acc, ftp_path, None


def download_via_datasets(asm_acc, workdir):
    """
    PT-BR: Baixa o assembly com a CLI ncbi-datasets. Retorna o caminho do .fna ou None.
    EN-US: Downloads the assembly with the ncbi-datasets CLI. Returns the .fna path or None.
    """
    if shutil.which("datasets") is None:
        print("  [INFO] CLI 'datasets' indisponivel / unavailable; usando FTP.")
        return None

    zip_path = os.path.join(workdir, "assembly.zip")
    cmd = ["datasets", "download", "genome", "accession", asm_acc,
           "--include", "genome", "--filename", zip_path]
    print(f"  [1/2] datasets download genome accession {asm_acc}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        print("  [AVISO/WARNING] 'datasets' excedeu o tempo limite / timed out.")
        return None

    if proc.returncode != 0 or not os.path.exists(zip_path):
        print(f"  [AVISO/WARNING] 'datasets' falhou / failed: {proc.stderr.strip()[:300]}")
        return None

    extract_dir = os.path.join(workdir, "datasets_out")
    try:
        subprocess.run(["unzip", "-o", "-q", zip_path, "-d", extract_dir], check=True, timeout=600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        print(f"  [AVISO/WARNING] unzip falhou / failed: {exc}")
        return None

    matches = glob.glob(os.path.join(extract_dir, "**", "*.fna"), recursive=True)
    if not matches:
        print("  [AVISO/WARNING] Nenhum .fna dentro do pacote / No .fna inside the package.")
        return None
    return matches[0]


def download_via_ftp(ftp_path, workdir):
    """
    PT-BR: Baixa o *_genomic.fna.gz do caminho FTP do assembly no NCBI.
    EN-US: Downloads *_genomic.fna.gz from the assembly's NCBI FTP path.
    """
    if not ftp_path:
        return None

    basename = ftp_path.rstrip('/').split('/')[-1]
    url = f"{ftp_path.replace('ftp://', 'https://')}/{basename}_genomic.fna.gz"
    print(f"  [2/2] FTP fallback: {url}")

    status, body, error = http_request(url, timeout=1800)
    if error is not None:
        print(f"  [AVISO/WARNING] Download FTP falhou / FTP download failed: {error}")
        return None

    gz_path = os.path.join(workdir, "genomic.fna.gz")
    # PT-BR: http_request decodifica texto; reabrimos em binario para o gzip.
    # EN-US: http_request decodes text; re-fetch in binary for gzip.
    import urllib.request
    from ncbi_http import USER_AGENT, get_ssl_context
    try:
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=1800, context=get_ssl_context()) as resp, \
                open(gz_path, 'wb') as out:
            shutil.copyfileobj(resp, out)
    except Exception as exc:
        print(f"  [AVISO/WARNING] Download binario falhou / binary download failed: {exc}")
        return None

    fna_path = os.path.join(workdir, "genomic.fna")
    try:
        with gzip.open(gz_path, 'rb') as src, open(fna_path, 'wb') as dst:
            shutil.copyfileobj(src, dst)
    except OSError as exc:
        print(f"  [AVISO/WARNING] Descompactacao falhou / decompression failed: {exc}")
        return None

    return fna_path


def validate_fasta(path):
    """
    PT-BR: Confere que o arquivo e um FASTA nucleotidico plausivel.
    EN-US: Checks that the file is a plausible nucleotide FASTA.

    PT-BR: Retorna (n_contigs, n_bases, erro).
    EN-US: Returns (n_contigs, n_bases, error).
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0, 0, "arquivo vazio ou inexistente / file empty or missing"

    n_contigs = 0
    n_bases = 0
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        first = handle.readline()
        if not first.startswith('>'):
            return 0, 0, "nao comeca com cabecalho FASTA / does not start with a FASTA header"
        n_contigs = 1
        for line in handle:
            if line.startswith('>'):
                n_contigs += 1
            else:
                n_bases += len(line.strip())

    if n_bases < MIN_GENOME_BASES:
        return n_contigs, n_bases, (f"apenas {n_bases} bases (< {MIN_GENOME_BASES}) / only {n_bases} bases; "
                                    "nao e um assembly bacteriano completo / not a complete bacterial assembly")

    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        sample = ''.join(line.strip().upper() for line in handle
                         if not line.startswith('>'))[:10000]
    invalid = set(sample) - VALID_BASES
    if invalid:
        return n_contigs, n_bases, f"caracteres invalidos na sequencia / invalid sequence characters: {sorted(invalid)[:5]}"

    return n_contigs, n_bases, None


def main():
    parser = argparse.ArgumentParser(description="Download a real genome assembly from NCBI")
    parser.add_argument("--accession", required=True, help="Accession GenBank (WGS master ou GCA/GCF)")
    parser.add_argument("--strain", required=True, help="Identificador da estirpe / Strain ID")
    parser.add_argument("--out", required=True, help="Caminho do FASTA de saida / Output FASTA path")
    args = parser.parse_args()

    print(f">>> Baixando genoma de / Downloading genome for {args.strain} ({args.accession})")

    asm_acc, ftp_path, error = resolve_assembly(args.accession)
    if error is not None:
        print(f"[ERRO/ERROR] Nao foi possivel resolver o assembly / Could not resolve assembly "
              f"for {args.accession}: {error}")
        return 1
    print(f"  [OK] {args.accession} -> assembly {asm_acc}")

    workdir = f".download_{args.strain}"
    os.makedirs(workdir, exist_ok=True)

    fna = download_via_datasets(asm_acc, workdir)
    if fna is None:
        fna = download_via_ftp(ftp_path, workdir)

    if fna is None:
        print(f"[ERRO/ERROR] Falha ao baixar o genoma de {args.strain} ({asm_acc}). "
              f"/ Failed to download the genome for {args.strain} ({asm_acc}).")
        print("             Nenhuma sequencia sera fabricada / No sequence will be fabricated.")
        return 1

    shutil.copyfile(fna, args.out)
    shutil.rmtree(workdir, ignore_errors=True)

    n_contigs, n_bases, error = validate_fasta(args.out)
    if error is not None:
        print(f"[ERRO/ERROR] Genoma baixado e invalido / Downloaded genome is invalid: {error}")
        os.remove(args.out)
        return 1

    print(f"[OK] {args.strain}: {n_contigs} contigs, {n_bases:,} bases -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
