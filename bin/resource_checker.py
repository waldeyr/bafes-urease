#!/usr/bin/env python3
"""
resource_checker.py
PT-BR: Validação prévia de recursos de rede e insumos de entrada (dry-run / pre-flight checks).
EN-US: Pre-flight check and resource connectivity validation script (dry-run / pre-flight checks).

PT-BR: Contrato de saída — 0 = tudo ok; 4 = recurso crítico indisponível após todas as tentativas.
EN-US: Exit contract — 0 = all good; 4 = critical resource unavailable after all retries.
"""

import sys
import os
import time
import shutil
import argparse
import json

from ncbi_http import (http_request, eutils_url, ncbi_api_key,
                       NCBI_DELAY_NO_KEY, NCBI_DELAY_WITH_KEY)


def check_accessions(accessions_file):
    """
    PT-BR: Valida se cada accession do GenBank existe no NCBI, respeitando o rate limit.
    EN-US: Validates that each GenBank accession exists on NCBI, respecting the rate limit.
    """
    print(">>> [PRE-CHECK 1/5] Verificando accessions GenBank no NCBI / Checking GenBank accessions on NCBI...")
    if not os.path.exists(accessions_file):
        print(f"  [ERRO/ERROR] Arquivo de accessions não encontrado / Accessions file not found: {accessions_file}")
        return False

    accessions = []
    with open(accessions_file, 'r', encoding='utf-8') as f:
        f.readline()  # PT-BR: cabeçalho / EN-US: header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[2]:
                accessions.append((parts[0], parts[1], parts[2]))
            elif len(parts) >= 1 and parts[0]:
                accessions.append((parts[0], "Unknown", parts[0]))

    if not accessions:
        print("  [ERRO/ERROR] Nenhum accession encontrado / No accessions found.")
        return False

    api_key = ncbi_api_key()
    delay = NCBI_DELAY_WITH_KEY if api_key else NCBI_DELAY_NO_KEY
    if api_key:
        print("  [INFO] Usando NCBI_API_KEY (10 req/s) / Using NCBI_API_KEY (10 req/s).")
    else:
        print("  [INFO] Sem NCBI_API_KEY: limitando a 3 req/s / No NCBI_API_KEY: throttling to 3 req/s.")

    missing = []      # PT-BR: accession inexistente / EN-US: accession does not exist
    unreachable = []  # PT-BR: NCBI fora do ar / EN-US: NCBI unreachable

    for index, (strain, species, acc) in enumerate(accessions):
        if index > 0:
            time.sleep(delay)

        found, error = _ncbi_accession_exists(acc, delay)
        if error is not None:
            print(f"  [ERRO/ERROR] NCBI inacessível para / NCBI unreachable for {acc}: {error}")
            unreachable.append(acc)
        elif found:
            print(f"  [OK] Strain {strain} ({species}): Accession {acc} verificado no NCBI / Verified on NCBI.")
        else:
            print(f"  [ERRO/ERROR] Strain {strain}: Accession {acc} não encontrado no NCBI / Not found on NCBI.")
            missing.append(acc)

    if unreachable:
        print(f"  [FALHA/FAILURE] {len(unreachable)} accession(s) sem resposta do NCBI após todas as tentativas.")
        print(f"                  {len(unreachable)} accession(s) got no NCBI response after all retries.")
        if not api_key:
            print("                  Dica / Hint: exporte NCBI_API_KEY para elevar o limite de requisições.")
            print("                  Export NCBI_API_KEY to raise the request rate limit.")
    if missing:
        print(f"  [FALHA/FAILURE] {len(missing)} accession(s) inexistente(s) no NCBI / not present on NCBI: {', '.join(missing)}")

    return not unreachable and not missing


def _ncbi_accession_exists(acc, delay):
    """
    PT-BR: Consulta esearch por [ACCN] e, se der 0, refaz com o termo cru (accessions WGS master).
    EN-US: Queries esearch by [ACCN] and, on 0 hits, retries with the raw term (WGS master accessions).

    PT-BR: Retorna (encontrado, erro) / EN-US: Returns (found, error)
    """
    bare = acc.split('.')[0]
    terms = [f"{bare}[ACCN]", acc]

    for position, term in enumerate(terms):
        if position > 0:
            time.sleep(delay)

        _, body, error = http_request(eutils_url("esearch.fcgi", db="nuccore",
                                                 term=term, retmode="json"))
        if error is not None:
            return False, error

        try:
            count = int(json.loads(body).get("esearchresult", {}).get("count", "0"))
        except (ValueError, json.JSONDecodeError) as exc:
            return False, f"resposta inválida do NCBI / invalid NCBI response ({exc})"

        if count > 0:
            return True, None

    return False, None


DB_VERSION_DIR = "db/.versions"


def recorded_version(key):
    """
    PT-BR: Lê a versão que o run.sh registrou ao baixar o banco. Ver 'db_version_set'.
    EN-US: Reads the version run.sh recorded when it downloaded the database. See 'db_version_set'.
    """
    path = os.path.join(DB_VERSION_DIR, key)
    try:
        with open(path, 'r', encoding='utf-8') as handle:
            return handle.readline().strip() or None
    except OSError:
        return None


def check_pfam_url(pfam_path=None):
    """
    PT-BR: Confirma o Pfam local ou testa o repositório do EMBL-EBI (não bloqueante).
    EN-US: Confirms the local Pfam file or probes the EMBL-EBI repository (non-blocking).
    """
    print("\n>>> [PRE-CHECK 2/5] Verificando Pfam HMM (EMBL-EBI FTP)...")
    if pfam_path and os.path.exists(pfam_path) and os.path.getsize(pfam_path) > 0:
        version = recorded_version("pfam")
        suffix = f" (release {version})" if version else ""
        print(f"  [OK] Pfam local encontrado em / Local Pfam found at: {pfam_path}{suffix}")
        if not version:
            print("  [INFO] Release não registrada / Release not recorded: baixado antes do registro de versões.")
        return True

    url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    status, _, error = http_request(url, method="HEAD", attempts=3, read_body=False)
    if error is None:
        print(f"  [OK] Repositório Pfam-A online no EMBL-EBI (Status {status}).")
    else:
        print(f"  [AVISO/WARNING] Verificação online do Pfam retornou aviso / Pfam online check notice: {error}")
    return True


def check_bakta_db(bakta_db_path):
    """
    PT-BR: Verifica o banco Bakta local ou a disponibilidade do registro no Zenodo.
    EN-US: Verifies the local Bakta database or the availability of the Zenodo record.
    """
    print("\n>>> [PRE-CHECK 3/5] Verificando banco do Bakta / Checking Bakta database...")
    version_json = os.path.join(bakta_db_path, "version.json") if bakta_db_path else None
    if version_json and os.path.exists(version_json):
        label = _bakta_version_label(version_json)
        print(f"  [OK] Banco local do Bakta encontrado em / Local Bakta DB found at: {bakta_db_path}{label}")
        return True

    if bakta_db_path and os.path.isdir(bakta_db_path):
        print(f"  [AVISO/WARNING] Diretório existe mas sem version.json / Directory exists but has no version.json: {bakta_db_path}")

    print("  [INFO] Testando acessibilidade ao Zenodo Bakta DB / Testing Zenodo Bakta DB accessibility...")
    # PT-BR: API do Zenodo; o resolvedor DOI é servido por CDN que bloqueia clientes não-browser (403).
    # EN-US: Zenodo API; the DOI resolver is CDN-fronted and blocks non-browser clients (403).
    url = "https://zenodo.org/api/records/4247252"
    status, _, error = http_request(url, method="HEAD", read_body=False)
    if error is not None:
        status, _, error = http_request(url, method="GET")

    if error is None:
        print(f"  [OK] Repositório Zenodo Bakta DB online (Status {status}).")
        return True

    print(f"  [ERRO/ERROR] Falha ao acessar Bakta DB Zenodo / Failed to access Bakta DB Zenodo: {error}")
    print("               Baixe o banco com / Download the database with: ./run.sh --bootstrap")
    return False


def _bakta_version_label(version_json):
    """
    PT-BR: Monta ' (vMAJOR.MINOR/tipo)' a partir do version.json do bakta e do tipo
           registrado pelo run.sh. Devolve '' se o arquivo não for legível.
    EN-US: Builds ' (vMAJOR.MINOR/type)' from bakta's version.json plus the type recorded
           by run.sh. Returns '' when the file is not readable.
    """
    try:
        with open(version_json, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return ""

    major, minor = data.get("major"), data.get("minor")
    if major is None or minor is None:
        return ""

    db_type = recorded_version("bakta_type")
    return f" (v{major}.{minor}/{db_type})" if db_type else f" (v{major}.{minor})"


def check_references(ref_fasta):
    """
    PT-BR: Valida existência e conteúdo do FASTA de referências curadas.
    EN-US: Validates presence and content of the curated references FASTA file.
    """
    print("\n>>> [PRE-CHECK 4/5] Verificando referências proteicas / Checking protein references...")
    if not (os.path.exists(ref_fasta) and os.path.getsize(ref_fasta) > 0):
        print(f"  [ERRO/ERROR] Arquivo de referências ausente / Reference file missing: {ref_fasta}")
        return False

    with open(ref_fasta, 'r', encoding='utf-8') as f:
        n_seqs = sum(1 for line in f if line.startswith('>'))

    if n_seqs == 0:
        print(f"  [ERRO/ERROR] Nenhuma sequência FASTA em / No FASTA sequence in: {ref_fasta}")
        return False

    print(f"  [OK] Arquivo de referências curadas válido / Valid reference file: {ref_fasta} ({n_seqs} seqs)")
    return True


def check_runtime_tools(in_container=False):
    """
    PT-BR: Checa as ferramentas que o pipeline realmente invoca.
    EN-US: Checks the tools the pipeline actually invokes.
    """
    print("\n>>> [PRE-CHECK 5/5] Verificando ferramentas no ambiente / Checking system tools...")
    tools = ['python3', 'nextflow', 'bakta', 'blastp', 'makeblastdb', 'hmmscan',
             'mafft', 'trimal', 'iqtree', 'quast.py']

    missing = [tool for tool in tools if shutil.which(tool) is None]
    for tool in tools:
        if tool not in missing:
            print(f"  [OK] Ferramenta / Tool '{tool}' ok.")

    if not missing:
        return True

    where = "no container" if in_container else "no host PATH"
    where_en = "in the container" if in_container else "on the host PATH"
    for tool in missing:
        level = "ERRO/ERROR" if in_container else "AVISO/WARNING"
        print(f"  [{level}] Ferramenta / Tool '{tool}' não localizada {where} / not found {where_en}.")

    if in_container:
        print("  [FALHA/FAILURE] A imagem do container está incompleta / The container image is incomplete.")
        print("                  Reconstrua com / Rebuild with: ./run.sh --bootstrap --force-rebuild")
        return False

    print("  [INFO] Execução containerizada não depende do PATH do host.")
    print("  [INFO] Containerized execution does not depend on the host PATH.")
    return True


def main():
    parser = argparse.ArgumentParser(description="Resource Checker & Pre-flight Validator (PT-BR / EN-US)")
    parser.add_argument("--accessions", default="data/accessions.tsv", help="Caminho do accessions.tsv / Path to accessions.tsv")
    parser.add_argument("--references", default="data/urease_references.fasta", help="Caminho das referencias FASTA / Path to references FASTA")
    parser.add_argument("--bakta-db", default="db/db", help="Caminho do banco Bakta / Path to Bakta DB")
    parser.add_argument("--pfam", default="Pfam-A.hmm", help="Caminho do Pfam-A.hmm / Path to Pfam-A.hmm")
    parser.add_argument("--in-container", action="store_true",
                        help="Indica execução dentro do container / Signals execution inside the container")
    args = parser.parse_args()

    print("==================================================")
    print("   PRE-FLIGHT RESOURCE & CONNECTIVITY CHECK       ")
    print("==================================================")

    ok1 = check_accessions(args.accessions)
    ok2 = check_pfam_url(args.pfam)
    ok3 = check_bakta_db(args.bakta_db)
    ok4 = check_references(args.references)
    ok5 = check_runtime_tools(args.in_container)

    print("==================================================")
    if ok1 and ok2 and ok3 and ok4 and ok5:
        print("[SUCESSO/SUCCESS] Todos os pre-checks passaram! / All pre-checks passed!")
        sys.exit(0)

    failed = [name for name, ok in (("1 (NCBI accessions)", ok1), ("2 (Pfam)", ok2),
                                    ("3 (Bakta DB)", ok3), ("4 (references)", ok4),
                                    ("5 (tools)", ok5)) if not ok]
    print(f"[FALHA/FAILURE] Falha nos pre-checks / Failed pre-checks: {', '.join(failed)}")
    sys.exit(4)


if __name__ == "__main__":
    main()
