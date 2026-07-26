#!/usr/bin/env python3
"""
resource_checker.py
PT-BR: Validação prévia de recursos de rede e insumos de entrada (dry-run / pre-flight checks).
EN-US: Pre-flight check and resource connectivity validation script (dry-run / pre-flight checks).
"""

import sys
import os
import argparse
import urllib.request
import json
import ssl

def get_ssl_context():
    """
    PT-BR: Retorna contexto SSL sem verificação rígida para testes de conectividade.
    EN-US: Returns an unverified SSL context for connectivity testing.
    """
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx

def check_accessions(accessions_file):
    """
    PT-BR: Valida conectividade HTTP com os accessions do NCBI GenBank.
    EN-US: Validates HTTP connectivity for NCBI GenBank accession numbers.
    """
    print(">>> [PRE-CHECK 1/5] Verificando accessions GenBank no NCBI / Checking GenBank accessions on NCBI...")
    if not os.path.exists(accessions_file):
        print(f"[ERRO/ERROR] Arquivo de accessions não encontrado / Accessions file not found: {accessions_file}")
        return False
    
    accessions = []
    with open(accessions_file, 'r', encoding='utf-8') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[2]:
                accessions.append((parts[0], parts[1], parts[2]))
            elif len(parts) >= 1 and parts[0]:
                accessions.append((parts[0], "Unknown", parts[0]))
                
    if not accessions:
        print("[ERRO/ERROR] Nenhum accession encontrado / No accessions found.")
        return False

    success = True
    ctx = get_ssl_context()
    for strain, species, acc in accessions:
        # PT-BR: Consulta via API Entrez NCBI / EN-US: Query via NCBI Entrez API
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={acc}&retmode=json"
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                count = data.get('esearchresult', {}).get('count', '0')
                if int(count) > 0:
                    print(f"  [OK] Strain {strain} ({species}): Accession {acc} verificado no NCBI / Verified on NCBI.")
                else:
                    print(f"  [ALERTA/WARNING] Strain {strain}: Accession {acc} retornou 0 resultados / Returned 0 results.")
                    success = False
        except Exception as e:
            print(f"  [ERRO/ERROR] Falha ao conectar ao NCBI para / Connection failed for {acc}: {e}")
            success = False

    return success

def check_pfam_url():
    """
    PT-BR: Testa se o repositório Pfam no EBI FTP está acessível.
    EN-US: Tests whether the Pfam repository on EBI FTP is accessible.
    """
    print("\n>>> [PRE-CHECK 2/5] Verificando Pfam HMM (EMBL-EBI FTP)...")
    url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    ctx = get_ssl_context()
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            print(f"  [OK] Repositório Pfam-A online no EMBL-EBI (Status {resp.status}).")
            return True
    except Exception as e:
        print(f"  [AVISO/WARNING] Verificação online do Pfam retornou aviso / Pfam online check notice: {e}")
        return True

def check_bakta_db(bakta_db_path):
    """
    PT-BR: Verifica presença do banco Bakta local ou conectividade com Zenodo.
    EN-US: Verifies presence of local Bakta database or Zenodo repository connectivity.
    """
    print("\n>>> [PRE-CHECK 3/5] Verificando banco do Bakta / Checking Bakta database...")
    if bakta_db_path and os.path.exists(bakta_db_path):
        print(f"  [OK] Banco local do Bakta encontrado em / Local Bakta DB found at: {bakta_db_path}")
        return True
    
    print("  [INFO] Testando acessibilidade ao Zenodo Bakta DB / Testing Zenodo Bakta DB accessibility...")
    url = "https://doi.org/10.5281/zenodo.4247252"
    ctx = get_ssl_context()
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            print(f"  [OK] Repositório Zenodo Bakta DB online (Status {resp.status}).")
            return True
    except Exception as e:
        print(f"  [ERRO/ERROR] Falha ao acessar Bakta DB Zenodo / Failed to access Bakta DB Zenodo: {e}")
        return False

def check_references(ref_fasta):
    """
    PT-BR: Valida existência do FASTA de referências curadas.
    EN-US: Validates presence of curated references FASTA file.
    """
    print("\n>>> [PRE-CHECK 4/5] Verificando referências proteicas / Checking protein references...")
    if os.path.exists(ref_fasta) and os.path.getsize(ref_fasta) > 0:
        print(f"  [OK] Arquivo de referências curadas válido / Valid reference file: {ref_fasta}")
        return True
    else:
        print(f"  [ERRO/ERROR] Arquivo de referências ausente / Reference file missing: {ref_fasta}")
        return False

def check_runtime_tools():
    """
    PT-BR: Checa a presença de ferramentas essenciais no sistema.
    EN-US: Checks presence of essential tools on the system PATH.
    """
    print("\n>>> [PRE-CHECK 5/5] Verificando ferramentas no ambiente / Checking system tools...")
    tools = ['python', 'nextflow']
    for tool in tools:
        res = os.system(f"command -v {tool} >/dev/null 2>&1" if os.name != 'nt' else f"where {tool} >nul 2>&1")
        if res == 0:
            print(f"  [OK] Ferramenta / Tool '{tool}' ok.")
        else:
            print(f"  [AVISO/WARNING] Ferramenta / Tool '{tool}' não localizada no host PATH.")
    return True

def main():
    parser = argparse.ArgumentParser(description="Resource Checker & Pre-flight Validator (PT-BR / EN-US)")
    parser.add_argument("--accessions", default="data/accessions.tsv", help="Caminho do accessions.tsv / Path to accessions.tsv")
    parser.add_argument("--references", default="data/urease_references.fasta", help="Caminho das referencias FASTA / Path to references FASTA")
    parser.add_argument("--bakta-db", default="db/db", help="Caminho do banco Bakta / Path to Bakta DB")
    args = parser.parse_args()

    print("==================================================")
    print("   PRE-FLIGHT RESOURCE & CONNECTIVITY CHECK       ")
    print("==================================================")
    
    ok1 = check_accessions(args.accessions)
    ok2 = check_pfam_url()
    ok3 = check_bakta_db(args.bakta_db)
    ok4 = check_references(args.references)
    ok5 = check_runtime_tools()

    print("==================================================")
    if ok1 and ok3 and ok4:
        print("[SUCESSO/SUCCESS] Todos os pre-checks passaram! / All pre-checks passed!")
        sys.exit(0)
    else:
        print("[FALHA/FAILURE] Falha nos pre-checks de recursos! / Resource pre-checks failed!")
        sys.exit(4)

if __name__ == "__main__":
    main()
