#!/usr/bin/env python3
"""
extract_urease_candidates.py
PT-BR: Consolidacao de evidencias em 3 camadas (Bakta, BLASTp, Pfam HMMER) por estirpe.
EN-US: 3-Layer evidence consolidation (Bakta, BLASTp, Pfam HMMER) per strain.

PT-BR: Os limiares abaixo sao os documentados no README. Antes eles nao eram aplicados:
       todo acerto de BLAST contava como evidencia, o que inflava artificialmente a
       confianca dos candidatos.
EN-US: The thresholds below are the ones documented in the README. They used not to be
       applied: every BLAST hit counted as evidence, which artificially inflated
       candidate confidence.
"""

import sys
import os
import argparse
import json
import pandas as pd

# PT-BR: Termos EC e genes alvo / EN-US: Target EC numbers and gene symbols
TARGET_EC = ["3.5.1.5", "6.3.4.6", "3.5.1.54"]
TARGET_GENES = ["ureC", "ureA", "ureB", "ureD", "ureE", "ureF", "ureG", "ureH",
                "urtA", "urtB", "urtC", "urtD", "urtE", "uc", "ah"]

# PT-BR: Limiares de evidencia / EN-US: Evidence thresholds
BLAST_MAX_EVALUE = 1e-5
BLAST_MIN_IDENTITY = 30.0
BLAST_MIN_COVERAGE = 50.0
HMMER_MIN_DOMAIN_SCORE = 20.0

# PT-BR: Colunas do -outfmt 6 usado no main.nf / EN-US: Columns of the -outfmt 6 used in main.nf
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
                 "qlen", "slen", "evalue", "bitscore"]


def parse_bakta_json(json_file):
    """
    PT-BR: Analisa saida JSON do Bakta filtrando genes e termos EC alvos (Camada A).
    EN-US: Parses Bakta JSON output filtering target genes and EC numbers (Layer A).
    """
    candidates = []
    if not os.path.exists(json_file):
        return candidates

    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    features = data.get('features', [])
    for feat in features:
        gene = feat.get('gene') or ''
        product = feat.get('product') or ''
        dbxrefs = feat.get('db_xrefs') or feat.get('dbxrefs') or []
        # PT-BR: O Bakta usa 'locus' no cabecalho do .faa; guardamos os dois para casar
        #        com os IDs vindos do BLASTp e do HMMER.
        # EN-US: Bakta uses 'locus' in the .faa header; we keep both so they match the
        #        IDs coming from BLASTp and HMMER.
        feat_id = feat.get('locus') or feat.get('id') or ''
        alt_id = feat.get('id') or ''

        is_match = False
        if any(g.lower() == str(gene).lower() for g in TARGET_GENES):
            is_match = True
        elif any(ec in str(dbxrefs) for ec in TARGET_EC):
            is_match = True
        elif any(term in product.lower() for term in ("urease", "urea", "allophanate")):
            is_match = True

        if is_match:
            candidates.append({
                'feature_id': feat_id,
                'alt_id': alt_id,
                'gene': gene or 'unnamed',
                'product': product,
                'bakta_hit': True
            })
    return candidates


def parse_blast_tsv(blast_file):
    """
    PT-BR: Extrai IDs de consulta do BLASTp que passam nos limiares (Camada B).
    EN-US: Extracts BLASTp query IDs that pass the thresholds (Layer B).
    """
    hits = {}
    if not (os.path.exists(blast_file) and os.path.getsize(blast_file) > 0):
        return hits

    df = pd.read_csv(blast_file, sep='\t', header=None, names=BLAST_COLUMNS)
    if df.empty:
        return hits

    # PT-BR: Cobertura da consulta coberta pelo alinhamento / EN-US: Query coverage by the alignment
    df['coverage'] = (df['qend'] - df['qstart'] + 1) / df['qlen'] * 100.0

    kept = df[(df['evalue'] <= BLAST_MAX_EVALUE) &
              (df['pident'] >= BLAST_MIN_IDENTITY) &
              (df['coverage'] >= BLAST_MIN_COVERAGE)]

    dropped = len(df) - len(kept)
    if dropped:
        print(f"  [FILTRO/FILTER] BLASTp: {dropped} acerto(s) descartado(s) pelos limiares "
              f"(E<={BLAST_MAX_EVALUE}, id>={BLAST_MIN_IDENTITY}%, cov>={BLAST_MIN_COVERAGE}%) "
              f"/ hit(s) dropped by thresholds.")

    # PT-BR: Guarda o melhor acerto por consulta / EN-US: Keep the best hit per query
    for _, row in kept.sort_values('bitscore', ascending=False).iterrows():
        qid = str(row['qseqid'])
        if qid not in hits:
            hits[qid] = {
                'subject': str(row['sseqid']),
                'identity': float(row['pident']),
                'coverage': float(row['coverage']),
                'evalue': float(row['evalue']),
            }
    return hits


def parse_hmmer_tbl(hmmer_file):
    """
    PT-BR: Extrai IDs de consulta do HMMER3 domtblout acima do score de dominio (Camada C).
    EN-US: Extracts HMMER3 domtblout query IDs above the domain score (Layer C).

    PT-BR: No domtblout do hmmscan a consulta (proteina) e a coluna 4 e o score do
           dominio e a coluna 14. O codigo anterior lia a coluna 1, que e o nome do
           perfil Pfam — por isso a Camada C nunca casava com nada.
    EN-US: In hmmscan's domtblout the query (protein) is column 4 and the domain score is
           column 14. The previous code read column 1, the Pfam profile name — which is
           why Layer C never matched anything.
    """
    hits = {}
    if not (os.path.exists(hmmer_file) and os.path.getsize(hmmer_file) > 0):
        return hits

    below = 0
    with open(hmmer_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            profile = parts[0]
            query = parts[3]
            try:
                domain_score = float(parts[13])
            except ValueError:
                continue

            if domain_score < HMMER_MIN_DOMAIN_SCORE:
                below += 1
                continue

            if query not in hits or domain_score > hits[query]['score']:
                hits[query] = {'profile': profile, 'score': domain_score}

    if below:
        print(f"  [FILTRO/FILTER] HMMER: {below} dominio(s) abaixo do score "
              f"{HMMER_MIN_DOMAIN_SCORE} descartado(s) / domain(s) below score dropped.")
    return hits


def main():
    parser = argparse.ArgumentParser(description="Extract and Consolidate Urease Candidates (PT-BR / EN-US)")
    parser.add_argument("--strain", required=True, help="Identificador da estirpe / Strain ID")
    parser.add_argument("--bakta-json", required=True, help="Bakta JSON output")
    parser.add_argument("--blast-tsv", required=True, help="BLASTp TSV output")
    parser.add_argument("--hmmer-tbl", required=True, help="HMMER domtblout")
    parser.add_argument("--out-tsv", required=True, help="Caminho do TSV de saida / Output TSV path")
    args = parser.parse_args()

    bakta_cands = parse_bakta_json(args.bakta_json)
    blast_hits = parse_blast_tsv(args.blast_tsv)
    hmmer_hits = parse_hmmer_tbl(args.hmmer_tbl)

    # PT-BR: Indexa os candidatos do Bakta pelos dois identificadores possiveis.
    # EN-US: Index Bakta candidates by both possible identifiers.
    by_id = {}
    for cand in bakta_cands:
        for key in (cand['feature_id'], cand['alt_id']):
            if key:
                by_id[key] = cand

    all_feat_ids = set(by_id.keys()) | set(blast_hits.keys()) | set(hmmer_hits.keys())

    results = []
    for fid in sorted(all_feat_ids):
        b_info = by_id.get(fid)
        in_bakta = b_info is not None
        blast = blast_hits.get(fid)
        hmmer = hmmer_hits.get(fid)

        # PT-BR: Regra de Consenso: Alta (>=2 camadas), Media (1 camada)
        # EN-US: Consensus Rule: High (>=2 layers), Medium (1 layer)
        score = sum([in_bakta, blast is not None, hmmer is not None])
        if score >= 2:
            confidence = "Alta/High"
        elif score == 1:
            confidence = "Media/Medium"
        else:
            continue  # PT-BR: sem evidencia nenhuma / EN-US: no evidence at all

        results.append({
            'strain': args.strain,
            'feature_id': fid,
            'gene': b_info['gene'] if b_info else '',
            'product': b_info['product'] if b_info else '',
            'layer_bakta': in_bakta,
            'layer_blast': blast is not None,
            'layer_hmmer': hmmer is not None,
            'blast_subject': blast['subject'] if blast else '',
            'blast_identity': round(blast['identity'], 1) if blast else '',
            'blast_coverage': round(blast['coverage'], 1) if blast else '',
            'blast_evalue': blast['evalue'] if blast else '',
            'hmmer_profile': hmmer['profile'] if hmmer else '',
            'hmmer_score': hmmer['score'] if hmmer else '',
            'confidence': confidence,
        })

    columns = ['strain', 'feature_id', 'gene', 'product', 'layer_bakta', 'layer_blast',
               'layer_hmmer', 'blast_subject', 'blast_identity', 'blast_coverage',
               'blast_evalue', 'hmmer_profile', 'hmmer_score', 'confidence']

    df_out = pd.DataFrame(results, columns=columns)
    df_out.to_csv(args.out_tsv, sep='\t', index=False)

    high = sum(1 for r in results if r['confidence'] == "Alta/High")
    print(f"[OK] Candidatos consolidados salvos em / Consolidated candidates saved to: {args.out_tsv} "
          f"({len(df_out)} candidatos, {high} de alta confianca / candidates, {high} high confidence)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
