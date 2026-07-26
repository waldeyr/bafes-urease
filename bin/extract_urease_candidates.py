#!/usr/bin/env python3
"""
extract_urease_candidates.py
PT-BR: Consolidação de evidências em 3 camadas (Bakta, BLASTp, Pfam HMMER) por estirpe.
EN-US: 3-Layer evidence consolidation (Bakta, BLASTp, Pfam HMMER) per strain.
"""

import sys
import os
import argparse
import json
import pandas as pd

# PT-BR: Termos EC e genes alvo / EN-US: Target EC numbers and gene symbols
TARGET_EC = ["3.5.1.5", "6.3.4.6", "3.5.1.54"]
TARGET_GENES = ["ureC", "ureA", "ureB", "ureD", "ureE", "ureF", "ureG", "urtA", "urtB", "urtC", "urtD", "urtE", "uc", "ah"]

def parse_bakta_json(json_file):
    """
    PT-BR: Analisa saída JSON do Bakta filtrando genes e termos EC alvos (Camada A).
    EN-US: Parses Bakta JSON output filtering target genes and EC numbers (Layer A).
    """
    candidates = []
    if not os.path.exists(json_file):
        return candidates
    
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
        
    features = data.get('features', [])
    for feat in features:
        gene = feat.get('gene', '')
        product = feat.get('product', '')
        dbxrefs = feat.get('dbxrefs', [])
        feat_id = feat.get('id', '')
        
        is_match = False
        if any(g.lower() in str(gene).lower() for g in TARGET_GENES):
            is_match = True
        elif any(ec in str(dbxrefs) for ec in TARGET_EC):
            is_match = True
        elif "urease" in product.lower() or "urea" in product.lower() or "allophanate" in product.lower():
            is_match = True
            
        if is_match:
            candidates.append({
                'feature_id': feat_id,
                'gene': gene or 'unnamed',
                'product': product,
                'bakta_hit': True
            })
    return candidates

def parse_blast_tsv(blast_file):
    """
    PT-BR: Extrai IDs de acertos do BLASTp (Camada B).
    EN-US: Extracts query IDs from BLASTp results (Layer B).
    """
    hits = set()
    if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
        df = pd.read_csv(blast_file, sep='\t', header=None)
        if not df.empty:
            hits = set(df[0].astype(str))
    return hits

def parse_hmmer_tbl(hmmer_file):
    """
    PT-BR: Extrai IDs de acertos do HMMER3 domtblout (Camada C).
    EN-US: Extracts query IDs from HMMER3 domtblout results (Layer C).
    """
    hits = set()
    if os.path.exists(hmmer_file) and os.path.getsize(hmmer_file) > 0:
        with open(hmmer_file, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split()
                    if parts:
                        hits.add(parts[0])
    return hits

def main():
    parser = argparse.ArgumentParser(description="Extract and Consolidate Urease Candidates (PT-BR / EN-US)")
    parser.add_argument("--strain", required=True, help="Identificador da estirpe / Strain ID")
    parser.add_argument("--bakta-json", required=True, help="Bakta JSON output")
    parser.add_argument("--blast-tsv", required=True, help="BLASTp TSV output")
    parser.add_argument("--hmmer-tbl", required=True, help="HMMER domtblout")
    parser.add_argument("--out-tsv", required=True, help="Caminho do TSV de saída / Output TSV path")
    args = parser.parse_args()

    bakta_cands = parse_bakta_json(args.bakta_json)
    blast_hits = parse_blast_tsv(args.blast_tsv)
    hmmer_hits = parse_hmmer_tbl(args.hmmer_tbl)

    all_feat_ids = set([c['feature_id'] for c in bakta_cands]).union(blast_hits).union(hmmer_hits)

    results = []
    for fid in all_feat_ids:
        b_info = next((c for c in bakta_cands if c['feature_id'] == fid), None)
        in_bakta = b_info is not None
        in_blast = fid in blast_hits
        in_hmmer = fid in hmmer_hits

        # PT-BR: Regra de Consenso: Alta (>=2 camadas), Media (1 camada)
        # EN-US: Consensus Rule: High (>=2 layers), Medium (1 layer)
        score = sum([in_bakta, in_blast, in_hmmer])
        confidence = "Alta/High" if score >= 2 else "Media/Medium" if score == 1 else "Baixa/Low"

        gene_symbol = b_info['gene'] if b_info else 'candidate'
        product_desc = b_info['product'] if b_info else 'Predicted Urease-related protein'

        results.append({
            'strain': args.strain,
            'feature_id': fid,
            'gene': gene_symbol,
            'product': product_desc,
            'layer_bakta': in_bakta,
            'layer_blast': in_blast,
            'layer_hmmer': in_hmmer,
            'confidence': confidence
        })

    df_out = pd.DataFrame(results)
    if df_out.empty:
        df_out = pd.DataFrame(columns=['strain', 'feature_id', 'gene', 'product', 'layer_bakta', 'layer_blast', 'layer_hmmer', 'confidence'])

    df_out.to_csv(args.out_tsv, sep='\t', index=False)
    print(f"[OK] Candidatos consolidados salvos em / Consolidated candidates saved to: {args.out_tsv} ({len(df_out)} cands)")

if __name__ == "__main__":
    main()
