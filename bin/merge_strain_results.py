#!/usr/bin/env python3
"""
merge_strain_results.py
PT-BR: Agregação dos resultados de todas as estirpes em matriz de presença/ausência e tabela unificada.
EN-US: Aggregation of all strain candidate results into a presence/absence matrix and unified summary table.
"""

import sys
import os
import argparse
import glob
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Merge Strain Results into Presence/Absence Matrix (PT-BR / EN-US)")
    parser.add_argument("--candidates-dir", required=True, help="Diretório de candidatos / Candidates directory")
    parser.add_argument("--out-matrix", required=True, help="Saída urease_presence_absence.tsv / Output matrix file")
    parser.add_argument("--out-summary", required=True, help="Saída summary_unified.tsv / Output summary file")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.candidates_dir, "*.tsv"))
    if not files:
        print("[AVISO/WARNING] Nenhum arquivo TSV encontrado / No candidate TSV files found.")
        pd.DataFrame(columns=['strain', 'ureC', 'ureA', 'ureB', 'ureD', 'ureE', 'ureF', 'ureG', 'urtA', 'uc', 'ah']).to_csv(args.out_matrix, sep='\t', index=False)
        pd.DataFrame().to_csv(args.out_summary, sep='\t', index=False)
        return

    all_dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty:
                all_dfs.append(df)
        except Exception as e:
            print(f"[ERRO/ERROR] Falha ao ler / Failed reading {f}: {e}")

    if not all_dfs:
        pd.DataFrame(columns=['strain', 'ureC', 'ureA', 'ureB', 'ureD', 'ureE', 'ureF', 'ureG', 'urtA', 'uc', 'ah']).to_csv(args.out_matrix, sep='\t', index=False)
        pd.DataFrame().to_csv(args.out_summary, sep='\t', index=False)
        return

    merged_df = pd.concat(all_dfs, ignore_index=True)
    merged_df.to_csv(args.out_summary, sep='\t', index=False)

    # PT-BR: Matriz binária de Presença (1) / Ausência (0) por estirpe
    # EN-US: Binary Presence (1) / Absence (0) matrix per strain
    strains = merged_df['strain'].unique()
    target_genes = ['ureC', 'ureA', 'ureB', 'ureD', 'ureE', 'ureF', 'ureG', 'urtA', 'uc', 'ah']
    
    matrix_rows = []
    for s in strains:
        sub = merged_df[merged_df['strain'] == s]
        row = {'strain': s}
        for g in target_genes:
            has_gene = any(g.lower() in str(gene).lower() for gene in sub['gene']) or any(g.lower() in str(prod).lower() for prod in sub['product'])
            row[g] = "1" if has_gene else "0"
        matrix_rows.append(row)

    matrix_df = pd.DataFrame(matrix_rows)
    matrix_df.to_csv(args.out_matrix, sep='\t', index=False)
    print(f"[OK] Matriz salva em / Matrix saved to: {args.out_matrix}")
    print(f"[OK] Sumário salvo em / Summary saved to: {args.out_summary}")

if __name__ == "__main__":
    main()
