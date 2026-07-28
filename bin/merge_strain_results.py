#!/usr/bin/env python3
"""
merge_strain_results.py
PT-BR: Agregação dos resultados de todas as estirpes em matriz de presença/ausência e tabela unificada.
EN-US: Aggregation of all strain candidate results into a presence/absence matrix and unified summary table.

PT-BR: A matriz comparava o símbolo do gene como substring, então 'ureA' casava com
       qualquer produto contendo 'urea': S6 era marcada com ureA=1 apenas por causa de
       "Aminoglycoside/hydroxyurea antibiotic resistance kinase". Agora o símbolo do gene
       é comparado por igualdade e o produto por padrões curados.
EN-US: The matrix compared the gene symbol as a substring, so 'ureA' matched any product
       containing 'urea': S6 was scored ureA=1 purely because of "Aminoglycoside/hydroxyurea
       antibiotic resistance kinase". The gene symbol is now compared by equality and the
       product by curated patterns.
"""

import os
import argparse
import glob
import re
import pandas as pd

# PT-BR: Genes da matriz, na ordem das colunas / EN-US: Matrix genes, in column order
TARGET_GENES = ['ureA', 'ureB', 'ureC', 'ureD', 'ureE', 'ureF', 'ureG', 'ureH',
                'urtA', 'uc', 'ah']

# PT-BR: Padrões de produto por gene. Só entram descrições inequívocas — o produto é texto
#        livre do Bakta, e um padrão frouxo aqui reintroduz exatamente o falso positivo que
#        esta correção elimina. 'uc' = urea carboxylase, 'ah' = allophanate hydrolase.
# EN-US: Per-gene product patterns. Only unambiguous descriptions qualify — the product is
#        free text from Bakta, and a loose pattern here reintroduces exactly the false
#        positive this fix removes. 'uc' = urea carboxylase, 'ah' = allophanate hydrolase.
PRODUCT_PATTERNS = {
    'ureA': [r'urease\s+subunit\s+gamma', r'urease\s+gamma\s+subunit'],
    'ureB': [r'urease\s+subunit\s+beta', r'urease\s+beta\s+subunit'],
    'ureC': [r'urease\s+subunit\s+alpha', r'urease\s+alpha\s+subunit',
             r'urea\s+amidohydrolase\s+subunit\s+alpha'],
    'ureD': [r'urease\s+accessory\s+protein\s+ured'],
    'ureE': [r'urease\s+accessory\s+protein\s+uree'],
    'ureF': [r'urease\s+accessory\s+protein\s+uref'],
    'ureG': [r'urease\s+accessory\s+protein\s+ureg'],
    'ureH': [r'urease\s+accessory\s+protein\s+ureh'],
    'urtA': [r'urea\s+(abc\s+)?transport', r'urea\s+transporter'],
    'uc': [r'urea\s+carboxylase'],
    'ah': [r'allophanate\s+hydrolase'],
}

COMPILED_PATTERNS = {gene: [re.compile(p, re.IGNORECASE) for p in pats]
                     for gene, pats in PRODUCT_PATTERNS.items()}


def strain_from_filename(path):
    """
    PT-BR: Deriva a estirpe do nome do arquivo <strain>_candidates.tsv.
    EN-US: Derives the strain from the <strain>_candidates.tsv filename.
    """
    base = os.path.basename(path)
    return base[:-len('_candidates.tsv')] if base.endswith('_candidates.tsv') else os.path.splitext(base)[0]


def has_gene(sub_df, gene):
    """
    PT-BR: Presença = símbolo do gene igual (sem distinção de caixa) OU produto casando com
           um dos padrões curados. Nunca substring do símbolo dentro do produto.
    EN-US: Presence = gene symbol equal (case-insensitive) OR product matching one of the
           curated patterns. Never a substring of the symbol inside the product.
    """
    for value in sub_df['gene']:
        if str(value).strip().lower() == gene.lower():
            return True
    for pattern in COMPILED_PATTERNS.get(gene, []):
        for value in sub_df['product']:
            if pattern.search(str(value)):
                return True
    return False


def main():
    parser = argparse.ArgumentParser(description="Merge Strain Results into Presence/Absence Matrix (PT-BR / EN-US)")
    parser.add_argument("--candidates-dir", required=True, help="Diretório de candidatos / Candidates directory")
    parser.add_argument("--out-matrix", required=True, help="Saída urease_presence_absence.tsv / Output matrix file")
    parser.add_argument("--out-summary", required=True, help="Saída summary_unified.tsv / Output summary file")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.candidates_dir, "*.tsv")))
    if not files:
        print("[AVISO/WARNING] Nenhum arquivo TSV encontrado / No candidate TSV files found.")
        pd.DataFrame(columns=['strain'] + TARGET_GENES).to_csv(args.out_matrix, sep='\t', index=False)
        pd.DataFrame().to_csv(args.out_summary, sep='\t', index=False)
        return

    # PT-BR: Toda estirpe com um TSV entra na matriz, mesmo sem candidato nenhum — ausência
    #        de urease é um resultado, e uma estirpe sumir da matriz seria indistinguível de
    #        uma estirpe que nunca foi processada.
    # EN-US: Every strain with a TSV enters the matrix, even with no candidate at all —
    #        absence of urease is a result, and a strain vanishing from the matrix would be
    #        indistinguishable from a strain that was never processed.
    all_strains = []
    all_dfs = []
    for f in files:
        strain = strain_from_filename(f)
        if strain not in all_strains:
            all_strains.append(strain)
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty:
                all_dfs.append(df)
        except Exception as e:
            print(f"[ERRO/ERROR] Falha ao ler / Failed reading {f}: {e}")

    if all_dfs:
        merged_df = pd.concat(all_dfs, ignore_index=True)
    else:
        merged_df = pd.DataFrame(columns=['strain', 'gene', 'product'])
    merged_df.to_csv(args.out_summary, sep='\t', index=False)

    # PT-BR: Matriz binária de Presença (1) / Ausência (0) por estirpe
    # EN-US: Binary Presence (1) / Absence (0) matrix per strain
    matrix_rows = []
    for s in all_strains:
        sub = merged_df[merged_df['strain'] == s] if not merged_df.empty else merged_df
        row = {'strain': s}
        for g in TARGET_GENES:
            row[g] = "1" if (not sub.empty and has_gene(sub, g)) else "0"
        matrix_rows.append(row)

    matrix_df = pd.DataFrame(matrix_rows, columns=['strain'] + TARGET_GENES)
    matrix_df.to_csv(args.out_matrix, sep='\t', index=False)
    print(f"[OK] Matriz salva em / Matrix saved to: {args.out_matrix} "
          f"({len(matrix_rows)} estirpes / strains)")
    print(f"[OK] Sumário salvo em / Summary saved to: {args.out_summary} "
          f"({len(merged_df)} candidatos / candidates)")


if __name__ == "__main__":
    main()
