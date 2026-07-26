#!/usr/bin/env python3
"""
build_phylogeny.py
PT-BR: Alinhamento de sequências UreC (MAFFT), trimming (trimAl) e reconstrução filogenética (IQ-TREE 2).
EN-US: UreC sequence alignment (MAFFT), trimming (trimAl), and phylogenetic reconstruction (IQ-TREE 2).
"""

import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="UreC Phylogeny Reconstruction (PT-BR / EN-US)")
    parser.add_argument("--fasta-in", required=True, help="FASTA de entrada / Input FASTA")
    parser.add_argument("--out-dir", required=True, help="Diretório de saída / Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    aln_out = os.path.join(args.out_dir, "ureC.aln")
    trim_out = os.path.join(args.out_dir, "ureC_trimmed.aln")
    tree_out_prefix = os.path.join(args.out_dir, "ureC")

    if not os.path.exists(args.fasta_in) or os.path.getsize(args.fasta_in) == 0:
        print("[AVISO/WARNING] Arquivo FASTA de UreC ausente / UreC FASTA file missing.")
        return

    # PT-BR: Contar sequências FASTA / EN-US: Count FASTA sequences
    count = 0
    with open(args.fasta_in, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1

    # PT-BR: Mínimo de 3 sequências para filogenia / EN-US: Minimum 3 sequences required for tree
    if count < 3:
        print(f"[AVISO/WARNING] Menos de 3 sequências UreC ({count} encontradas) / Less than 3 UreC sequences found.")
        return

    # PT-BR: Alinhamento, Trimming e IQ-TREE 2 / EN-US: Alignment, Trimming, and IQ-TREE 2
    print(f"Alinhando {count} sequências UreC com MAFFT...")
    os.system(f"mafft --auto {args.fasta_in} > {aln_out}")

    print("Limpando alinhamento com trimAl...")
    os.system(f"trimal -in {aln_out} -out {trim_out} -automated1")

    print("Inferindo árvore filogenética com IQ-TREE 2...")
    os.system(f"iqtree -s {trim_out} --prefix {tree_out_prefix} -m AUTO -bb 1000 -redo")

    print(f"[OK] Reconstrução filogenética concluída / Phylogeny completed: {tree_out_prefix}.treefile")

if __name__ == "__main__":
    main()
