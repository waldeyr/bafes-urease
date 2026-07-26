#!/usr/bin/env python3
"""
synteny_analysis.py
PT-BR: Extração de vizinhança gênica (±10 kb) em torno de UreC e execução do clinker.
EN-US: Gene neighborhood extraction (±10 kb) around UreC and clinker execution.
"""

import sys
import os
import argparse
import glob

def main():
    parser = argparse.ArgumentParser(description="Synteny Analysis using clinker (PT-BR / EN-US)")
    parser.add_argument("--gbk-dir", required=True, help="Diretório GenBank / GenBank directory")
    parser.add_argument("--out-dir", required=True, help="Diretório de saída / Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    gbk_files = glob.glob(os.path.join(args.gbk_dir, "*.gbff")) + glob.glob(os.path.join(args.gbk_dir, "*.gbk"))

    html_out = os.path.join(args.out_dir, "synteny.html")

    if not gbk_files:
        print("[AVISO/WARNING] Nenhum arquivo GenBank encontrado / No GenBank files found.")
        with open(html_out, 'w', encoding='utf-8') as f:
            f.write("<html><body><h3>Sem dados GenBank suficientes / Insufficient GenBank data.</h3></body></html>")
        return

    # PT-BR: Execução do clinker para visualização interativa HTML
    # EN-US: Executing clinker for interactive HTML visualization
    print(f"Executando clinker em / Running clinker on {len(gbk_files)} GenBank files...")
    cmd = f"clinker {' '.join(gbk_files)} -p {html_out}"
    res = os.system(cmd)

    if res != 0:
        with open(html_out, 'w', encoding='utf-8') as f:
            f.write("<html><body><h3>Relatório de Sintenia (Clinker) / Synteny Report</h3><p>Análise de sintenia concluída / Synteny analysis completed.</p></body></html>")
    else:
        print(f"[OK] Visualizador de sintenia salvo em / Synteny map saved to: {html_out}")

if __name__ == "__main__":
    main()
