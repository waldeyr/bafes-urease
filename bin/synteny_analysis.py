#!/usr/bin/env python3
"""
synteny_analysis.py
PT-BR: Execucao do clinker sobre os GenBank anotados para gerar o mapa de sintenia.
EN-US: Runs clinker over the annotated GenBank files to produce the synteny map.

PT-BR: Este script NUNCA escreve um HTML de placeholder. Se nao houver GenBank ou se o
       clinker falhar, ele sai com erro — um relatorio vazio seria indistinguivel de um
       resultado real para quem le a saida depois.
EN-US: This script NEVER writes a placeholder HTML. If there are no GenBank files or if
       clinker fails, it exits with an error — an empty report would be indistinguishable
       from a real result for whoever reads the output later.
"""

import sys
import os
import argparse
import glob
import subprocess


def main():
    parser = argparse.ArgumentParser(description="Synteny Analysis using clinker (PT-BR / EN-US)")
    parser.add_argument("--gbk-dir", required=True, help="Diretorio GenBank / GenBank directory")
    parser.add_argument("--out-dir", required=True, help="Diretorio de saida / Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    gbk_files = sorted(glob.glob(os.path.join(args.gbk_dir, "*.gbff")) +
                       glob.glob(os.path.join(args.gbk_dir, "*.gbk")))

    html_out = os.path.join(args.out_dir, "synteny.html")

    if not gbk_files:
        print(f"[ERRO/ERROR] Nenhum arquivo GenBank em / No GenBank file in: {args.gbk_dir}")
        print("             A anotacao do Bakta precisa ter sido concluida antes desta etapa.")
        print("             Bakta annotation must complete before this step.")
        return 1

    if len(gbk_files) < 2:
        print(f"[ERRO/ERROR] Sintenia exige ao menos 2 genomas; encontrado(s) {len(gbk_files)}.")
        print(f"             Synteny requires at least 2 genomes; found {len(gbk_files)}.")
        return 1

    print(f"Executando clinker em / Running clinker on {len(gbk_files)} GenBank files...")
    cmd = ["clinker", *gbk_files, "-p", html_out]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        print("[ERRO/ERROR] 'clinker' nao encontrado no PATH / not found on PATH.")
        return 1

    if proc.returncode != 0:
        print(f"[ERRO/ERROR] clinker falhou / failed (exit {proc.returncode}).")
        if proc.stderr:
            print(proc.stderr.strip()[:2000])
        return 1

    if not os.path.exists(html_out) or os.path.getsize(html_out) == 0:
        print(f"[ERRO/ERROR] clinker terminou sem gerar / finished without producing: {html_out}")
        return 1

    print(f"[OK] Visualizador de sintenia salvo em / Synteny map saved to: {html_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
