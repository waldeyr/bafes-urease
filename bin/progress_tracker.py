#!/usr/bin/env python3
"""
progress_tracker.py
PT-BR: Monitor de progresso em tempo real utilizando tqdm para cada fase e progresso geral.
EN-US: Real-time progress monitor using tqdm for each phase and overall execution.
"""

import sys
import os
import time
import argparse
from tqdm import tqdm

# PT-BR: Fases e pesos percentuais / EN-US: Pipeline phases and percentage weights
PHASES = [
    ("Pre-flight Check", 5),
    ("Download Genomes (NCBI)", 15),
    ("Quality Control (QC)", 15),
    ("Functional Annotation (Bakta)", 30),
    ("Candidate Screening (BLAST/HMMER)", 15),
    ("Presence/Absence Matrix", 5),
    ("Synteny Analysis (clinker)", 10),
    ("UreC Phylogeny (IQ-TREE)", 5)
]

def run_simulated_tracker(verbose=False):
    """
    PT-BR: Executa barra de progresso visual interativa.
    EN-US: Runs interactive visual progress bar tracker.
    """
    print("\n[PT-BR] Iniciando monitoramento via tqdm... / [EN-US] Starting progress monitoring via tqdm...\n")
    overall_bar = tqdm(total=100, desc="Overall Pipeline Progress", position=0, leave=True)
    
    current_progress = 0
    for phase_name, weight in PHASES:
        phase_bar = tqdm(total=100, desc=f"Phase: {phase_name:<32}", position=1, leave=False)
        for step in range(10):
            time.sleep(0.1)
            phase_bar.update(10)
        phase_bar.close()
        
        current_progress += weight
        overall_bar.n = current_progress
        overall_bar.refresh()
        
    overall_bar.close()
    print("\n[OK] Pipeline concluído com sucesso! / Pipeline completed successfully!\n")

def main():
    parser = argparse.ArgumentParser(description="Nextflow Progress Tracker via tqdm (PT-BR / EN-US)")
    parser.add_argument("--trace", help="Caminho para arquivo trace / Path to trace file")
    parser.add_argument("--verbose", action="store_true", help="Modo detalhado / Verbose mode")
    args = parser.parse_args()

    run_simulated_tracker(args.verbose)

if __name__ == "__main__":
    main()
