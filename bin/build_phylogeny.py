#!/usr/bin/env python3
"""
build_phylogeny.py
PT-BR: Seleciona as sequencias de UreC, alinha (MAFFT), apara (trimAl) e infere a arvore (IQ-TREE 2).
EN-US: Selects UreC sequences, aligns (MAFFT), trims (trimAl), and infers the tree (IQ-TREE 2).

PT-BR: Duas garantias — (1) a arvore e construida apenas com as sequencias de UreC
       identificadas na triagem, nunca com o proteoma inteiro; (2) falha de ferramenta
       aborta com erro, em vez de deixar para tras um arquivo parcial que passaria por
       resultado. Quando ha menos de 3 UreC, isso e registrado como resultado legitimo
       (nao ha arvore possivel) e nenhum arquivo de arvore e escrito.
EN-US: Two guarantees — (1) the tree is built only from the UreC sequences identified in
       the screening, never from the whole proteome; (2) a tool failure aborts with an
       error instead of leaving behind a partial file that would pass for a result. With
       fewer than 3 UreC sequences that is recorded as a legitimate outcome (no tree is
       possible) and no tree file is written.
"""

import sys
import os
import glob
import argparse
import subprocess

# PT-BR: Termos que identificam a subunidade alfa da urease / EN-US: Terms identifying the urease alpha subunit
UREC_GENE_TERMS = ("urec",)
UREC_PRODUCT_TERMS = ("urease subunit alpha", "urease alpha", "urea amidohydrolase subunit alpha")


def read_fasta(path):
    """
    PT-BR: Le um FASTA e devolve {id: (cabecalho, sequencia)} usando o 1o token como id.
    EN-US: Reads a FASTA and returns {id: (header, sequence)} using the 1st token as id.
    """
    records = {}
    header = None
    chunks = []
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    records[header.split()[0]] = (header, ''.join(chunks))
                header = line[1:]
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records[header.split()[0]] = (header, ''.join(chunks))
    return records


def collect_urec_ids(candidates_dir):
    """
    PT-BR: Varre os TSVs de candidatos e devolve os feature_id anotados como UreC.
    EN-US: Scans the candidate TSVs and returns the feature_ids annotated as UreC.
    """
    wanted = set()
    for tsv in sorted(glob.glob(os.path.join(candidates_dir, "*_candidates.tsv"))):
        with open(tsv, 'r', encoding='utf-8') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            try:
                i_id = header.index('feature_id')
                i_gene = header.index('gene')
                i_product = header.index('product')
            except ValueError:
                continue
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_id, i_gene, i_product):
                    continue
                gene = parts[i_gene].strip().lower()
                product = parts[i_product].strip().lower()
                if any(t in gene for t in UREC_GENE_TERMS) or any(t in product for t in UREC_PRODUCT_TERMS):
                    wanted.add(parts[i_id].strip())
    return wanted


def run(cmd, label):
    """
    PT-BR: Executa um comando externo e aborta se ele falhar.
    EN-US: Runs an external command and aborts if it fails.
    """
    print(f"  $ {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        print(f"[ERRO/ERROR] '{cmd[0]}' nao encontrado no PATH / not found on PATH.")
        return False
    if proc.returncode != 0:
        print(f"[ERRO/ERROR] {label} falhou / failed (exit {proc.returncode}).")
        if proc.stderr:
            print(proc.stderr.strip()[:2000])
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="UreC Phylogeny Reconstruction (PT-BR / EN-US)")
    parser.add_argument("--fasta-in", required=True, help="FASTA com todas as proteinas / FASTA with all proteins")
    parser.add_argument("--candidates-dir", required=True, help="Diretorio dos TSVs de candidatos / Candidate TSV directory")
    parser.add_argument("--out-dir", required=True, help="Diretorio de saida / Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    status_path = os.path.join(args.out_dir, "phylogeny_status.txt")
    urec_fasta = os.path.join(args.out_dir, "ureC.faa")
    aln_out = os.path.join(args.out_dir, "ureC.aln")
    trim_out = os.path.join(args.out_dir, "ureC_trimmed.aln")
    tree_prefix = os.path.join(args.out_dir, "ureC")

    def write_status(text):
        with open(status_path, 'w', encoding='utf-8') as handle:
            handle.write(text + '\n')
        print(text)

    if not os.path.exists(args.fasta_in) or os.path.getsize(args.fasta_in) == 0:
        print(f"[ERRO/ERROR] FASTA de proteinas ausente ou vazio / protein FASTA missing or empty: {args.fasta_in}")
        return 1

    proteins = read_fasta(args.fasta_in)
    urec_ids = collect_urec_ids(args.candidates_dir)
    print(f"Candidatos UreC na triagem / UreC candidates in screening: {len(urec_ids)}")

    selected = [(pid, proteins[pid]) for pid in sorted(urec_ids) if pid in proteins]
    missing = len(urec_ids) - len(selected)
    if missing:
        print(f"[AVISO/WARNING] {missing} id(s) de UreC sem sequencia correspondente no FASTA "
              f"/ UreC id(s) with no matching sequence in the FASTA.")

    if len(selected) < 3:
        write_status(
            f"SEM ARVORE / NO TREE: apenas {len(selected)} sequencia(s) de UreC recuperada(s); "
            f"IQ-TREE exige no minimo 3. / only {len(selected)} UreC sequence(s) recovered; "
            f"IQ-TREE requires at least 3. Nenhuma arvore foi gerada e nenhum dado foi fabricado. "
            f"/ No tree was produced and no data was fabricated."
        )
        return 0

    # PT-BR: IQ-TREE trunca nomes de sequencia em '|' e recusa nomes repetidos, entao os
    #        identificadores sao saneados. O mapeamento e gravado para que cada folha da
    #        arvore continue rastreavel ate o feature_id original.
    # EN-US: IQ-TREE truncates sequence names at '|' and rejects duplicates, so the
    #        identifiers are sanitised. The mapping is written so every leaf of the tree
    #        stays traceable back to its original feature_id.
    labels_path = os.path.join(args.out_dir, "ureC_labels.tsv")
    used = set()
    mapping = []

    with open(urec_fasta, 'w', encoding='utf-8') as handle:
        for pid, (header, sequence) in selected:
            safe = ''.join(ch if (ch.isalnum() or ch in '_-.') else '_' for ch in pid)
            candidate = safe
            suffix = 1
            while candidate in used:
                suffix += 1
                candidate = f"{safe}_{suffix}"
            used.add(candidate)
            mapping.append((candidate, pid, header))

            handle.write(f">{candidate}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + '\n')

    with open(labels_path, 'w', encoding='utf-8') as handle:
        handle.write("tree_label\tfeature_id\toriginal_header\n")
        for safe, pid, header in mapping:
            handle.write(f"{safe}\t{pid}\t{header}\n")

    print(f"Alinhando {len(selected)} sequencias UreC com MAFFT / Aligning {len(selected)} UreC sequences with MAFFT...")
    try:
        with open(aln_out, 'w', encoding='utf-8') as out:
            proc = subprocess.run(["mafft", "--auto", urec_fasta], stdout=out,
                                  stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        print("[ERRO/ERROR] 'mafft' nao encontrado no PATH / not found on PATH.")
        return 1
    if proc.returncode != 0 or os.path.getsize(aln_out) == 0:
        print(f"[ERRO/ERROR] MAFFT falhou / failed (exit {proc.returncode}).")
        if proc.stderr:
            print(proc.stderr.strip()[:2000])
        return 1

    if not run(["trimal", "-in", aln_out, "-out", trim_out, "-automated1"], "trimAl"):
        return 1

    if not run(["iqtree", "-s", trim_out, "--prefix", tree_prefix,
                "-m", "MFP", "-bb", "1000", "-redo"], "IQ-TREE 2"):
        return 1

    treefile = f"{tree_prefix}.treefile"
    if not os.path.exists(treefile) or os.path.getsize(treefile) == 0:
        print(f"[ERRO/ERROR] IQ-TREE terminou sem gerar / finished without producing: {treefile}")
        return 1

    write_status(
        f"ARVORE GERADA / TREE BUILT: {len(selected)} sequencias UreC, 1000 ultrafast bootstraps. "
        f"/ {len(selected)} UreC sequences, 1000 ultrafast bootstraps. Arquivo / File: {treefile}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
